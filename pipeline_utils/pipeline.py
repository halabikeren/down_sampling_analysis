import shutil
import typing as t
from pipeline_utils import (
    SequenceDataType,
    AlignmentMethod,
    TreeReconstructionMethod,
    PipelineInput,
)
from programs import *
import re
from dataclasses import dataclass
import time
from Bio.Phylo import NewickIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from ete3 import Tree
from copy import deepcopy
import subprocess
from Bio import SeqIO, AlignIO
import os
from samplers import *

import logging

logger = logging.getLogger(__name__)


@dataclass
class Pipeline:
    pipeline_dir: str
    unaligned_sequence_data_path: str
    unaligned_sequence_data: t.List[SeqIO.SeqRecord]
    aligned_sequence_data_path: str
    tree_path: str
    samples_info: t.Dict[
        float, t.Dict[SamplingMethod, t.Dict[str, t.Any]]
    ]  # will map sample fractions to maps of sampling methods to their info (input path, output paths,
    # analysis results ect.)

    def __init__(self, pipeline_input: PipelineInput):

        # set initial parameters
        dataset_name_regex = re.compile("(.*)\.")
        dataset_name = dataset_name_regex.search(
            os.path.basename(pipeline_input.unaligned_sequence_data_path)
        ).group(1)
        self.pipeline_dir = pipeline_input.pipeline_dir
        processed_data_dir = f"{self.pipeline_dir}/input_data/"
        os.makedirs(processed_data_dir)
        logger.info(f"Setting input for pipeline at {processed_data_dir}")

        self.unaligned_sequence_data_path = (
            f"{processed_data_dir}{dataset_name}_unaligned.fasta"
        )
        shutil.copyfile(
            pipeline_input.unaligned_sequence_data_path,
            self.unaligned_sequence_data_path,
        )
        logger.info(
            f"Unaligned sequence data saved at {self.unaligned_sequence_data_path}"
        )

        self.aligned_sequence_data_path = (
            f"{processed_data_dir}{dataset_name}_aligned.fasta"
        )
        if pipeline_input.aligned_sequence_data_path:
            shutil.copyfile(
                pipeline_input.aligned_sequence_data_path,
                self.aligned_sequence_data_path,
            )
            logger.info(f"Aligned data saved at {self.aligned_sequence_data_path}")

        self.tree_path = f"{processed_data_dir}{dataset_name}_tree.nwk"
        if pipeline_input.tree_path:
            shutil.copyfile(pipeline_input.tree_path, self.tree_path)
            logger.info(f"Tree saved at {self.tree_path}")

        # fill in available parameters
        self.sequence_data_type = pipeline_input.sequence_data_type

        # align, if needed
        if not os.path.exists(self.aligned_sequence_data_path):
            self.aligned_sequence_data_path = Pipeline.align(
                self.unaligned_sequence_data_path,
                self.aligned_sequence_data_path,
                pipeline_input.sequence_data_type,
                pipeline_input.alignment_method,
                pipeline_input.alignment_params,
            )
            logger.info(
                f"Alignment generated successfully using {pipeline_input.alignment_method.value}"
            )

        self.unaligned_sequence_data = list(
            SeqIO.parse(self.unaligned_sequence_data_path, "fasta")
        )
        logger.info(
            f"Processed sequence data of size {len(self.unaligned_sequence_data)}"
        )

        # build tree, if needed
        if not os.path.exists(self.tree_path):
            self.build_tree(
                self.aligned_sequence_data_path,
                self.tree_path,
                pipeline_input.tree_reconstruction_method,
                pipeline_input.tree_reconstruction_params,
            )
            logger.info(
                f"Tree reconstructed successfully at {self.tree_path} using {pipeline_input.tree_reconstruction_method.value}"
            )

        # set sampling info structure
        self.samples_info = dict()
        for fraction in pipeline_input.sampling_fractions:
            self.samples_info[fraction] = dict()
            for method in pipeline_input.sampling_methods:
                self.samples_info[fraction][method.value] = {
                    "unaligned_sequence_data_path": None,
                    "aligned_sequence_data_path": None,
                    "programs_performance": dict(),
                }
                for program_name in pipeline_input.programs:
                    self.samples_info[fraction][method.value]["programs_performance"][
                        program_name.value
                    ] = {"input_path": None, "output_path": None, "result": None}

    @staticmethod
    def translate(sequence_records: t.List[SeqIO.SeqRecord], output_path: str):
        """
        :param sequence_records: list of coding sequence records
        :param output_path: the output path to which the translated sequences corresponding to the coding sequence
        records will be written to in fasta format
        :return: None
        """
        aa_records = deepcopy(sequence_records)
        for aa_record in aa_records:
            aa_record.seq = aa_record.seq.translate(table=1)
        SeqIO.write(aa_records, output_path, "fasta")

    @staticmethod
    def reverse_translate(
        unaligned_codon_path: str, aligned_aa_path: str, aligned_codon_path: str
    ):
        """
        :param unaligned_codon_path: path to non-aligned coding sequences
        :param aligned_aa_path: path to the aligned corresponding translated sequences
        :param aligned_codon_path: path to write the aligned coding sequences to
        :return: None
        """
        unaligned_codon_records = list(SeqIO.parse(unaligned_codon_path, "fasta"))
        aligned_aa_records = list(SeqIO.parse(aligned_aa_path, "fasta"))
        aligned_codon_records = []
        for aligned_aa_record in aligned_aa_records:
            unaligned_codon_record = [
                sequence_record
                for sequence_record in unaligned_codon_records
                if sequence_record.description == aligned_aa_record.description
            ][0]
            aligned_codon_record = deepcopy(unaligned_codon_record)
            aa_index = 0
            codon_index = 0
            while aa_index < len(aligned_aa_record.seq):
                if aligned_aa_record.seq[aa_index] != "-":
                    aligned_codon_record.seq += unaligned_codon_record.seq[
                        codon_index * 3 : codon_index * 3 + 3
                    ]
                    codon_index += 1
                else:
                    aligned_codon_record.seq += "---"
                aa_index += 1
            aligned_codon_records.append(aligned_codon_record)
        SeqIO.write(aligned_codon_records, aligned_codon_path, "fasta")

    @staticmethod
    def align(
        input_path: str,
        output_path: str,
        sequence_data_type: SequenceDataType,
        alignment_method: AlignmentMethod,
        alignment_params: t.Optional[t.Dict[str, t.Any]] = None,
    ) -> str:
        """
        Aligns sequence data according to given method
        input_path: unaligned sequence data path
        output_path: the path to write the alignment to
        alignment_method: the method that the sequence data should be aligned with
        alignment_params: the parameters that the method should be run with
        :return:
        """
        if os.path.exists(output_path):
            alternative_alignment_path = output_path.replace(
                ".fasta", f"_{time.time()}.fasta"
            )
            logger.warning(
                f"{output_path} already exists. The program will create the alignment in the "
                f"alternative path {alternative_alignment_path} "
            )
            output_path = alternative_alignment_path

        sequence_records = list(SeqIO.parse(input_path, "fasta"))
        alignment_input_path = input_path
        alignment_output_path = output_path
        if (
            sequence_data_type == SequenceDataType.CODON
            and alignment_method == AlignmentMethod.MAFFT
        ):
            alignment_input_path = output_path.replace(".fasta", "_translated.fasta")
            Pipeline.translate(
                sequence_records,
                input_path.replace(".fasta", "_translated.fasta"),
            )
            alignment_output_path = output_path.replace(".fasta", "_translated.fasta")

        cmd = ""
        if not alignment_method or alignment_method == AlignmentMethod.MAFFT:
            cmd = f"mafft --localpair --maxiterate 1000 {alignment_input_path} > {alignment_output_path}"
        elif alignment_method == AlignmentMethod.PRANK:
            cmd = f"prank -d={alignment_input_path} -o={alignment_output_path} -f=fasta -support {'-codon' if sequence_data_type == SequenceDataType.CODON else ''} -iterate=100 -showtree "
        if alignment_params:
            cmd += " ".join(
                [
                    f"{param_name} {alignment_params[param_name]}"
                    for param_name in alignment_params
                ]
            )
        process = subprocess.run(cmd, shell=True, capture_output=True)
        if process.returncode != 0:
            raise IOError(
                f"failed to align {output_path} with {alignment_method.value} execution output is {process.stdout} "
            )
        if (
            alignment_method == AlignmentMethod.MAFFT
            and sequence_data_type == SequenceDataType.CODON
        ):
            Pipeline.reverse_translate(
                input_path,
                alignment_output_path,
                output_path,
            )
            os.remove(alignment_input_path)
            os.remove(alignment_output_path)
        return output_path

    @staticmethod
    def build_tree(
        input_path: str,
        output_path: str,
        tree_reconstruction_method: TreeReconstructionMethod,
        tree_reconstruction_params: t.Optional[t.Dict[str, t.Any]] = None,
    ):
        """
        :param input_path path to aligned sequence data in a fasta format
        :param output_path path in which the tree should be written in newick format
        :param tree_reconstruction_method: enum representing the tree reconstruction method
        :param tree_reconstruction_params: map of parameter names to parameter values
        :return: None
        """
        if tree_reconstruction_method in [
            TreeReconstructionMethod.UPGMA,
            TreeReconstructionMethod.NJ,
        ]:
            alignment = list(AlignIO.parse(input_path, "fasta"))[0]
            calculator = DistanceCalculator("identity")
            dm = calculator.get_distance(alignment)
            constructor = DistanceTreeConstructor()
            if tree_reconstruction_method == TreeReconstructionMethod.UPGMA:
                tree = constructor.upgma(dm)
            else:
                tree = constructor.nj(dm)
            with open(output_path, "w") as output_handle:
                NewickIO.write([tree], output_handle)
            tree = Tree(output_path, format=1)
            tree.write(outfile=output_path, format=5)

        elif tree_reconstruction_method == TreeReconstructionMethod.ML:
            output_dir, output_filename = os.path.split(output_path)
            aux_dir = f"{output_dir}/raxml_aux/"
            os.mkdir(aux_dir)
            os.chdir(aux_dir)
            model = (
                tree_reconstruction_params["-m"]
                if tree_reconstruction_params and "-m" in tree_reconstruction_params
                else "GTRGAMMA"
            )
            num_of_categories = (
                tree_reconstruction_params["-n"]
                if tree_reconstruction_params and "-n" in tree_reconstruction_params
                else 4
            )
            cmd = f"raxmlHPC -s {input_path} -n out -m {model} -c {num_of_categories}"
            process = subprocess.run(cmd, shell=True, capture_output=True)
            if process.returncode != 0:
                raise IOError(
                    f"failed to reconstruct ML tree with raxml. Execution output is {process.stdout}"
                )
            os.rename(f"{aux_dir}RAxML_bestTree.out", output_path)
            os.remove(aux_dir)

    def generate_samples(self, pipeline_input: PipelineInput):
        """
        :param pipeline_input instance holding parameters for pipeline execution
        :return: None. generated samples and saves them in respective directories under self.pipeline_dir
        """
        samples_dir = f"{self.pipeline_dir}/samples/"
        os.makedirs(samples_dir, exist_ok=True)
        # set sampling environment
        fraction_to_samples_dir = {}
        for fraction in pipeline_input.sampling_fractions:
            fraction_to_samples_dir[fraction] = f"{samples_dir}fraction_{fraction}/"
            os.makedirs(fraction_to_samples_dir[fraction], exist_ok=True)

        # generate the samples
        tree = Tree(self.tree_path)
        for method in pipeline_input.sampling_methods:
            sampler_instance = method_to_callable[method.value](
                sequences_path=self.unaligned_sequence_data_path,
                aux_dir=f"{samples_dir}method_{method.value}_aux/",
                tree=tree,
            )
            for fraction in pipeline_input.sampling_fractions:
                self.samples_info[fraction][method.value][
                    "unaligned_sequence_data_path"
                ] = f"{fraction_to_samples_dir[fraction]}unaligned_method_{method.value}.fasta"
                sample_size = int(fraction * len(self.unaligned_sequence_data))
                logger.info(f"Sampling data of size {sample_size} using {method.value}")
                try:
                    if method.value == "pda" and pipeline_input.weight_pda:
                        sampler_instance.compute_taxon_weights(
                            self.aligned_sequence_data_path
                        )

                    sampler_instance.write_sample(
                        sample_size,
                        output_path=self.samples_info[fraction][method.value][
                            "unaligned_sequence_data_path"
                        ],
                        is_weighted=pipeline_input.weight_pda,
                    )
                except Exception as e:
                    logger.error(
                        f"Failed to sample {sample_size} sequences with {method.value} due to error {e}"
                    )
                    raise IOError(
                        f"Failed to sample {sample_size} sequences with {method.value} due to error {e}"
                    )

                # align the sample
                self.samples_info[fraction][method.value][
                    "aligned_sequence_data_path"
                ] = f"{fraction_to_samples_dir[fraction]}aligned_method_{method.value}.fasta"
                self.samples_info[fraction][method.value][
                    "aligned_sequence_data_path"
                ] = Pipeline.align(
                    self.samples_info[fraction][method.value][
                        "unaligned_sequence_data_path"
                    ],
                    self.samples_info[fraction][method.value][
                        "aligned_sequence_data_path"
                    ],
                    pipeline_input.sequence_data_type,
                    pipeline_input.samples_alignment_method,
                    pipeline_input.samples_alignment_params,
                )

        logger.info("Completed samples generation")

    def execute_programs(self, pipeline_input: PipelineInput):
        """
        :param pipeline_input: instance holding arguments for pipeline execution
        :return: None. The input, output and result of each program will be be
        saved to self.samples_info
        """
        programs_dir = f"{self.pipeline_dir}/programs_results/"
        program_name_to_instance = dict()
        completion_validators = (
            []
        )  # this list will hold paths to touch files that should be created upon completion of the jobs holding the
        # programs executions
        for program_name in pipeline_input.programs:
            program_dir = f"{programs_dir}{program_name.value}/"
            subprocess.run(f"mkdir -p {program_dir}", shell=True, capture_output=True)
            program_to_exec = program_to_callable[program_name.value]()
            program_name_to_instance[program_name] = program_to_exec
            for fraction in self.samples_info:
                for method_name in self.samples_info[fraction]:

                    # set input and output paths
                    program_exec_info = self.samples_info[fraction][method_name][
                        "programs_performance"
                    ][program_name.value]
                    program_exec_info["input_path"] = self.samples_info[fraction][
                        method_name
                    ]["aligned_sequence_data_path"]
                    program_exec_info[
                        "output_path"
                    ] = f"{program_dir}fraction_{fraction}_sampling_method_{method_name}_output.txt"

                    # execute the program - its the same issue as with the SamplingMethod enum. Here the enum is
                    # called ProgramName
                    program_params = None
                    if (
                        pipeline_input.programs_params
                        and program_name.value in pipeline_input.programs_params
                    ):
                        program_params = pipeline_input.programs_params[
                            program_name.value
                        ]
                    if pipeline_input.parallelize:
                        completion_validator_path = program_to_exec.exec(
                            program_exec_info["input_path"],
                            program_exec_info["output_path"],
                            additional_params=program_params,
                            parallelize=pipeline_input.parallelize,
                            aux_dir=program_dir,
                            priority=pipeline_input.priority,
                            queue=pipeline_input.queue,
                            wait_until_complete=False,
                            get_completion_validator=True,
                        )
                        completion_validators.append(completion_validator_path)
                    else:
                        program_to_exec.exec(
                            program_exec_info["input_path"],
                            program_exec_info["output_path"],
                            additional_params=program_params,
                        )

        # wait for programs to finish
        if pipeline_input.parallelize:
            for validator in completion_validators:
                if not os.path.exists(validator):
                    time.sleep(5)

        # parse programs outputs
        for program_name in pipeline_input.programs:
            program_instance = program_name_to_instance[program_name]
            for fraction in self.samples_info:
                for method_name in self.samples_info[fraction]:
                    program_output_path = self.samples_info[fraction][method_name][
                        "programs_performance"
                    ][program_name.value]["output_path"]
                    self.samples_info[fraction][method_name]["programs_performance"][
                        program_name.value
                    ]["result"] = program_instance.parse_output(program_output_path)
