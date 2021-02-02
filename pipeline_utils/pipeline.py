import typing as t
from pipeline_utils import (
    SequenceDataType,
    AlignmentMethod,
    TreeReconstructionMethod,
    PipelineInput,
)
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
from pydantic import FilePath
from samplers import SamplingMethod

import logging

logger = logging.getLogger(__name__)


@dataclass
class Pipeline:
    pipeline_dir: str
    unaligned_sequence_data_path: str
    unaligned_sequence_data: t.List[SeqIO.SeqRecord]
    aligned_sequence_data_path: str
    tree_path: str
    sequence_data_type: SequenceDataType
    samples_info: t.Dict[
        float, t.Dict[SamplingMethod, t.Dict[str, t.Any]]
    ]  # will map sample fractions to maps of sampling methods to their info (input path, output paths,
    # analysis results ect.)
    weight_pda: bool = False
    parallelize: bool = True
    priority: int = 0
    queue: str = "itaym"

    def __init__(self, pipeline_input: PipelineInput):

        # set initial parameters
        dataset_name_regex = re.compile("(.*)\.")
        dataset_name = dataset_name_regex.search(
            os.path.basename(pipeline_input.sequence_data_path)
        ).group(1)
        self.pipeline_dir = pipeline_input.pipeline_dir
        processed_data_dir = f"{self.pipeline_dir}/input_data/"
        subprocess.run(
            f"mkdir -p {processed_data_dir}", shell=True, capture_output=True
        )

        logger.info(f"Setting input for pipeline at {processed_data_dir}")

        self.unaligned_sequence_data_path = (
            f"{processed_data_dir}{dataset_name}_unaligned.fasta"
        )
        logger.info(
            f"Unaligned sequence data saved at {self.unaligned_sequence_data_path}"
        )
        self.aligned_sequence_data_path = (
            f"{processed_data_dir}{dataset_name}_aligned.fasta"
        )
        logger.info(f"Aligned sequence data saved at {self.aligned_sequence_data_path}")

        self.tree_path = f"{processed_data_dir}{dataset_name}_tree.nwk"
        if pipeline_input.tree_path:
            subprocess.run(f"cp -r {pipeline_input.tree_path} {self.tree_path}")
        logger.info(f"Tree saved at {self.tree_path}")

        # fill in available parameters
        self.sequence_data_type = pipeline_input.sequence_data_type
        if Pipeline.is_aligned(pipeline_input.sequence_data_path):
            subprocess.run(
                f"cp -r {pipeline_input.sequence_data_path} {self.aligned_sequence_data_path}",
                shell=True,
                capture_output=True,
            )
            self.un_align(str(pipeline_input.sequence_data_path))
        else:
            subprocess.run(
                f"cp -r {pipeline_input.sequence_data_path} {self.unaligned_sequence_data_path}",
                shell=True,
                capture_output=True,
            )
            self.align(pipeline_input.alignment_method, pipeline_input.alignment_params)

        logger.info(
            f"Alignment generated successfully using {pipeline_input.alignment_method.value}"
        )

        self.unaligned_sequence_data = list(
            SeqIO.parse(self.unaligned_sequence_data_path, "fasta")
        )
        logger.info(
            f"Processed sequence data of size {len(self.unaligned_sequence_data)}"
        )

        if not os.path.exists(self.tree_path):
            self.build_tree(
                pipeline_input.tree_reconstruction_method,
                pipeline_input.tree_reconstruction_params,
            )
        logger.info(
            f"Tree reconstructed successfully using {pipeline_input.tree_reconstruction_method.value}"
        )

        # set sampling info structure
        self.samples_info = dict()
        for fraction in pipeline_input.sampling_fractions:
            self.samples_info[fraction] = dict()
            for method in pipeline_input.sampling_methods:
                self.samples_info[fraction][method.value[0]] = {
                    "path": None,
                    "program_performance": dict(),
                }
                for program in pipeline_input.programs:
                    self.samples_info[fraction][method.value[0]]["program_performance"][
                        program
                    ] = {"input_dir": None, "output_dir": None, "result": None}

    @staticmethod
    def is_aligned(sequence_data: t.Union[FilePath, t.List[SeqIO.SeqRecord]]) -> bool:
        """
        :param sequence_data sequence data in the form of a filepath or sequence records
        :return: a boolean indicating weather the records are aligned or not
        """
        if type(sequence_data) is not list:
            sequence_data = list(SeqIO.parse(sequence_data, "fasta"))
        for record in sequence_data:
            if "-" in record.seq:
                return True
        return False

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

    def un_align(self, input_path: str):
        """
        :param input_path: sequence data in the form of a filepath
        :return: None. writes unaligned data to the respective file
        """
        if not self.unaligned_sequence_data:
            self.unaligned_sequence_data = list(SeqIO.parse(input_path, "fasta"))
        for record in self.unaligned_sequence_data:
            record.seq = str(record.seq).replace("-", "")
        SeqIO.write(
            self.unaligned_sequence_data, self.unaligned_sequence_data_path, "fasta"
        )

    def align(
        self,
        alignment_method: AlignmentMethod,
        alignment_params: t.Dict[str, str] = None,
    ):
        """
        Aligns sequence data according to given method
        :return:
        """
        if os.path.exists(self.aligned_sequence_data_path):
            alternative_alignment_path = self.aligned_sequence_data_path.replace(
                ".fasta", f"_{time.time()}.fasta"
            )
            logger.warning(
                f"{self.aligned_sequence_data_path} already exists. The program will create the alignment in the "
                f"alternative path {alternative_alignment_path} "
            )
            self.aligned_sequence_data_path = alternative_alignment_path

        sequence_records = list(SeqIO.parse(self.unaligned_sequence_data_path, "fasta"))
        alignment_input_path = self.unaligned_sequence_data_path
        alignment_output_path = self.aligned_sequence_data_path
        if (
            self.sequence_data_type == SequenceDataType.CODON
            and alignment_method == AlignmentMethod.MAFFT
        ):
            alignment_input_path = self.unaligned_sequence_data_path.replace(
                ".fasta", "_translated.fasta"
            )
            self.translate(
                sequence_records,
                self.unaligned_sequence_data_path.replace(
                    ".fasta", "_translated.fasta"
                ),
            )
            alignment_output_path = self.aligned_sequence_data_path.replace(
                ".fasta", "_translated.fasta"
            )

        cmd = ""
        if alignment_method == AlignmentMethod.MAFFT:
            cmd = f"mafft --localpair --maxiterate 1000 {alignment_input_path} > {alignment_output_path}"
        elif alignment_method == AlignmentMethod.PRANK:
            cmd = f"prank -d={alignment_input_path} -o={alignment_output_path} -f=fasta -support {'-codon' if self.sequence_data_type == SequenceDataType.CODON else ''} -iterate=100 -showtree "
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
                f"failed to align {self.unaligned_sequence_data_path} with {alignment_method.value} execution output is {process.stdout} "
            )
        if (
            alignment_method == AlignmentMethod.MAFFT
            and self.sequence_data_type == SequenceDataType.CODON
        ):
            self.reverse_translate(
                self.unaligned_sequence_data_path,
                alignment_output_path,
                self.aligned_sequence_data_path,
            )
            os.remove(alignment_input_path)
            os.remove(alignment_output_path)

    def build_tree(
        self,
        tree_reconstruction_method: TreeReconstructionMethod,
        tree_reconstruction_params: t.Dict[str, str] = None,
    ):
        """
        :param tree_reconstruction_method: enum representing the tree reconstruction method
        :param tree_reconstruction_params: map of parameter names to parameter values
        :return: None
        """
        if tree_reconstruction_method in [
            TreeReconstructionMethod.UPGMA,
            TreeReconstructionMethod.NJ,
        ]:
            alignment = list(AlignIO.parse(self.aligned_sequence_data_path, "fasta"))[0]
            calculator = DistanceCalculator("identity")
            dm = calculator.get_distance(alignment)
            constructor = DistanceTreeConstructor()
            if tree_reconstruction_method == TreeReconstructionMethod.UPGMA:
                tree = constructor.upgma(dm)
            else:
                tree = constructor.nj(dm)
            with open(self.tree_path, "w") as output_handle:
                NewickIO.write([tree], output_handle)
            tree = Tree(self.tree_path, format=1)
            tree.write(outfile=self.tree_path, format=5)

        elif tree_reconstruction_method == TreeReconstructionMethod.ML:
            output_dir, output_filename = os.path.split(self.tree_path)
            aux_dir = f"{output_dir}/raxml_aux/"
            os.mkdir(aux_dir)
            os.chdir(aux_dir)
            model = (
                tree_reconstruction_params["-m"]
                if "-m" in tree_reconstruction_params
                else "GTRGAMMA"
            )
            num_of_categories = (
                tree_reconstruction_params["-n"]
                if "-n" in tree_reconstruction_params
                else 4
            )
            cmd = f"raxmlHPC -s {self.aligned_sequence_data_path} -n out -m {model} -c {num_of_categories}"
            process = subprocess.run(cmd, shell=True, capture_output=True)
            if process.returncode != 0:
                raise IOError(
                    f"failed to reconstruct ML tree with raxml. Execution output is {process.stdout}"
                )
            os.rename(f"{aux_dir}RAxML_bestTree.out", self.tree_path)
            os.remove(aux_dir)

    def generate_samples(
        self,
        sampling_fractions: t.List[float],
        sampling_methods: t.List[SamplingMethod],
    ):
        """
        :param sampling_fractions: the fractions of data that should be sampled
        :param sampling_methods: the methods that the data should be sampled with
        :return: None. generated samples and saves them in respective directories under self.pipeline_dir
        """
        samples_dir = f"{self.pipeline_dir}/samples/"
        subprocess.run(f"mkdir -p {samples_dir}", shell=True, capture_output=True)
        # set sampling environment
        fraction_to_samples_dir = {}
        for fraction in sampling_fractions:
            fraction_to_samples_dir[fraction] = f"{samples_dir}fraction_{fraction}/"
            subprocess.run(
                f"mkdir -p {fraction_to_samples_dir[fraction]}",
                shell=True,
                capture_output=True,
            )

        # generate the samples
        tree = Tree(self.tree_path)
        for method in sampling_methods:
            sampler = method.value[1](
                sequences_path=self.unaligned_sequence_data_path,
                aux_dir=f"{samples_dir}method_{method.value[0]}_aux/",
                tree=tree,
            )
            for fraction in sampling_fractions:
                self.samples_info[fraction][method.value[0]][
                    "path"
                ] = f"{fraction_to_samples_dir[fraction]}method_{method.value[0]}.fasta"
                sample_size = int(fraction * len(self.unaligned_sequence_data))
                logger.info(
                    f"Sampling data of size {sample_size} using {method.value[0]}"
                )
                try:
                    sampler.write_sample(
                        sample_size,
                        output_path=self.samples_info[fraction][method.value[0]][
                            "path"
                        ],
                    )
                except Exception as e:
                    logger.error(
                        f"Failed to sample {sample_size} sequences with {method.value[0]} due to error {e}"
                    )
                    raise IOError(
                        f"Failed to sample {sample_size} sequences with {method.value[0]} due to error {e}"
                    )

        logger.info("Completed samples generation")

    # def execute_programs(
    #         self,
    #         programs: t.List[str],
    #         programs_params: t.Optional[t.Dict[str, t.Dict[str, t.Any]]] = None,
    #         parallelize: bool = True,
    #         priority: int = 0,
    #         queue: str = "itaym",
    # ):
    #     """
    #     :param programs: list of program executables or aliases of them, that should be executed on the samples in
    #     self.samples_info
    #     :param programs_params:
    #     :param parallelize: a boolean indicating weather execution of the programs on the samples
    #     should be parallelized (i.e., for each combo of program, fraction and method of sampling, a job will be
    #     generated) or should be run in turns
    #     :param priority:
    #     :param queue:
    #     :return: None. The input, output and result of each program will be be
    #     saved to self.samples_info
    #     """
    #
    #     for fraction in self.samples_info:
    #         for method in self.samples_info[fraction]:
    #             sample_path = self.samples_info[fraction][method.value[0]]["path"]
    #             for program in programs:
    #                 params = programs_params[program]
    #                 # set input from the program
    #
    #                 # construct command line
    #
    #                 # execute the program

    # def run_program( self, program: str, input_parameters: t.Dict[str, str], working_dir: str = None, ): """ :param
    # program: path to te executable of a program :param input_parameters: dictionary that maps input argument name
    # to input argument value :param input_parameters dictionary that maps input parameter names to values for the
    # program :param working_dir: a directory to create and run the jobs in, if needed :return: the output of the
    # program in a string format, if exists, or none (in case the output is written to a file, for example """
    # command = f"{program}" + " ".join( [ f"{param_name} {input_parameters[param_name]}" for param_name in
    # input_parameters ] ) if not self.parallelize: process = subprocess.run(command, shell=True,
    # capture_output=True) if process.returncode != 0: raise RuntimeError( f"Program {program} failed to execute with
    # error {process.stderr} and output {process.stdout}" ) else: subprocess.run(f"mkdir -p {working_dir}",
    # shell=True, capture_output=True) self.submit_as_job(command, working_dir, self.priority, self.queue)

    def submit_as_job(
        self, command: str, working_dir: str, priority: int = 0, queue="itaym"
    ):
        """
        executes a command via a job
        :param command: string of a command to run as shell
        :param working_dir: directory to hold the job and job output files
        :param priority:
        :param queue:
        :return: None
        """
        error_filepath = f"{working_dir}job_output"
        job_filepath = f"{working_dir}job.sh"
        job_name = command.split(" ")[0]
        self.create_job_file(
            job_name, [command], error_filepath, job_filepath, queue=queue
        )
        subprocess.run(f"qsub -p {priority} {job_filepath}")

    @staticmethod
    def create_job_file(
        job_name,
        commands,
        error_filepath,
        job_filepath,
        cpus_number=1,
        queue="itaym",
        mem_alloc="2gb",
    ):
        with open(job_filepath, "w") as handle:
            handle.write(
                f"# !/bin/bash\n\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -r y\n#PBS -q {queue}\n"
            )
            handle.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n")
            handle.write(
                f"#PBS -N {job_name}\n#PBS -e {error_filepath}\n#PBS -o {error_filepath}\n"
            )
            handle.write(f"#PBS -l select=ncpus={cpus_number}:mem={mem_alloc}\n")
            handle.write("\n".join(commands))
