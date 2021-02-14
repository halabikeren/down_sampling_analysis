import json
import pickle
import shutil
import typing as t
from .pipeline_input import PipelineInput
from programs import *
import re
from dataclasses import dataclass
import time
from ete3 import Tree
from Bio import SeqIO
import os
from samplers import *
from utils import BaseTools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())

import logging

logger = logging.getLogger(__name__)


@dataclass
class Pipeline:
    pipeline_dir: str
    unaligned_sequence_data_path: str
    unaligned_sequence_data: t.List[SeqIO.SeqRecord]
    aligned_sequence_data_path: str
    tree_path: str
    new_to_orig_names_map: t.Dict[str, str]
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
        os.makedirs(processed_data_dir, exist_ok=True)
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

        # simplify input sequences names
        self.new_to_orig_names_map = BaseTools.simplify_names(input_path=self.unaligned_sequence_data_path, output_path=self.unaligned_sequence_data_path)
        new_to_orig_names_map_path = f"{pipeline_input.pipeline_dir}/new_to_orig_names_map.pickle"
        with open(new_to_orig_names_map_path, "wb") as outfile:
            pickle.dump(self.new_to_orig_names_map, outfile)

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
        logger.info(
            f"Aligning sequence data. Output will be saved to {self.aligned_sequence_data_path}"
        )
        BaseTools.align(
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
        logger.info(
            f"Reconstructing phylogenetic tree. Output will be saved to {self.tree_path}"
        )
        BaseTools.build_tree(
            self.aligned_sequence_data_path,
            self.tree_path,
            self.sequence_data_type,
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
                    "aux_dir": None,
                    "programs_performance": dict(),
                }
                for program_name in pipeline_input.programs:
                    self.samples_info[fraction][method.value]["programs_performance"][
                        program_name.value
                    ] = {
                        "input_path": None,
                        "output_path": None,
                        "aux_dir": None,
                        "result": dict(),
                        "full_data_result": None,
                    }

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
                tree=tree,
            )
            for fraction in pipeline_input.sampling_fractions:
                self.samples_info[fraction][method.value][
                    "unaligned_sequence_data_path"
                ] = f"{fraction_to_samples_dir[fraction]}unaligned_method_{method.value}.fasta"
                self.samples_info[fraction][method.value][
                    "aux_dir"
                ] = f"{fraction_to_samples_dir[fraction]}method_{method.value}_aux/"
                sample_size = int(fraction * len(self.unaligned_sequence_data))
                logger.info(
                    f"Sampling data of fraction {fraction} and size {sample_size} using {method.value}"
                )

                if (
                        os.path.exists(
                            self.samples_info[fraction][method.value][
                                "unaligned_sequence_data_path"
                            ]
                        )
                        and os.path.getsize(
                    self.samples_info[fraction][method.value][
                        "unaligned_sequence_data_path"
                    ]
                )
                        > 0
                ):
                    logger.info(
                        f"Sample of fraction {fraction} and size {sample_size} using method {method.value} already "
                        f"exists. "
                    )
                else:
                    if method.value == "pda" and pipeline_input.weight_pda:
                        sampler_instance.compute_taxon_weights(
                            self.aligned_sequence_data_path
                        )
                    sampler_instance.write_sample(
                        sample_size,
                        output_path=self.samples_info[fraction][method.value][
                            "unaligned_sequence_data_path"
                        ],
                        aux_dir=self.samples_info[fraction][method.value]["aux_dir"],
                        is_weighted=pipeline_input.weight_pda,
                        use_external=pipeline_input.use_external_pda,
                    )

                # align the sample
                self.samples_info[fraction][method.value][
                    "aligned_sequence_data_path"
                ] = f"{fraction_to_samples_dir[fraction]}aligned_method_{method.value}.fasta"
                if os.path.exists(
                        self.samples_info[fraction][method.value][
                            "aligned_sequence_data_path"
                        ]
                ) and os.path.getsize(
                    self.samples_info[fraction][method.value][
                        "aligned_sequence_data_path"
                    ]
                ) >= os.path.getsize(
                    self.samples_info[fraction][method.value][
                        "unaligned_sequence_data_path"
                    ]
                ):
                    logger.info(
                        f"Alignment of sample of fraction {fraction} and size {sample_size} using method {method.value} already "
                        f"exists. "
                    )
                else:
                    logger.info(
                        f"Aligning the sampled data to {self.samples_info[fraction][method.value]['aligned_sequence_data_path']}"
                    )
                    BaseTools.align(
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
        completion_validators = []  # this list will hold paths to touch files that should be created upon completion of the jobs holding the

        # programs executions
        program_to_full_data_output = dict()
        full_data_duration = np.nan
        for program_name in pipeline_input.programs:
            program_dir = f"{programs_dir}{program_name.value}/"
            os.makedirs(program_dir, exist_ok=True)
            program_aux_dir = f"{program_dir}/aux/"
            os.makedirs(program_aux_dir, exist_ok=True)
            program_to_exec = program_to_callable[program_name.value]()
            program_name_to_instance[program_name] = program_to_exec
            program_params = None
            if (
                    pipeline_input.programs_params
                    and program_name.value in pipeline_input.programs_params
            ):
                program_params = pipeline_input.programs_params[program_name.value]

            # execute on the full dataset
            if pipeline_input.exec_on_full_data:
                input_path = self.aligned_sequence_data_path
                output_path = f"{program_dir}/full_data_{program_name.value}.out"
                program_to_full_data_output[program_name.value] = output_path
                full_data_program_aux_dir = f"{os.path.dirname(self.aligned_sequence_data_path)}/{program_name.value}_aux/"

                if pipeline_input.parallelize:
                    completion_validator_path = program_to_exec.exec(
                        input_path=input_path,
                        output_path=output_path,
                        aux_dir=full_data_program_aux_dir,
                        additional_params=program_params,
                        parallelize=pipeline_input.parallelize,
                        cluster_data_dir=pipeline_input.cluster_data_dir,
                        priority=pipeline_input.priority,
                        queue=pipeline_input.queue,
                        wait_until_complete=False,
                        get_completion_validator=True,
                        control_file_path=f"{full_data_program_aux_dir}/input.ctl",
                        input_tree_path=f"{full_data_program_aux_dir}/tree.nwk"
                    )
                    completion_validators.append(completion_validator_path)
                else:
                    full_data_duration = program_to_exec.exec(
                        input_path=input_path,
                        output_path=output_path,
                        aux_dir=full_data_program_aux_dir,
                        additional_params=program_params,
                        control_file_path=f"{full_data_program_aux_dir}/input.ctl",
                        input_tree_path=f"{full_data_program_aux_dir}/tree.nwk"
                    )

            # execute program on each sample
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
                    program_exec_info[
                        "aux_dir"
                    ] = f"{program_dir}fraction_{fraction}_sampling_method_{method_name}_aux/"

                    if pipeline_input.parallelize:
                        completion_validator_path = program_to_exec.exec(
                            program_exec_info["input_path"],
                            program_exec_info["output_path"],
                            program_exec_info["aux_dir"],
                            additional_params=program_params,
                            parallelize=pipeline_input.parallelize,
                            cluster_data_dir=pipeline_input.cluster_data_dir,
                            priority=pipeline_input.priority,
                            queue=pipeline_input.queue,
                            wait_until_complete=False,
                            get_completion_validator=True,
                        )
                        completion_validators.append(completion_validator_path)
                    else:
                        duration = program_to_exec.exec(
                            program_exec_info["input_path"],
                            program_exec_info["output_path"],
                            program_exec_info["aux_dir"],
                            additional_params=program_params,
                        )
                        self.samples_info[fraction][method_name]["programs_performance"][program_name.value][
                            "result"].update({"duration(minutes)": duration})

        # wait for programs to finish
        if pipeline_input.parallelize:
            for validator in completion_validators:
                if validator:
                    while not os.path.exists(validator):
                        time.sleep(5)

        # parse programs outputs
        for program_name in pipeline_input.programs:
            program_instance = program_name_to_instance[program_name]

            # parse output of the full program
            if pipeline_input.exec_on_full_data:
                full_program_output = program_to_full_data_output[program_name.value]
                full_data_result = program_instance.parse_output(output_path=full_program_output, job_output_dir=(
                    full_data_program_aux_dir if pipeline_input.parallelize else None))
                if "duration(minutes)" not in full_data_result:
                    full_data_result["duration(minutes)"] = full_data_duration

            for fraction in self.samples_info:
                for method_name in self.samples_info[fraction]:
                    program_output_path = self.samples_info[fraction][method_name][
                        "programs_performance"
                    ][program_name.value]["output_path"]
                    job_output_dir = \
                        self.samples_info[fraction][method_name]["programs_performance"][program_name.value]["aux_dir"]
                    self.samples_info[fraction][method_name]["programs_performance"][
                        program_name.value
                    ]["result"].update(program_instance.parse_output(output_path=program_output_path,
                                                                     job_output_dir=job_output_dir if pipeline_input.parallelize else None))
                    if pipeline_input.exec_on_full_data:
                        self.samples_info[fraction][method_name]["programs_performance"][
                            program_name.value
                        ]["full_data_result"] = full_data_result

    def write_results(self, output_path: str):
        if os.path.exists(output_path):
            logger.info(f"written output already exists at {output_path}")
            return
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, "w") as outfile:
            json.dump(self.samples_info, outfile)

    def analyze_results(self, pipeline_input: PipelineInput):
        """
        :param pipeline_input: pipeline input instance
        :return: nothing. analyses programs outputs and writes figures with the result to output dir
        """
        output_dir = f"{pipeline_input.pipeline_dir}/figures"
        os.makedirs(output_dir, exist_ok=True)
        programs = pipeline_input.programs
        sampling_methods = pipeline_input.sampling_methods
        sampling_fractions = pipeline_input.sampling_fractions
        for program_name in programs:
            figure_path = f"{output_dir}/{program_name.value}.svg"
            plt.grid(False)
            fig = plt.figure(figsize=[1 * 8.5 + 2, 1 * 7.58 + 2], frameon=True)
            accuracy_dfs = []
            for fraction in sampling_fractions:
                for method in sampling_methods:
                    result_data = self.samples_info[fraction][method.value]["programs_performance"][program_name.value]
                    full_data_result = result_data["full_data_result"]
                    sampled_data_result = result_data["result"]
                    comparison_df = pd.DataFrame(columns=["sampling_fraction", "sampling_method", "accuracy"])
                    comparison_df["accuracy"] = program_to_callable[program_name.value].get_accuracy(reference_data=full_data_result, test_data=sampled_data_result)
                    comparison_df["sampling_fraction"] = fraction
                    comparison_df["sampling_method"] = method.value
                    accuracy_dfs.append(comparison_df)
            accuracy_df = pd.concat(accuracy_dfs)
            sns.boxplot(y="accuracy", x="sampling_fraction", data=accuracy_df, palette="colorblind", hue="sampling_method")
            fig.subplots_adjust()
            fig.tight_layout()
            plt.savefig(figure_path, bbox_inches="tight", transparent=True)
            plt.clf()




