import itertools
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
from copy import deepcopy
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
        self.sequence_data_type = pipeline_input.sequence_data_type

        # prepare input for pipeline
        processed_data_dir = f"{self.pipeline_dir}/input_data/"
        os.makedirs(processed_data_dir, exist_ok=True)
        logger.info(f"Setting input for pipeline at {processed_data_dir}")

        # save unaligned sequence data
        self.unaligned_sequence_data_path = (
            f"{processed_data_dir}{dataset_name}_unaligned.fasta"
        )
        # simplify input sequences names
        self.new_to_orig_names_map = BaseTools.simplify_names(input_path=pipeline_input.unaligned_sequence_data_path,
                                                              output_path=self.unaligned_sequence_data_path)
        new_to_orig_names_map_path = f"{pipeline_input.pipeline_dir}/new_to_orig_names_map.pickle"
        with open(new_to_orig_names_map_path, "wb") as outfile:
            pickle.dump(self.new_to_orig_names_map, outfile)
        logger.info(
            f"Unaligned sequence data saved at {self.unaligned_sequence_data_path}"
        )
        self.unaligned_sequence_data = list(
            SeqIO.parse(self.unaligned_sequence_data_path, "fasta")
        )
        logger.info(
            f"Processed sequence data of size {len(self.unaligned_sequence_data)}"
        )

        # save aligned sequence data
        self.aligned_sequence_data_path = (
            f"{processed_data_dir}{dataset_name}_aligned.fasta"
        )
        if pipeline_input.aligned_sequence_data_path:
            BaseTools.simplify_names(pipeline_input.aligned_sequence_data_path, self.aligned_sequence_data_path,
                                     self.new_to_orig_names_map)
            logger.info(f"Aligned data saved at {self.aligned_sequence_data_path}")
        else:
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

        # save tree
        self.tree_path = f"{processed_data_dir}{dataset_name}_tree.nwk"
        if pipeline_input.tree_path:
            BaseTools.simplify_names(pipeline_input.tree_path, self.tree_path, self.new_to_orig_names_map)
            logger.info(f"Tree saved at {self.tree_path}")
        else:
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
                    "tree_path": None,
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
                        "reference_data": None
                    }

    def align_sampled_data(self, pipeline_input: PipelineInput, fraction: float, method: SamplingMethod,
                           fraction_to_samples_dir: t.Dict[t.Any, t.Any], sample_member_names: t.List[str]):
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
                f"Alignment of sample of fraction {fraction} using method {method.value} already exists."
            )
        else:
            if pipeline_input.use_full_alignment_in_sample:
                logger.info("Trimming full alignment to include only the sampled sequences")

                all_aligned_sequence_records = list(SeqIO.parse(self.aligned_sequence_data_path, "fasta"))
                sampled_aligned_sequence_records = [record for record in all_aligned_sequence_records if
                                                    record.name in sample_member_names]
                SeqIO.write(sampled_aligned_sequence_records, self.samples_info[fraction][method.value][
                    "aligned_sequence_data_path"
                ], "fasta")
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
                logger.info(
                    f"Aligned sample records written to {self.samples_info[fraction][method.value]['aligned_sequence_data_path']}"
                )

    def build_sampled_tree(self, pipeline_input: PipelineInput, fraction: float, method: SamplingMethod,
                           fraction_to_samples_dir: t.Dict[t.Any, t.Any], sample_member_names: t.List[str]):
        self.samples_info[fraction][method.value][
            "tree_path"
        ] = f"{fraction_to_samples_dir[fraction]}method_{method.value}_tree.nwk"
        if os.path.exists(
                self.samples_info[fraction][method.value][
                    "tree_path"
                ]
        ) and os.path.getsize(
            self.samples_info[fraction][method.value][
                "tree_path"
            ]) > 0:
            logger.info(
                f"Tree of sample of fraction {fraction} using method {method.value} already exists."
            )
        else:
            if pipeline_input.use_full_tree_in_sample:
                logger.info("Pruning full tree to include only the sampled sequences")
                full_tree = Tree(self.tree_path)
                full_tree.prune(sample_member_names)
                full_tree.write(outfile=self.samples_info[fraction][method.value]["tree_path"])
            else:
                logger.info(
                    f"Building tree based on sampled data to {self.samples_info[fraction][method.value]['tree_path']}"
                )
                BaseTools.build_tree(
                    self.samples_info[fraction][method.value][
                        "aligned_sequence_data_path"
                    ],
                    self.samples_info[fraction][method.value][
                        "tree_path"
                    ],
                    pipeline_input.sequence_data_type,
                    pipeline_input.tree_reconstruction_method,
                    pipeline_input.tree_reconstruction_params,
                )
                logger.info(
                    f"Tree of sample records written to {self.samples_info[fraction][method.value]['tree_path']}"
                )

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
        for method in pipeline_input.sampling_methods:
            sampler_instance = method_to_callable[method.value](
                sequence_data_path=self.unaligned_sequence_data_path,
                tree_path=self.tree_path, exclude_a_ref_sequence=pipeline_input.exclude_ref_seq
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
                    sample_members_names = [record.name for record in
                                            list(SeqIO.parse(self.samples_info[fraction][method.value][
                                                                 "unaligned_sequence_data_path"], "fasta"))]
                else:
                    if method.value == "pda" and pipeline_input.weight_pda:
                        sampler_instance.compute_taxon_weights(
                            self.aligned_sequence_data_path
                        )
                    sample = sampler_instance.write_sample(
                        sample_size,
                        output_path=self.samples_info[fraction][method.value][
                            "unaligned_sequence_data_path"
                        ],
                        aux_dir=self.samples_info[fraction][method.value]["aux_dir"],
                        is_weighted=pipeline_input.weight_pda,
                        use_external=pipeline_input.use_external_pda,
                    )
                    if type(sample) is str:
                        sample = list(SeqIO.parse(sample, "fasta"))
                    sample_members_names = [record.name for record in sample]

                # align the sample
                self.align_sampled_data(pipeline_input=pipeline_input, fraction=fraction, method=method,
                                        fraction_to_samples_dir=fraction_to_samples_dir,
                                        sample_member_names=sample_members_names)

                # write the tree
                self.build_sampled_tree(pipeline_input=pipeline_input, fraction=fraction, method=method,
                                        fraction_to_samples_dir=fraction_to_samples_dir,
                                        sample_member_names=sample_members_names)

        logger.info("Completed samples generation")

        # compare % overlap between samples of the same size
        df = pd.DataFrame(columns=["sampling_fraction", "method_1", "method_2", "overlap_fraction"])
        for fraction in pipeline_input.sampling_fractions:
            method_to_sample = dict()
            for method in pipeline_input.sampling_methods:
                sample_path = self.samples_info[fraction][method.value]["unaligned_sequence_data_path"]
                method_to_sample[method.value] = [record.name for record in list(SeqIO.parse(sample_path, "fasta"))]
            for (method_1, method_2) in [pair for pair in itertools.combinations(list(method_to_sample.keys()), 2)]:
                overlap_fraction = len([record for record in method_to_sample[method_1] if record in method_to_sample[method_2]]) / len(method_to_sample[method_1])
                df_record = {"sampling_fraction": fraction, "method_1": method_1, "method_2": method_2, "overlap_fraction": overlap_fraction}
                df = df.append(df_record, ignore_index=True)
            df.to_csv(f"{samples_dir}/samples_overlap.csv")

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
            program_params = dict()
            if (
                    pipeline_input.programs_params
                    and program_name.value in pipeline_input.programs_params
            ):
                program_params = pipeline_input.programs_params[program_name.value]
            if program_name.value == "rate4site" and pipeline_input.exclude_ref_seq:
                program_params["-a"] = self.unaligned_sequence_data[0].name

            # execute on the full dataset
            if pipeline_input.exec_on_full_data:
                input_path = self.aligned_sequence_data_path
                output_path = f"{program_dir}/full_data_{program_name.value}.out"
                program_to_full_data_output[program_name.value] = output_path
                full_data_program_aux_dir = f"{program_dir}/full_data_{program_name.value}_aux/"
                full_data_program_params = dict()
                if program_params:
                    full_data_program_params = deepcopy(program_params)
                if program_name.value == "paml":
                    full_data_program_params["input_tree_path"] = self.tree_path

                if pipeline_input.parallelize:
                    completion_validator_path = program_to_exec.exec(
                        input_path=input_path,
                        output_path=output_path,
                        aux_dir=full_data_program_aux_dir,
                        additional_params=full_data_program_params,
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
                        additional_params=full_data_program_params,
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
                try:
                    full_program_output = program_to_full_data_output[program_name.value]
                    full_data_result = program_instance.parse_output(output_path=full_program_output, job_output_dir=
                    full_data_program_aux_dir)
                    if "duration(minutes)" not in full_data_result:
                        full_data_result["duration(minutes)"] = full_data_duration
                except Exception as e:
                    logger.error(
                        f"Failed to parse results on full data due to error {e}. check why it failed in {full_data_program_aux_dir}")
                    full_data_result = None

            for fraction in self.samples_info:
                for method_name in self.samples_info[fraction]:
                    program_output_path = self.samples_info[fraction][method_name][
                        "programs_performance"
                    ][program_name.value]["output_path"]
                    job_output_dir = \
                        self.samples_info[fraction][method_name]["programs_performance"][program_name.value]["aux_dir"]
                    if os.path.exists(job_output_dir):
                        self.samples_info[fraction][method_name]["programs_performance"][
                            program_name.value
                        ]["result"].update(program_instance.parse_output(output_path=program_output_path,
                                                                         job_output_dir=job_output_dir))
                    if pipeline_input.exec_on_full_data:
                        self.samples_info[fraction][method_name]["programs_performance"][
                            program_name.value
                        ]["full_data_result"] = full_data_result

                    if pipeline_input.reference_data_paths and program_name.value in pipeline_input.reference_data_paths:
                        self.samples_info[fraction][method_name]["programs_performance"][
                            program_name.value
                        ]["reference_data"] = program_instance.parse_reference_data(
                            pipeline_input.reference_data_paths[program_name.value])

    def write_results(self, output_path: str):
        if os.path.exists(output_path):
            logger.info(f"written output already exists at {output_path}")
            return
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, "w") as outfile:
            json.dump(self.samples_info, outfile)

    @staticmethod
    def get_result_df(program_instance: Program, reference_data: t.Dict[str, t.Any], test_data: t.Dict[str, t.Any],
                      fraction: float, method: SamplingMethod, use_normalized_rates: bool = False,
                      use_relative_error: bool = False):
        """
        :param program_instance: program instance to get its error computation
        :param reference_data: data of reference set
        :param test_data: data of test set
        :param fraction: sampling fraction
        :param method: sampling method
        :param use_normalized_rates: for rate4site. indicates weather normalized rates should be used or not
        :param use_relative_error: boolean indicating weather relative or absolute error should be used
        :return:
        """
        df = pd.DataFrame(columns=["sampling_fraction", "sampling_method", "error", "result", "reference"])
        df["error"] = program_instance.get_error(
            reference_data=reference_data, test_data=test_data, use_normalized_rates=use_normalized_rates,
            use_relative_error=use_relative_error)
        df["result"] = program_instance.get_result(data=test_data)
        df["reference"] = program_instance.get_result(data=reference_data)
        df["sampling_fraction"] = fraction
        df["sampling_method"] = method.value
        return df

    def plot_error(self, pipeline_input, program_name: ProgramName, output_path: str, use_relative_error: bool = False):
        plt.grid(False)
        ncols = 1 if not pipeline_input.reference_data_paths[program_name.value] else 2
        fig, axis = plt.subplots(
            nrows=1,
            ncols=ncols,
            sharex="none",
            sharey="none",
            figsize=[ncols * 8.5 + 2 + 2, 7.58 + 2],
            frameon=True,
        )
        full_error_dfs = []
        ref_error_dfs = []
        for fraction in pipeline_input.sampling_fractions:
            for method in pipeline_input.sampling_methods:
                result_data = self.samples_info[fraction][method.value]["programs_performance"][program_name.value]
                full_data_result = result_data["full_data_result"]
                reference_data_result = result_data["reference_data"]
                sampled_data_result = result_data["result"]
                if full_data_result:
                    full_error_dfs.append(
                        Pipeline.get_result_df(program_instance=program_to_callable[program_name.value],
                                               reference_data=full_data_result, test_data=sampled_data_result,
                                               fraction=fraction, method=method, use_relative_error=use_relative_error))
                if reference_data_result:
                    ref_error_dfs.append(
                        Pipeline.get_result_df(program_instance=program_to_callable[program_name.value],
                                               reference_data=reference_data_result,
                                               test_data=sampled_data_result,
                                               fraction=fraction, method=method, use_normalized_rates=False,
                                               use_relative_error=use_relative_error))
        if len(full_error_dfs) > 0:
            full_error_df = pd.concat(full_error_dfs)
            sns.boxplot(ax=axis[0], y="error", x="sampling_fraction", data=full_error_df,
                        palette="colorblind",
                        hue="sampling_method")
            axis[0].set_title("reference: full data")
        if len(ref_error_dfs) > 0:
            ref_error_df = pd.concat(ref_error_dfs)
            sns.boxplot(ax=axis[1], y="error", x="sampling_fraction", data=ref_error_df, palette="colorblind",
                        hue="sampling_method")
            axis[1].set_title("reference: simulated data")
        fig.subplots_adjust()
        fig.tight_layout()
        plt.savefig(output_path, bbox_inches="tight", transparent=True)
        plt.clf()

    def plot_bias(self, pipeline_input, program_name: ProgramName, output_path: str):
        plt.grid(False)
        ncols = 1 if not pipeline_input.reference_data_paths[program_name.value] else 2
        fig, axis = plt.subplots(
            nrows=1,
            ncols=ncols,
            sharex="none",
            sharey="none",
            figsize=[ncols * 8.5 + 2 + 2, 7.58 + 2],
            frameon=True,
        )
        full_error_dfs = []
        ref_error_dfs = []
        for fraction in pipeline_input.sampling_fractions:
            for method in pipeline_input.sampling_methods:
                result_data = self.samples_info[fraction][method.value]["programs_performance"][program_name.value]
                full_data_result = result_data["full_data_result"]
                reference_data_result = result_data["reference_data"]
                sampled_data_result = result_data["result"]
                if full_data_result:
                    full_error_dfs.append(
                        Pipeline.get_result_df(program_instance=program_to_callable[program_name.value],
                                               reference_data=full_data_result, test_data=sampled_data_result,
                                               fraction=fraction, method=method, use_relative_error=True))
                if reference_data_result:
                    ref_error_dfs.append(
                        Pipeline.get_result_df(program_instance=program_to_callable[program_name.value],
                                               reference_data=reference_data_result, test_data=sampled_data_result,
                                               fraction=fraction, method=method, use_relative_error=True))
        if len(full_error_dfs) > 0:
            full_error_df = pd.concat(full_error_dfs)
            full_error_df["bias"] = full_error_df["reference"] - full_error_df["result"]
            sns.boxplot(ax=axis[0], y="bias", x="sampling_fraction", data=full_error_df,
                        palette="colorblind",
                        hue="sampling_method")
            axis[0].set_ylabel("full_rate-sampled_rate")
            axis[0].set_title("reference: full")
        if len(ref_error_dfs) > 0:
            ref_error_df = pd.concat(ref_error_dfs)
            ref_error_df["bias"] = ref_error_df["reference"] - ref_error_df["result"]
            sns.boxplot(ax=axis[1], y="bias", x="sampling_fraction", data=ref_error_df,
                        palette="colorblind",
                        hue="sampling_method")
            axis[1].set_ylabel("simulated_rate-sampled_rate")
            axis[1].set_title("reference: simulated")
        fig.subplots_adjust()
        fig.tight_layout()
        plt.savefig(output_path, bbox_inches="tight", transparent=True)
        plt.clf()

    def plot_results(self, pipeline_input: PipelineInput):
        """
        :param pipeline_input: pipeline input instance
        :return: nothing. analyses programs outputs and writes figures with the result to output dir
        """
        output_dir = f"{pipeline_input.pipeline_dir}/figures"
        os.makedirs(output_dir, exist_ok=True)
        for program_name in pipeline_input.programs:
            self.plot_error(pipeline_input=pipeline_input, program_name=program_name,
                            output_path=f"{output_dir}/{program_name.value}_absolute_error.svg",
                            use_relative_error=False)
            self.plot_error(pipeline_input=pipeline_input, program_name=program_name,
                            output_path=f"{output_dir}/{program_name.value}_relative_error.svg",
                            use_relative_error=True)
            self.plot_bias(pipeline_input=pipeline_input, program_name=program_name,
                            output_path=f"{output_dir}/{program_name.value}_bias.svg",)

    def analyse_results(self, pipeline_input: PipelineInput):
        """
        :param pipeline_input: pipeline input instance
        :return: none. writes error data of program to output csv
        """
        output_dir = f"{pipeline_input.pipeline_dir}/tables"
        os.makedirs(output_dir, exist_ok=True)
        for program_name in pipeline_input.programs:
            program_class = program_to_callable[program_name.value]
            output_path = f"{output_dir}/{program_name.value}_summary.csv"
            output_dfs = []
            for fraction in pipeline_input.sampling_fractions:
                for method in pipeline_input.sampling_methods:
                    df = pd.DataFrame(columns=["sampling_fraction", "sampling_method", "relative_error_to_ref",
                                               "relative_error_to_full", "absolute_error_to_ref",
                                               "absolute_error_to_full"])
                    sample_result = \
                        self.samples_info[fraction][method.value]["programs_performance"][program_name.value]["result"]
                    full_result = self.samples_info[fraction][method.value]["programs_performance"][program_name.value][
                        "full_data_result"]
                    reference_result = \
                        self.samples_info[fraction][method.value]["programs_performance"][program_name.value][
                            "reference_data"]
                    df["result"] = program_class.get_result(sample_result)
                    df["simulated"] = program_class.get_result(reference_result)
                    df["full_result"] = program_class.get_result(full_result)
                    try:
                        df["relative_error_to_ref"] = program_class.get_error(reference_data=reference_result,
                                                                              test_data=sample_result,
                                                                              use_relative_error=True)
                        df["absolute_error_to_ref"] = program_class.get_error(reference_data=reference_result,
                                                                              test_data=sample_result,
                                                                              use_relative_error=False)
                    except Exception as e:
                        logger.error(
                            f"Could not compute error of sample {fraction}_{method.value} relative to reference due to error {e}")
                    try:
                        df["relative_error_to_full"] = program_class.get_error(reference_data=full_result,
                                                                               test_data=sample_result,
                                                                               use_relative_error=True)
                        df["absolute_error_to_full"] = program_class.get_error(reference_data=full_result,
                                                                               test_data=sample_result,
                                                                               use_relative_error=False)
                    except Exception as e:
                        logger.error(
                            f"Could not compute error of sample {fraction}_{method.value} relative to full due to error {e}")
                    df["sampling_fraction"] = fraction
                    df["sampling_method"] = method.value
                    output_dfs.append(df)
            output_df = pd.concat(output_dfs)
            output_df.to_csv(output_path)
