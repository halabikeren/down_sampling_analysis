import json
import typing as t
from dataclasses import dataclass
import os
import re
import pandas as pd
from .program import Program
from utils import SequenceDataType, SimulationInput, BaseTools
from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())

import logging

logger = logging.getLogger(__name__)


@dataclass
class PhyML(Program):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name = "phyml"
        self.program_exe = os.environ["phyml"]
        self.cluster_program_exe = os.environ["cluster_phyml"]
        self.input_param_name = "-i"
        self.output_param_name = ""
        self.module_to_load = ""

    def set_command(
            self,
            input_path: str,
            output_path: str,
            additional_params: t.Optional[t.Dict[str, str]],
            parallelize: bool,
            cluster_data_dir: str,
            sequence_data_type: SequenceDataType,
            output_dir: str
    ) -> str:
        """
        :param input_path: path to the input of the program
        :param output_path: path to the output of the program
        :param additional_params: additional parameters to run the program with (maps parameter name to value)
        :param parallelize: boolean indicating weather program execution should be parallelized
        :param cluster_data_dir: directory to concat to directory arguments in case of parallelization on the cluster
        :param sequence_data_type: indicates the type of
        :param output_dir: directory in which a control file will be generated
        :return: a string representing the command
        """
        program_input_path = (
            input_path
            if not parallelize
            else input_path.replace(os.environ["container_data_dir"], cluster_data_dir)
        )
        program_output_dir = output_dir if not parallelize else output_dir.replace(os.environ["container_data_dir"],
                                                                                   cluster_data_dir)
        default_model = "WAG" if sequence_data_type == SequenceDataType.AA else "GTR"
        cmd = f"cd {program_output_dir}\nphyml {self.input_param_name} {program_input_path} -d {'aa' if sequence_data_type == SequenceDataType.AA else 'nt'} -m {additional_params['model'] if 'model' in additional_params else default_model} -c {additional_params['ncat'] if 'ncat' in additional_params else 16}"
        return cmd

    @staticmethod
    def parse_phyml_stats(input_path: str) -> t.Dict[str, t.Any]:
        """
        :param input_path: phyml stats file path
        :return: dictionary of the parsed output
        """
        with open(input_path, "r") as input_file:
            input_content = input_file.read()
        stats = dict()
        ntaxa_regex = re.compile("Number of taxa\:\s*(\d*)\n", re.DOTALL)
        stats["ntaxa"] = ntaxa_regex.search(input_content).group(1)
        tree_size_regex = re.compile("Tree size\:\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL)
        stats["tree_length"] = float(tree_size_regex.search(input_content).group(1))
        alpha_param_regex = re.compile("Gamma shape parameter\:\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL)
        stats["alpha"] = float(alpha_param_regex.search(input_content).group(1))
        relative_rates_regex = re.compile("Relative rate in class\s(\d*)\:\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL)
        stats["relative_rates"] = {int(match.group(1)): float(match.group(2)) for match in
                                   relative_rates_regex.finditer(input_content)}
        states_frequencies_regex = re.compile("\s-\sf\((\w)\)=\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL)
        stats["states_frequencies"] = {match.group(1): float(match.group(2)) for match in
                                       states_frequencies_regex.finditer(input_content)}
        substitution_model_regex = re.compile("Model of nucleotides substitution\:\s*(.*?)\n", re.DOTALL)
        stats["substitution_model"] = substitution_model_regex.search(input_content).group(1)
        substitution_rates_regex = re.compile("\s*(\w)\s<->\s(\w)\s*(\d*\.?\d*)", re.MULTILINE | re.DOTALL)
        stats["substitution_model_params"] = {str(match.group(1), match.group(2)): float(match.group(3)) for match in
                                       substitution_rates_regex.finditer(input_content)}
        seed_regex = re.compile("Random seed\:\s*(\d*)", re.MULTILINE | re.DOTALL)
        stats["random_seed"] = float(seed_regex.search(input_content).group(1))
        version_regex = re.compile("Version\:\s*(\d*\.?\d*\.?\d*)", re.MULTILINE | re.DOTALL)
        stats["phyml_version"] = version_regex.search(input_content).group(1)
        duration_regex = re.compile("Time used\:\s*(\d*)h(\d*)m(\d*)s", re.MULTILINE | re.DOTALL)
        stats["duration"] = float(duration_regex.search(input_content).group(1)) * 60 + float(
            duration_regex.search(input_content).group(2)) + float(duration_regex.search(input_content).group(3)) / 60
        return stats

    @staticmethod
    def parse_output(output_path: str, job_output_dir: t.Optional[str] = None) -> t.Dict[str, t.Any]:
        """
        :param output_path: directory holding the output files of the program
        :param job_output_dir: directory holding the output of the job, in case parallelization was chosen
        :return: a dictionary holding the parsed result of the program execution
        """
        output = dict()
        tree_path = stats_path = None
        for path in os.listdir(output_path):
            if "phy_phyml_stats.txt" in path:
                stats_path = f"{output_path}/{path}"
            elif "phy_phyml_tree.txt" in path:
                tree_path = f"{output_path}/{path}"
        with open(tree_path, "r") as tree_file:
            tree_str = tree_file.read()
        output["tree_path"] = tree_path
        output["stats_path"] = stats_path
        output["tree"] = tree_str
        output.update(PhyML.parse_phyml_stats(input_path=stats_path))
        return output

    @staticmethod
    def parse_reference_data(input_path: str) -> t.Dict[str, t.Any]:
        pass # TO DO

    @staticmethod
    def get_error(reference_data: t.Dict[str, t.Any], test_data: t.Dict[str, t.Any], **kwargs) -> pd.Series:
        """
        :param reference_data: reference data to compute results by reference to
        :param test_data: test data to compare to the reference data
        :return: the output of pd series with indices as the members for which error it assessed (be it positions in a sequence of sequences) and the values are the error values computed for them
        """
        pass  # TO DO

    @staticmethod
    def get_result(data: t.Dict[str, t.Any], **kwargs) -> pd.Series:
        """
        :param data: dictionary mapping results
        :return: the relevant data to compute error for
        """
        pass  # TO DO

    @staticmethod
    def write_output_to_simulation_pipeline_json(program_output: t.Dict[str, t.Any], output_path: str, additional_simulation_parameters: t.Dict[str, t.Any]):
        """
        :param program_output: output of the program to translate to simulation params
        :param output_path: output path for simulation pipeline input json
        :param additional_simulation_parameters: additional parameters
        :return:
        """
        simulation_input_parameters = {parameter: program_output[parameter] for parameter in ["substitution_model",
                                                                                              "substitution_model_params",
                                                                                              "states_frequencies",
                                                                                              "tree_length",
                                                                                              "alpha"]}

        simulation_input_parameters.update({"tree_random": False,
                                            "simulation_tree_path": program_output["tree_path"],
                                             "pinv": 0,
                                             "ngamcat": len(program_output["relative_rates"].keys())})
        simulation_input_parameters.update(additional_simulation_parameters)
        simulation_input = BaseTools.jsonable_encoder(SimulationInput(**simulation_input_parameters))
        clean_simulation_input = {k: v for k, v in simulation_input.items() if v is not None}
        with open(output_path, "w") as output_file:
            json.dump(output_file, clean_simulation_input)
