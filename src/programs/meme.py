import re
import typing as t
from dataclasses import dataclass
import os
from ete3 import Tree
import json
import pandas as pd
from utils import SequenceDataType, TreeReconstructionMethod, BaseTools, SimulationInput
from .program import Program

from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())

import logging

logger = logging.getLogger(__name__)


@dataclass
class Meme(Program):
    def __init__(self):
        super().__init__()
        self.name = "hyphy"
        self.program_exe = os.environ["hyphy"]
        self.cluster_program_exe = os.environ["cluster_hyphy"]
        self.tree_reconstruction_method: TreeReconstructionMethod = (
            TreeReconstructionMethod.FASTTREE
        )
        self.tree_reconstruction_prams: t.Optional[t.Dict[str, str]] = None

    def set_command(
        self,
        input_path: str,
        output_path: str,
        additional_params: t.Optional[t.Dict[str, str]],
        parallelize: bool,
        cluster_data_dir: str,
        sequence_data_type: SequenceDataType = SequenceDataType.CODON,
        input_tree_path: t.Optional[str] = None,
    ) -> t.List[str]:
        """
        :param input_path: path to the input of the program
        :param output_path: path to the output of the program
        :param additional_params: additional parameters to run the program with (maps parameter name to value)
        :param parallelize: boolean indicating weather program execution should be parallelized
        :param cluster_data_dir: directory to concat to directory arguments in case of parallelization on the cluster
        :param sequence_data_type: indicates the type of
        :param input_tree_path: path in which the input tree for paml will be generated
        :return: a string representing the command
        """
        program_input_path = (
            input_path
            if not parallelize
            else input_path.replace(os.environ["container_data_dir"], cluster_data_dir)
        )
        program_output_path = (
            output_path
            if not parallelize
            else output_path.replace(os.environ["container_data_dir"], cluster_data_dir)
        )
        if not input_tree_path:
            input_tree_path = f"{os.path.dirname(program_input_path)}/hyphy_tree.nwk"
        # program_output_path = program_output_path if os.path.isdir(program_output_path) else os.path.dirname(program_output_path)
        self.set_additional_params(
            additional_params, parallelize, cluster_data_dir, return_as_str=False
        )

        if additional_params and "input_tree_path" in additional_params:
            input_tree_path = additional_params["input_tree_path"]
        BaseTools.build_tree(
            input_path,
            input_tree_path,
            sequence_data_type,
            self.tree_reconstruction_method,
            self.tree_reconstruction_prams,
        )

        cmd = f"printf '1\\n1\\n{program_input_path}\\n{input_tree_path}\\n' | hyphy"
        return [cmd, f"cp -r {program_input_path}.MEME.json {program_output_path}"]

    @staticmethod
    def parse_reference_data(input_path: str) -> t.Dict[str, t.Any]:
        """
        :param input_path: path to the reference data
        :return: a dictionary with the parsed reference data
        """
        parameter_to_regex = {
            "kappa": re.compile("\[submodel\]\s*(\d*\.?\d*)", re.DOTALL),
            "omegas": re.compile("\[submodel\].*?\n.*?\n(.*?)\n", re.DOTALL),
            "omegas_proportions": re.compile("\[submodel\].*?\n(.*?)\n", re.DOTALL),
            "states_frequencies": re.compile("\[statefreq\](.*?)\n"),
            "tree": re.compile("\[TREE\][^\(]*(.*?)\n", re.DOTALL),
        }
        reference_data = dict()
        with open(input_path, "r") as infile:
            content = infile.read()
        reference_data["kappa"] = float(
            parameter_to_regex["kappa"].search(content).group(1)
        )
        reference_data["omegas_to_props"] = dict()
        omegas = [
            float(item)
            for item in parameter_to_regex["omegas_proportions"]
            .search(content)
            .group(1)
            .split("\s*")
        ]
        omegas_proportions = [
            float(item)
            for item in parameter_to_regex["omegas_proportions"]
            .search(content)
            .group(1)
            .split("\s*")
        ]
        prop_last_omega = 1 - sum(omegas_proportions)
        omegas_proportions.append(prop_last_omega)
        for i in range(len(omegas)):
            reference_data["omegas_to_props"][omegas[i]] = omegas_proportions[i]
        reference_data["tree"] = parameter_to_regex["tree"].search().group(1)
        return reference_data

    @staticmethod
    def write_inferred_tree(
        inference_results: t.Dict[str, t.Any],
        output_path: str,
        by_model: str = "unconstrained",
    ):
        """
        :param inference_results: dictionary holding the output json of the program
        :param output_path: path to write the tree with the inferred branch lengths to
        :param by_model: the model by which the inferred branch lengths should be extracted
        :return: none. writes the tree to the given file
        """
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        branch_lengths_data = inference_results["branch attributes"]["0"]
        tree_str = inference_results["input"]["trees"]["0"]
        tree = Tree(f"{tree_str};", format=1)
        tree_length = 0
        for node in tree.traverse():
            if node.name != "":
                required_branch_length = branch_lengths_data[node.name][by_model]
                node.dist = required_branch_length
                tree_length += required_branch_length
        inference_results["tree_length"] = tree_length
        tree.write(outfile=output_path, format=5)

    @staticmethod
    def parse_output(
        output_path: str, job_output_dir: t.Optional[str] = None
    ) -> t.Dict[str, t.Any]:
        """
        the parser is currently compatible only with site-models
        :param output_path: path to json output file of the program
        :param job_output_dir: path from which the job was submitted
        :return:
        """
        if not "BUSTED" in output_path:
            try:
                output_path = [
                    f"{output_path}/{path}"
                    for path in os.listdir(output_path)
                    if "BUSTED" in path
                ][0]
            except:
                logger.error(f"No output file for busted in {output_path}")
                exit(1)
        with open(output_path, "r") as json_file:
            results = json.load(json_file)
        codon_freq_vector = results["fits"][
            "MG94xREV with separate rates for branch sets"
        ]["Equilibrium frequencies"]
        nuc = ["A", "C", "G", "T"]
        stop_codons = ["TAG", "TAA", "TGA"]
        codon_hyphy_order = [
            f"{i}{j}{k}"
            for i in nuc
            for j in nuc
            for k in nuc
            if f"{i}{j}{k}" not in stop_codons
        ]
        results["fits"]["MG94xREV with separate rates for branch sets"][
            "Equilibrium frequencies"
        ] = {
            codon_hyphy_order[i]: codon_freq_vector[i][0]
            for i in range(len(codon_hyphy_order))
        }
        for stop_codon in stop_codons:
            results["fits"]["MG94xREV with separate rates for branch sets"][
                "Equilibrium frequencies"
            ][stop_codon] = 0
        results["tree_path"] = f"{os.path.dirname(output_path)}/busted_tree.nwk"
        Meme.write_inferred_tree(
            inference_results=results, output_path=results["tree_path"]
        )
        return results

    @staticmethod
    def write_output_to_simulation_pipeline_json(
        program_output: t.Dict[str, t.Any],
        output_path: str,
        additional_simulation_parameters: t.Dict[str, t.Any],
    ):
        """
        :param program_output: output of the program to translate to simulation params
        :param output_path: output path for simulation pipeline input json
        :param additional_simulation_parameters:  additional parameters
        :return: none, writes simulation input parameters to json
        """
        transitions_rates = sum(
            [
                program_output["fits"]["Nucleotide GTR"]["Rate Distributions"][
                    f"Substitution rate from nucleotide {transition[0]} to nucleotide {transition[1]}"
                ]
                for transition in [("A", "G"), ("C", "T")]
            ]
        )
        transversions_rates = sum(
            [
                program_output["fits"]["Nucleotide GTR"]["Rate Distributions"][
                    f"Substitution rate from nucleotide {transversion[0]} to nucleotide {transversion[1]}"
                ]
                for transversion in [("A", "C"), ("A", "T"), ("C", "G"), ("G", "T")]
            ]
        )
        kappa = transversions_rates / transitions_rates
        hyphy_selection_parameters = program_output["fits"]["Unconstrained model"][
            "Rate Distributions"
        ]["Test"]
        simulation_selection_parameters = {
            int(cat): {
                "prop": hyphy_selection_parameters[cat]["proportion"],
                "w": hyphy_selection_parameters[cat]["omega"],
            }
            for cat in hyphy_selection_parameters
        }
        if simulation_selection_parameters[2]["prop"] > 0.05:
            logger.info(
                f"{output_path} is a fitting candidate for a simulation study involving positive selection"
            )
        simulation_input_parameters = {
            "substitution_model": "",
            "substitution_model_params": {
                "kappa": kappa,
                "selection_parameters": simulation_selection_parameters,
            },
            "states_frequencies": program_output["fits"][
                "MG94xREV with separate rates for branch sets"
            ]["Equilibrium frequencies"],
            "tree_random": False,
            "tree_length": program_output["tree_length"],
            "simulation_tree_path": program_output["tree_path"],
            "pinv": 0,
            "alpha": 0,
            "ngamcat": 0,
        }
        simulation_input_parameters.update(additional_simulation_parameters)
        simulation_input = BaseTools.jsonable_encoder(
            SimulationInput(**simulation_input_parameters)
        )
        clean_simulation_input = {
            k: v for k, v in simulation_input.items() if v is not None
        }
        with open(output_path, "w") as output_file:
            json.dump(obj=clean_simulation_input, fp=output_file)

    @staticmethod
    def get_error(
        reference_data: t.Dict[str, t.Any],
        test_data: t.Dict[str, t.Any],
        use_relative_error: bool = True,
        use_normalized_rates: bool = False,
        penalize_by_std: bool = False,
    ) -> pd.Series:
        """
        :param reference_data: reference data to compute results by reference to
        :param test_data: test data to compare to the reference data
        :param use_relative_error: indicates weather absolute or relative error should be used
        :param use_normalized_rates: indicates weather normalized rates should be used for error computation or denormalized rates
        :param penalize_by_std: boolean indicating weather to penalize by std or not
        :return: the output of pd series with indices as the members for which error it assessed (be it positions in a sequence of sequences) and the values are the error values computed for them
        """

    pass

    @staticmethod
    def get_result(
        data: t.Dict[str, t.Any], use_normalized_rates: bool = False
    ) -> pd.Series:
        """
        :param data: dictionary mapping results
        :param use_normalized_rates: indicates weather normalized rates should be used for error computation or denormalized rates
        :return: the relevant data to compute error for
        """
        pass
