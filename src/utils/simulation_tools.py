import pathlib
import subprocess
import typing as t
from dataclasses import dataclass

from ete3 import Tree

from .base_tools import BaseTools
from .types import SequenceDataType
from .simulation_input import SimulationInput
from enum import Enum
import socket
import re
import json
import os

from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())

import logging

logger = logging.getLogger(__name__)

nucleotide_substitution_models = [
                "JC",
                "F81",
                "K80",
                "HKY",
                "TrNef",
                "TrN",
                "K81",
                "K81uf",
                "TIMef",
                "TIM",
                "TVMef",
                "TVM",
                "SYM",
                "GTR",
                "F84ef",
                "F84",
            ]
protein_substitution_models = [
                "Poisson",
                "JTT",
                "JTT-dcmut",
                "Dayhoff",
                "Dayhoff-dcmut",
                "WAG",
                "mtMAM",
                "mtART",
                "mtREV",
                "rtREV",
                "cpREV",
                "Vt",
                "Blosum",
                "LG",
                "HIVb",
                "HIVw",
                "USER",
            ]
nucleotides = ["T", "C", "A", "G"]
amino_acids = [
                "A",
                "R",
                "N",
                "D",
                "C",
                "Q",
                "E",
                "G",
                "H",
                "I",
                "L",
                "K",
                "M",
                "F",
                "P",
                "S",
                "T",
                "W",
                "Y",
                "V",
            ]
codons = [
                "TTT",
                "TTC",
                "TTA",
                "TTG",
                "TCT",
                "TCC",
                "TCA",
                "TCG",
                "TAT",
                "TAC",
                "TAA",
                "TAG",
                "TGT",
                "TGC",
                "TGA",
                "TGG",
                "CTT",
                "CTC",
                "CTA",
                "CTG",
                "CCT",
                "CCC",
                "CCA",
                "CCG",
                "CAT",
                "CAC",
                "CAA",
                "CAG",
                "CGT",
                "CGC",
                "CGA",
                "CGG",
                "ATT",
                "ATC",
                "ATA",
                "ATG",
                "ACT",
                "ACC",
                "ACA",
                "ACG",
                "AAT",
                "AAC",
                "AAA",
                "AAG",
                "AGT",
                "AGC",
                "AGA",
                "AGG",
                "GTT",
                "GTC",
                "GTA",
                "GTG",
                "GCT",
                "GCC",
                "GCA",
                "GCG",
                "GAT",
                "GAC",
                "GAA",
                "GAG",
                "GGT",
                "GGC",
                "GGA",
                "GGG",
            ]

@dataclass
class SimulationTools:
    @staticmethod
    def check_model_legality(simulation_input: SimulationInput):
        if (
            simulation_input.sequence_data_type == SequenceDataType.NUC
            and simulation_input.substitution_model
            not in nucleotide_substitution_models
        ):
            logger.error(
                f"provided model {simulation_input.substitution_model} is not in nucleotide substitution model set {nucleotide_substitution_models}"
            )
            raise ValueError(
                f"provided model {simulation_input.substitution_model} is not in nucleotide substitution model set {nucleotide_substitution_models}"
            )

        elif (
            simulation_input.sequence_data_type == SequenceDataType.AA
            and simulation_input.substitution_model
            not in protein_substitution_models
        ):
            logger.error(
                f"provided model {simulation_input.substitution_model} is not in protein substitution model set {protein_substitution_models}"
            )
            raise ValueError(
                f"provided model {simulation_input.substitution_model} is not in protein substitution model set {protein_substitution_models}"
            )
        elif (
            simulation_input.sequence_data_type == SequenceDataType.AA
            and simulation_input.substitution_model == "USER"
            and (
                "parameters_file" not in simulation_input.substitution_model_params
                or not os.path.exists(
                    simulation_input.substitution_model_params["parameters_file"]
                )
            )
        ):
            logger.error(
                "chosen aa model is USER but parameters file is either not provided or does not exist"
            )
            raise ValueError(
                "chosen aa model  is USER but parameters file is either not provided or does not exist"
            )

    @staticmethod  # to do: add parsing of codon and aa models
    def get_substitution_parameter(
        state_1: str, state_2: str, parameters: t.Dict[str, t.Union[float,str]]
    ) -> float:
        return parameters.get(f"{state_1.upper()},{state_2.upper()}") or parameters.get(f"{state_2.upper()},{state_1.upper()}")


    @staticmethod
    def parse_simulation_substitution_parameters(
        simulation_input: SimulationInput,
    ) -> str:
        """
        :param simulation_input: SimulationInput instance
        :return: string representing the parsed parameters in indelible format
        """
        substitution_params_str = ""
        if simulation_input.substitution_model in [
            "K80",
            "HKY",
            "TrNef",
            "TrN",
            "TIMef",
            "TIM",
            "SYM",
            "GTR",
        ]:  # add a=T<->C
            substitution_params_str += str(
                SimulationTools.get_substitution_parameter(
                    state_1="T",
                    state_2="C",
                    parameters=simulation_input.substitution_model_params,
                )
            )
        if simulation_input.substitution_model in ["TrNef", "TrN"]:  # add f=A<->G
            substitution_params_str += f" {SimulationTools.get_substitution_parameter(state_1='A', state_2='G', parameters=simulation_input.substitution_model_params)}"
        if simulation_input.substitution_model in [
            "TIMef",
            "TIM",
            "SYM",
            "GTR",
            "K81",
            "K81uf",
            "TVMef",
            "TVM",
        ]:  # add b=A<->T, c=G<->T
            substitution_params_str += f" {SimulationTools.get_substitution_parameter(state_1='A', state_2='T', parameters=simulation_input.substitution_model_params)}"
            substitution_params_str += f" {SimulationTools.get_substitution_parameter(state_1='G', state_2='T', parameters=simulation_input.substitution_model_params)}"
        if simulation_input.substitution_model in [
            "TVMef",
            "TVM",
            "SYM",
            "GTR",
        ]:  # add d=A<->C, e=G<->C
            substitution_params_str += f" {SimulationTools.get_substitution_parameter(state_1='A', state_2='C', parameters=simulation_input.substitution_model_params)}"
            substitution_params_str += f" {SimulationTools.get_substitution_parameter(state_1='G', state_2='C', parameters=simulation_input.substitution_model_params)}"
        if simulation_input.substitution_model in [
            "F84ef",
            "F84",
        ]:  # add k s.t. a=(1+k/Y), f=(1+k/R)
            a = SimulationTools.get_substitution_parameter(
                state_1="T",
                state_2="C",
                parameters=simulation_input.substitution_model_params,
            )
            k = (
                a
                * (
                    simulation_input.states_frequencies["T"]
                    + simulation_input.states_frequencies["C"]
                )
                - 1
            )
            substitution_params_str += f" {k}"
        elif simulation_input.substitution_model == "":  # codon model
            substitution_params_str += (
                f"{simulation_input.substitution_model_params['kappa']}\n"
            )
            ncat = len(
                simulation_input.substitution_model_params[
                    "selection_parameters"
                ].keys()
            )
            selection_parameters_str = (
                " ".join(
                    [
                        str(
                            simulation_input.substitution_model_params[
                                "selection_parameters"
                            ][str(cat)]["prop"]
                        )
                        for cat in range(ncat - 1)
                    ]
                )
                + "\n"
                + " ".join(
                    [
                        str(
                            simulation_input.substitution_model_params[
                                "selection_parameters"
                            ][str(cat)]["w"]
                        )
                        for cat in range(ncat)
                    ]
                )
            )
            substitution_params_str += selection_parameters_str
        return substitution_params_str

    @staticmethod
    def parse_states_frequencies(simulation_input: SimulationInput) -> str:
        if simulation_input.sequence_data_type == SequenceDataType.NUC:
            ordered_states = nucleotides
        elif simulation_input.sequence_data_type == SequenceDataType.AA:
            ordered_states = amino_acids
        elif simulation_input.sequence_data_type == SequenceDataType.CODON:
            ordered_states = codons
        state_to_frequency = {
            state: (
                simulation_input.states_frequencies[state]
                if state in simulation_input.states_frequencies
                else 1 / len(ordered_states)
            )
            for state in ordered_states
        }
        sum_of_frequencies = sum(list(state_to_frequency.values()))
        scaling_factor = 1 / sum_of_frequencies
        for state in state_to_frequency:
            state_to_frequency[state] *= scaling_factor
        states_frequencies_str = " ".join(
            [str(state_to_frequency[state]) for state in ordered_states]
        )
        return states_frequencies_str

    @staticmethod
    def parse_simulation_tree(simulation_input: SimulationInput) -> str:
        simulation_tree_str = ""
        if not simulation_input.tree_random:
            with open(simulation_input.simulation_tree_path, "r") as tree_file:
                simulation_tree = Tree(tree_file.read())
            if simulation_input.tree_length:
                BaseTools.scale_tree(
                    tree=simulation_tree, required_size=simulation_input.tree_length
                )
            for node in simulation_tree.traverse():
                node.dist = round(node.dist, 4)
            simulation_tree_str = simulation_tree.write(format=5)
        tree_settings = simulation_tree_str
        if simulation_input.tree_random:
            tree_settings = (
                tree_settings
                + f"""
            {'[rooted]' if simulation_input.tree_rooted else '[unrooted]'} {simulation_input.ntaxa} {simulation_input.birth_rate} {simulation_input.death_rate} {simulation_input.sample_rate}  {simulation_input.mutation_rate}
            [treelength] {simulation_input.tree_length}"""
            )
        return tree_settings

    @staticmethod
    def write_control_file(simulation_input: SimulationInput, output_path: str):
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        SimulationTools.check_model_legality(simulation_input)
        rates_str = (
            f"[rates]     {simulation_input.pinv} {simulation_input.alpha} {simulation_input.ngamcat}"
            if simulation_input.sequence_data_type != SequenceDataType.CODON
            else ""
        )
        control_content = f"""[TYPE] {'NUCLEOTIDE' if simulation_input.sequence_data_type == SequenceDataType.NUC else ('CODON' if simulation_input.sequence_data_type == SequenceDataType.CODON else 'AMINOACID')} 1
                                        [SETTINGS]
                                            [ancestralprint]    NEW
                                            [fastaextension]    fasta
                                            [output]          	FASTA
                                            [fileperrep]      	TRUE
                                            [printrates]        TRUE

                                        [MODEL]    model
                                          [submodel]  {simulation_input.substitution_model} {SimulationTools.parse_simulation_substitution_parameters(simulation_input=simulation_input)}
                                          [statefreq] {SimulationTools.parse_states_frequencies(simulation_input=simulation_input)}
                                          {rates_str}

                                        [TREE] tree {SimulationTools.parse_simulation_tree(simulation_input=simulation_input)}

                                        [PARTITIONS] basic_partition [tree model {simulation_input.seq_len}]

                                        [EVOLVE]     basic_partition  1  seq_data
                    """
        with open(output_path, "w") as control_file:
            control_file.write(control_content)

    @staticmethod
    def write_pipeline_input(simulation_input: SimulationInput, output_path: str):
        nwk_tree_pattern = re.compile("TREE STRING.*?(\(.*?\;)", re.MULTILINE | re.DOTALL)
        pipeline_json_input = json.loads(simulation_input.json())
        pipeline_json_input["pipeline_dir"] = f"{os.getcwd()}/pipeline_dir/"
        pipeline_json_input["cluster_data_dir"] = os.getcwd()
        pipeline_json_input["unaligned_sequence_data_path"] = f"{os.getcwd()}/seq_data_1.fasta"
        if simulation_input.use_simulated_alignment:
            pipeline_json_input[
                "aligned_sequence_data_path"
            ] = f"{os.getcwd()}/seq_data_TRUE_1.fasta"
        if simulation_input.use_simulated_tree:
            tree_path = f"{os.getcwd()}/simulated_tree.nwk"
            if simulation_input.tree_random:
                with open(f"{os.getcwd()}/trees.txt", "r") as trees_file:
                    tree_str = nwk_tree_pattern.search(trees_file.read()).group(1)
                    tree = Tree(tree_str, format=1)
                    tree.write(outfile=tree_path, format=1)
            else:
                with open(simulation_input.simulation_tree_path, "r") as tree_file:
                    simulation_tree = Tree(tree_file.read())
                simulation_tree.write(outfile=tree_path, format=1)
            pipeline_json_input["tree_path"] = tree_path
        from programs import ProgramName

        pipeline_json_input["reference_data_paths"] = {
            ProgramName.RATE4SITE.value: {
                "per_position_reference": f"{os.getcwd()}/seq_data_RATES.txt"
            },
            ProgramName.BUSTED.value: {
                "parameters_reference": f"{os.getcwd()}/control.txt"
            },
            ProgramName.PAML.value: {
                "parameters_reference": f"{os.getcwd()}/control.txt",
                "per_position_reference": f"{os.getcwd()}/seq_data_RATES.txt",
            },
        }
        with open(output_path, "w") as json_file:
            json.dump(pipeline_json_input, json_file)

    @staticmethod
    def simulate(simulation_input: SimulationInput) -> t.List[str]:
        number_of_repeats = simulation_input.nrep
        simulation_output_dir = simulation_input.simulations_output_dir
        os.makedirs(simulation_output_dir, exist_ok=True)
        pipeline_input_paths = []
        for n in range(1, number_of_repeats + 1):
            # create control file
            output_dir = f"{simulation_output_dir}/rep_{n}/"
            control_file_path = f"{output_dir}/control.txt"
            SimulationTools.write_control_file(
                simulation_input=simulation_input, output_path=control_file_path
            )
            # simulate
            os.chdir(output_dir)
            cmd = f"{os.environ['cluster_indelible'] if 'tau' in socket.gethostname() else os.environ['indelible']}"
            res = os.system(cmd)
            if res != 0 or not os.path.exists(
                f"{os.getcwd()}/seq_data_1.fasta"
            ):
                if not os.path.exists(control_file_path):
                    logger.error(f"control file does not exists at {output_dir}. The only available content is {os.listdir(output_dir)}")
                raise RuntimeError(
                    f"INDELible failed to properly execute and provide an output file with error {res} and command {cmd}"
                )
            # prepare pipeline input
            json_path = f"{os.getcwd()}/input.json"
            SimulationTools.write_pipeline_input(
                simulation_input=simulation_input, output_path=json_path
            )
            pipeline_input_paths.append(json_path)

        return pipeline_input_paths
