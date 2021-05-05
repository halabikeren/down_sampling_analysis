import typing as t
from dataclasses import dataclass
import os
import re
from io import StringIO
import pandas as pd
from ete3 import Tree
import itertools
import json
from utils import SequenceDataType, TreeReconstructionMethod, BaseTools, SimulationInput
from .program import Program

from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())


@dataclass
class Paml(Program):
    def __init__(self):
        super().__init__()
        self.name = "paml"
        self.program_exe = os.environ["paml"]
        self.cluster_program_exe = os.environ["cluster_paml"]
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
        control_file_path: str = f"{os.getcwd()}/paml.ctl",
        input_tree_path: str = f"{os.getcwd()}/paml_tree.nwk",
    ) -> t.List[str]:
        """
        :param input_path: path to the input of the program
        :param output_path: path to the output of the program
        :param additional_params: additional parameters to run the program with (maps parameter name to value)
        :param parallelize: boolean indicating weather program execution should be parallelized
        :param cluster_data_dir: directory to concat to directory arguments in case of parallelization on the cluster
        :param sequence_data_type: indicates the type of
        :param control_file_path: path in which a control file will be generated
        :param input_tree_path: path in which the input tree for paml will be generated
        :return: a list of strings representing the command
        """
        program_input_path = (
            input_path
            if not parallelize or not os.environ["in_container"]
            else input_path.replace(os.environ["container_data_dir"], cluster_data_dir)
        )
        program_output_path = (
            output_path
            if not parallelize or not os.environ["in_container"]
            else output_path.replace(os.environ["container_data_dir"], cluster_data_dir)
        )
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
        # shorted file paths for paml (gives up upon receiving file paths > 120 chars)
        shared_dir = os.path.commonprefix(
            [
                program_input_path,
                input_tree_path,
                program_output_path,
                control_file_path,
            ]
        )
        program_input_path = program_input_path.replace(shared_dir, "./").replace(
            "//", "/"
        )
        input_tree_path = input_tree_path.replace(shared_dir, "./").replace("//", "/")
        program_output_path = program_output_path.replace(shared_dir, "./").replace(
            "//", "/"
        )

        self.set_control_file(
            program_input_path,
            input_tree_path,
            program_output_path,
            control_file_path,
            additional_params,
            sequence_data_type=sequence_data_type,
        )
        control_file_path = control_file_path.replace(shared_dir, "./").replace(
            "//", "/"
        )
        return [f"cd {shared_dir}", f"{self.program_exe} {control_file_path}"]

    @staticmethod
    def set_control_file(
        input_aln_path: str,
        input_tree_path: str,
        output_path: str,
        control_file_path: str,
        additional_params: t.Optional[t.Dict[str, str]] = None,
        sequence_data_type: SequenceDataType = SequenceDataType.CODON,
    ):
        """
        :param input_aln_path: alignment path
        :param input_tree_path: tree path
        :param output_path: program output path
        :param control_file_path: path in which the control file will be generated
        :param additional_params: additional parameters for PAML's control file
        :param sequence_data_type: sequence data type (needs to be specified in the input file)
        :return: none. writes the data to the control file
        """
        control_file_parameters = {
            "seqfile": input_aln_path,
            "treefile": input_tree_path,
            "outfile": output_path,
            "noisy": 0,  # 0,1,2,3,9: how much rubbish on the screen
            "verbose": 0,  # 1: detailed output, 0: concise output
            "runmode": 0,  # 0: user tree;  1: semi-automatic;  2: automatic; 3: StepwiseAddition; (4,5):PerturbationNNI
            "seqtype": 1
            if sequence_data_type == SequenceDataType.CODON
            else 2,  # 1:codons; 2:AAs; 3:codons-->AAs
            "CodonFreq": 2,  # 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
            "clock": 0,  # 0: no clock, unrooted tree, 1: clock, rooted tree
            "model": 0,  # use a single branch category (site model)
            "NSsites": 8,  # M8 model. for more options, see http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf
            "icode": 0,  # 0:standard genetic code; 1:mammalian mt; 2-10:see below
            "fix_kappa": 0,  # 1: kappa fixed, 0: kappa to be estimated
            "kappa": 2,  # initial or fixed kappa
            "fix_omega": 0,  # 1: omega or omega_1 fixed, 0: estimate
            "omega": 1,  # initial or fixed omega, for codons or codon-transltd AAs
            "fix_alpha": 1,  # 0: estimate gamma shape parameter; 1: fix it at alpha
            "alpha": 0.0,  # initial or fixed alpha, 0:infinity (constant rate)
            "Malpha": 0,  # different alphas for genes
            "ncatG": 4,  # # of categories in the dG or AdG models of rates
            "getSE": 0,  # 0: don't want them, 1: want S.E.s of estimates
            "RateAncestor": 0,  # (1/0): rates (alpha>0) or ancestral states (alpha=0)
            "fix_blength": 1,  # 0: ignore, -1: random, 1: initial, 2: fixed
            "method": 0,  # 0: simultaneous; 1: one branch at a time
        }
        if additional_params:
            for field in control_file_parameters:
                if field in additional_params:
                    control_file_parameters[field] = additional_params[field]

        control_file_content = "\n".join(
            [
                f"{field} = {control_file_parameters[field]}"
                for field in control_file_parameters
            ]
        )
        with open(control_file_path, "w") as control_file:
            control_file.write(control_file_content)

    @staticmethod
    def parse_reference_data(input_path: str) -> t.Dict[str, t.Any]:
        """
        :param input_path: path to the reference data
        :return: a dictionary with the parsed reference data
        """
        rates_data_regex = re.compile("(Site\s*Class.*)", re.MULTILINE | re.DOTALL)
        with open(input_path, "r") as input_file:
            rates_data = rates_data_regex.search(input_file.read()).group(1)
        f = StringIO(rates_data)
        rates_df = pd.read_csv(f, sep="\t")
        return rates_df.to_dict()

    @staticmethod
    def parse_positive_selection_analysis(
        data: str, data_regex: re, column_names
    ) -> t.Dict[t.Any, t.Any]:
        """
        :param data: string or the output file content
        :param data_regex: regex to catch the result of the analysis of interest
        :param column_names: names of columns in the generated df
        :return: a dataframe of the parsed result
        """
        if "p(w>1)" not in column_names:
            raise ValueError(
                f"Failed to process data for regex {data_regex.pattern} due to missing column name p(w>1) on {column_names}"
            )
        data_content = data_regex.search(data).group(1)
        data_content = data_content.lstrip()
        f = StringIO(data_content)
        df = pd.read_csv(f, names=column_names, delimiter=r"\s+")
        df["is_significant"] = df["p(w>1)"].apply(lambda x: "*" in str(x))
        df["p(w>1)"] = df["p(w>1)"].apply(lambda x: float(str(x).replace("*", "")))
        return df.to_dict()

    @staticmethod
    def parse_codon_frequencies(content: str) -> t.Dict[str, float]:
        states_frequencies_regex = re.compile(
            "Codon\s*position\s*x\s*base\s*\(3x4\)\s*table\,\s*overall.*?position\s*1\:\s*T\:\s*(\d*\.?\d*)\s*C\:\s*(\d*\.?\d*)\s*A\:\s*(\d*\.?\d*)\s*G\:\s*(\d*\.?\d*)\s*.*?position\s*2\:\s*T\:\s*(\d*\.?\d*)\s*C\:\s*(\d*\.?\d*)\s*A\:\s*(\d*\.?\d*)\s*G\:\s*(\d*\.?\d*)\s*position\s*3\:\s*T\:\s*(\d*\.?\d*)\s*C\:\s*(\d*\.?\d*)\s*A\:\s*(\d*\.?\d*)\s*G\:\s*(\d*\.?\d*)\s*",
            re.MULTILINE | re.DOTALL,
        )
        match = states_frequencies_regex.search(content)
        # extract F3X4 parameters
        f3x4_parameters = {
            "pos1": {
                "T": float(match.group(1)),
                "C": float(match.group(2)),
                "A": float(match.group(3)),
                "G": float(match.group(4)),
            },
            "pos2": {
                "T": float(match.group(5)),
                "C": float(match.group(6)),
                "A": float(match.group(7)),
                "G": float(match.group(8)),
            },
            "pos3": {
                "T": float(match.group(9)),
                "C": float(match.group(10)),
                "A": float(match.group(11)),
                "G": float(match.group(12)),
            },
        }
        # compute codon frequencies according to them
        nucleotides = ["T", "C", "A", "G"]
        stop_codons = ["TAG", "TAA", "TGA"]
        codons = [
            "".join([nuc1, nuc2, nuc3])
            for nuc1 in nucleotides
            for nuc2 in nucleotides
            for nuc3 in nucleotides
        ]
        codon_frequencies = {
            codon: (
                0
                if codon in stop_codons
                else f3x4_parameters["pos1"][codon[0]]
                * f3x4_parameters["pos2"][codon[1]]
                * f3x4_parameters["pos3"][codon[2]]
            )
            for codon in codons
        }
        return codon_frequencies

    @staticmethod
    def parse_output(
        output_path: str, job_output_dir: t.Optional[str] = None
    ) -> t.Dict[str, t.Any]:
        """
        the parser is currently compatible only with site-models
        :param output_path:
        :param job_output_dir:
        :return:
        """
        result = super(Paml, Paml).parse_output(
            output_path=output_path, job_output_dir=job_output_dir
        )
        with open(output_path, "r") as outfile:
            content = outfile.read()

        neb_positive_selection_regex = re.compile(
            ".*Naive Empirical Bayes \(NEB\) analysis.*?Pr\(w>1\).*?for w\n\n(.*?)\n\n",
            re.MULTILINE | re.DOTALL,
        )
        result[
            "NEB_positive_selection_analysis"
        ] = Paml.parse_positive_selection_analysis(
            content,
            neb_positive_selection_regex,
            column_names=["position", "sequence", "p(w>1)", "mean(w)"],
        )
        try:
            beb_positive_selection_regex = re.compile(
                ".*Bayes Empirical Bayes \(BEB\) analysis.*?Pr\(w>1\).*?for w\n\n(.*?)\n\n",
                re.MULTILINE | re.DOTALL,
            )
            result[
                "BEB_positive_selection_analysis"
            ] = Paml.parse_positive_selection_analysis(
                content,
                beb_positive_selection_regex,
                column_names=[
                    "position",
                    "sequence",
                    "p(w>1)",
                    "mean(w)",
                    "sign",
                    "standard_error",
                ],
            )
        except Exception as e:  # no bayes empirical bayes analysis is available
            pass

        inferred_w_regex = re.compile(
            "MLEs of dN/dS \(w\) for site classes.*?\n\n(.*?)\n\n",
            re.MULTILINE | re.DOTALL,
        )
        inferred_w_content = inferred_w_regex.search(content).group(1)
        [props, ws] = [res.split("  ")[1:] for res in inferred_w_content.split("\n")]
        result["selection_parameters"] = {
            cat: {"prop": float(props[cat]), "w": float(ws[cat])}
            for cat in range(len(ws))
        }

        duration_regex = re.compile("Time used:\s*(\d*\:\d*)", re.MULTILINE | re.DOTALL)
        result["duration(minutes)"] = float(
            (duration_regex.search(content).group(1).replace(":", "."))
        )
        kappa_regex = re.compile("kappa\s*\(ts\/tv\)\s= \s(\d*\.?\d*)")
        result["kappa"] = float(kappa_regex.search(content).group(1))
        tree_length_regex = re.compile("tree\s*length\s*=\s*(\d*\.?\d*)")
        result["tree_length"] = float(tree_length_regex.search(content).group(1))
        tree_str_regex = re.compile("(\(.*?\)\;)")
        tree_str = [match for match in tree_str_regex.finditer(content)][0].group(1)
        tree = Tree(tree_str)
        tree_path = f"{os.path.dirname(output_path)}/tree.nwk"
        tree.write(outfile=tree_path, format=5)
        result["tree_path"] = tree_path
        result["states_frequencies"] = Paml.parse_codon_frequencies(content=content)
        return result

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
        :return:
        """
        substitution_model_params_values = {
            "kappa": program_output["kappa"],
            "selection_parameters": program_output["selection_parameters"],
        }
        simulation_input_parameters = {
            "substitution_model": "",
            "substitution_model_params": substitution_model_params_values,
            "states_frequencies": program_output["states_frequencies"],
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
