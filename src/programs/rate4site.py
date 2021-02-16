import typing as t
from dataclasses import dataclass
import os
import re
from io import StringIO
import pandas as pd
from .program import Program

from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())


import logging

logger = logging.getLogger(__name__)


@dataclass
class Rate4Site(Program):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.name = "rate4site"
        self.program_exe = os.environ["rate4site"]
        self.cluster_program_exe = os.environ["cluster_rate4site"]
        self.input_param_name = "-s"
        self.output_param_name = "-o"
        self.module_to_load = "Rate4Site/Rate4Site-3.0"

    @staticmethod
    def parse_rates(rates_content: str) -> t.Dict[t.Any, t.Any]:
        """
        :param rates_content: a string given from rate4site output file
        :return: a parsed rates content in the form of a json
        """
        rates_content = rates_content.lstrip()
        rates_content = re.sub(",\s*", ",", rates_content)
        rates_content = re.sub("\[\s*", "[", rates_content)
        rates_content = re.sub("\s*\]", "]", rates_content)
        f = StringIO(rates_content)
        rates_df = pd.read_csv(
            f,
            names=["position", "sequence", "rate", "qq_interval", "std", "msa_data"],
            delimiter=r"\s+",
        )
        return rates_df.to_dict()

    @staticmethod
    def parse_output(
        output_path: str, job_output_dir: t.Optional[str] = None
    ) -> t.Dict[str, t.Any]:
        """
        :param output_path
        :param job_output_dir
        :return: None. parses the output file into a json form and saves it into self.result
        """
        result = super(Rate4Site, Rate4Site).parse_output(output_path=output_path, job_output_dir=job_output_dir)
        with open(output_path, "r") as output_file:
            output_content = output_file.read()
        output_regex = re.compile(
            "#POS\s*SEQ\s*SCORE\s*QQ-INTERVAL\s*STD\s*MSA DATA.*?The alpha parameter (.\d*\.?\d*).*?LL=(-?\d*\.?\d*)(.*?)#Average",
            re.MULTILINE | re.DOTALL,
        )
        output_match = output_regex.search(output_content)
        result["alpha"] = float(output_match.group(1))
        result["log_likelihood"] = float(output_match.group(2))
        result["rate_by_position"] = Rate4Site.parse_rates(output_match.group(3))
        return result

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
        rates_df.rename({"Site": "position", "Rate": "rate"})
        reference_data = {"rate_by_position": rates_df.to_dict()}
        return reference_data

    @staticmethod
    def get_accuracy(reference_data: t.Dict[str, t.Any], test_data: t.Dict[str, t.Any]) -> pd.Series:
        """
        :param reference_data: reference data to compute results by reference to
        :param test_data: test data to compare to the reference data
        :return: the output of pd series with indices as the members for which accuracy it assessed (be it positions in a sequence of sequences) and the values are the accuracy values computed for them
        """
        reference_df = pd.DataFrame.from_dict(reference_data["rate_by_position"])
        test_df = pd.DataFrame.from_dict(test_data["rate_by_position"])
        test_positions = list(test_df["position"].values)
        reference_positions = list(reference_df["position"].values)
        if len(test_positions) < len(reference_positions):
            logger.error(f"Number of positions in test data is {len(test_positions)} and is inconsistent with the number of positions in the reference data {len(reference_positions)}")
            raise ValueError(f"Number of positions in test data is {len(test_positions)} and is inconsistent with the number of positions in the reference data {len(reference_positions)}")
        reference_df = reference_df.loc[reference_df["position"].isin(test_positions)]
        absolute_error = abs(reference_df["rate"]-test_df["rate"])
        denominator = abs(reference_df["rate"]) + abs(test_df["rate"])
        relative_error = absolute_error/denominator
        penalized_error_by_std = relative_error * (abs(reference_df["std"]-test_df["std"])/test_df["std"])  # will punish error with low test std more than one without
        penalized_accuracy_by_std = 1-penalized_error_by_std
        return penalized_accuracy_by_std
