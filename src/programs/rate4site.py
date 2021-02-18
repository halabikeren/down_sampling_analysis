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
        # TO DO: parse non normalized rates (required for comparison to simulated rates)
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
        denormalized_rates_by_position_path = f"{job_output_dir}/r4sOrig.res"
        with open(denormalized_rates_by_position_path, "r") as output_file:
            output_content = output_file.read()
        result["denormalized_rate_by_position"] = Rate4Site.parse_rates(output_regex.search(output_content).group(3))
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
        rates_df.rename(columns={"Site": "position", "Rate": "rate"}, inplace=True)
        rates_df["std"] = rates_df["rate"].std()
        reference_data = {"denormalized_rate_by_position": rates_df.to_dict()}
        return reference_data

    @staticmethod
    def get_accuracy(reference_data: t.Dict[str, t.Any], test_data: t.Dict[str, t.Any], use_normalized_rates: bool = False, penalize_by_std: bool = False) -> pd.Series:
        """
        :param reference_data: reference data to compute results by reference to
        :param test_data: test data to compare to the reference data
        :param use_normalized_rates: indicates weather normalized rates should be used for accuracy computation or denormalized rates
        :param penalize_by_std: boolean indicating weather to penalize by std or not
        :return: the output of pd series with indices as the members for which accuracy it assessed (be it positions in a sequence of sequences) and the values are the accuracy values computed for them
        """
        rate_type = "rate_by_position"
        if not use_normalized_rates:
            rate_type = "denormalized_rate_by_position"
        reference_df = pd.DataFrame.from_dict(reference_data[rate_type])
        test_df = pd.DataFrame.from_dict(test_data[rate_type])
        test_positions = list(test_df["position"].values)
        reference_positions = list(reference_df["position"].values)
        if len(test_positions) < len(reference_positions):
            logger.error(f"Number of positions in test data is {len(test_positions)} and is inconsistent with the number of positions in the reference data {len(reference_positions)}")
            raise ValueError(f"Number of positions in test data is {len(test_positions)} and is inconsistent with the number of positions in the reference data {len(reference_positions)}")
        reference_df = reference_df.loc[reference_df["position"].isin(test_positions)]
        absolute_error = abs(reference_df["rate"]-test_df["rate"])
        denominator = 2 if use_normalized_rates else (max(list(reference_df["rate"])+list(test_df["rate"]))-min(list(reference_df["rate"])+list(test_df["rate"])))  # in the case of normalized rates: corresponds to the range of rates [-1,1] which is also the maximal value of difference between rates
        relative_error = absolute_error/denominator
        denominator = max(list(reference_df["std"])+list(test_df["std"]))-min(list(reference_df["std"])+list(test_df["std"]))
        relative_accuracy = 1-relative_error
        if penalize_by_std:
            std_normalizer = (abs(reference_df["std"] - test_df["std"])) / denominator
            penalized_error_by_std = relative_error * std_normalizer  # will punish error with low test std more than one without
            penalized_accuracy_by_std = 1 - penalized_error_by_std
            return penalized_accuracy_by_std
        return relative_accuracy  # for now, use relative accuracy (17.2.2021)

    @staticmethod
    def get_result(data: t.Dict[str, t.Any], use_normalized_rates: bool = False) -> pd.Series:
        """
        :param data: dictionary mapping results
        :param use_normalized_rates: indicates weather normalized rates should be used for accuracy computation or denormalized rates
        :return: the relevant data to compute accuracy for
        """
        df = data["denormalized_rate_by_position"]
        if use_normalized_rates:
            df = data["rates_by_position"]
        return df["rate"]