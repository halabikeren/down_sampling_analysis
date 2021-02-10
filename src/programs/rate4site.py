import typing as t
from dataclasses import dataclass
import os
import re
from io import StringIO
import pandas as pd
from .program import Program

from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())


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
    def parse_rates(rates_content: str) -> pd.DataFrame:
        """
        :param rates_content: a string given from rate4site output file
        :return: a parsed rates content in the form of a json
        """
        rates_content = rates_content.lstrip()
        rates_content = re.sub(",\s*", ",", rates_content)
        f = StringIO(rates_content)
        rates_df = pd.read_csv(
            f,
            names=["position", "sequence", "rate", "qq_interval", "std", "msa_data"],
            delimiter=r"\s+",
        )
        return rates_df

    @staticmethod
    def parse_output(
        output_path: str, aux_dir: t.Optional[str] = None
    ) -> t.Dict[str, t.Any]:
        """
        :return: None. parses the output file into a json form and saves it into self.result
        """
        with open(output_path, "r") as output_file:
            output_content = output_file.read()
        output_regex = re.compile(
            "#POS\s*SEQ\s*SCORE\s*QQ-INTERVAL\s*STD\s*MSA DATA.*?The alpha parameter (.\d*\.?\d*).*?LL=(-?\d*\.?\d*)(.*?)#Average",
            re.MULTILINE | re.DOTALL,
        )
        output_match = output_regex.search(output_content)
        result = dict()
        result["alpha"] = float(output_match.group(1))
        result["log_likelihood"] = float(output_match.group(2))
        result["rate_by_position"] = Rate4Site.parse_rates(output_match.group(3))
        if aux_dir:
            job_path = [jpath for jpath in os.listdir(aux_dir) if ".OU" in jpath][0]
            with open(job_path, "r") as outfile:
                job_content = outfile.readlines()
                start_time = float(job_content[0])
                end_time = float(job_content[-1])
                result["duration"] = end_time - start_time
        return result
