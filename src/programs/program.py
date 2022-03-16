import os
import typing as t
from dataclasses import dataclass
from datetime import datetime
from time import time
import pandas as pd
import re

from utils import Job
import json

from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())

import logging

logger = logging.getLogger(__name__)


@dataclass
class Program:
    name: str
    program_exe: str
    cluster_program_exe: str
    input_param_name: str = ""
    output_param_name: str = ""  # maps parameter name to parameter value

    def __init__(
        self,
    ):  # set to allow inheriting classes to receive additional arguments for setting additional
        pass

    @staticmethod
    def set_additional_params(
        additional_params: t.Optional[t.Dict[str, str]],
        parallelize: bool,
        cluster_data_dir: str,
        return_as_str: bool = True,
    ) -> t.Optional[str]:
        """
        :param additional_params: dictionary mapping names to values f additional parameters for the program
        :param parallelize: indicates weather the execution
        of the program should be parallelized. In such case, any additional parameter whose value is a directory will
        be translated ot an absolute cluster directory
        :param cluster_data_dir the cluster data directory that should convert the container or relative directories of
        paths to absolute ones in the cluster
        :param return_as_str: indicates weather the additional parameters should be returned as a string
        :return: if the addition of parameters is expressed via the command line, a string that should be added to it
        will be returned
        """
        additional_params_str = ""
        if additional_params:
            for field in additional_params:
                if (
                    parallelize
                    and os.path.isdir(additional_params[field])
                    and not additional_params[field].startswith(cluster_data_dir)
                ):
                    additional_params[
                        field
                    ] = f"{cluster_data_dir}/{additional_params[field]}"
                elif (
                    not parallelize
                    and os.path.isdir(additional_params[field])
                    and not additional_params[field].startswith(
                        os.environ["container_data_dir"]
                    )
                ):
                    additional_params[
                        field
                    ] = f"{os.environ['container_data_dir']}/{additional_params[field]}"
                additional_params_str = " ".join(
                    [
                        f"{param_name} {additional_params[param_name]}"
                        for param_name in additional_params
                    ]
                )
            if return_as_str:
                return additional_params_str

    def set_command(
        self,
        input_path: str,
        output_path: str,
        additional_params: t.Optional[t.Dict[str, str]],
        parallelize: bool,
        cluster_data_dir: str,
        **kwargs,
    ) -> t.List[str]:
        """
        :param input_path: path to the input of the program
        :param output_path: path to the output of the program
        :param additional_params: additional parameters to run the program with (maps parameter name to value)
        :param parallelize: boolean indicating weather program execution should be parallelized
        :param cluster_data_dir: directory to concat to directory arguments in case of parallelization on the cluster
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
        input_str = f"{self.input_param_name} {program_input_path}"
        output_str = f"{self.output_param_name} {program_output_path}"
        if additional_params:
            additional_params_str = self.set_additional_params(
                additional_params=additional_params,
                parallelize=parallelize,
                cluster_data_dir=cluster_data_dir,
                return_as_str=True,
            )
        else:
            additional_params_str = ""
        command = f"{self.program_exe if not parallelize else self.cluster_program_exe} {input_str} {output_str} {additional_params_str} "
        return [command]

    def exec(
        self,
        input_path: str,
        output_path: str,
        aux_dir: str,
        additional_params: t.Optional[t.Dict[str, str]] = None,
        parallelize: bool = False,
        cluster_data_dir: t.Optional[str] = None,
        priority: int = 0,
        queue: str = "itaym",
        wait_until_complete: bool = False,
        get_completion_validator: bool = True,
    ) -> t.Union[float, str]:
        """
        :param input_path: path to alignment file
        :param output_path: path in which the program should write its output
        :param additional_params: additional parameters unique to the program
        :param parallelize: boolean indicating weather execution of the program should be parallelized in the cluster or not
        :param cluster_data_dir: cluster directory that is mounted to the container data directory. Must be provided with parallleize is True
        :param aux_dir: directory in which auxiliary files should be generated by the job submission process
        :param priority: priority of the jobs
        :param queue: queue to submit the jobs to
        :param wait_until_complete: indicator weather the main program should wait until completion of all jobs (recommended: True)
        :param get_completion_validator: boolean indicating weather a validator file should be generated upon job completion (recommended: True)
        :return: either the duration of the command in minutes, if no parallelization was selected, or the path to the touch file that is used for validation of job completion in case of parallelization
        """
        additional_args = dict()
        from .paml import Paml
        from .busted import Busted
        from.phyml import PhyML

        if type(self) in [Paml, Busted]:
            additional_args["input_tree_path"] = re.sub(
                "\.fas[^.]*", "_tree.nwk", input_path
            )
        if type(self) is PhyML:
            additional_args["output_dir"] = os.path.dirname(output_path)
        if type(self) is Paml:
            additional_args["control_file_path"] = re.sub(
                "\.fas[^.]*", "_paml.ctl", input_path
            )
        command = self.set_command(
            input_path=input_path,
            output_path=output_path,
            additional_params=additional_params,
            parallelize=parallelize,
            cluster_data_dir=cluster_data_dir,
            **additional_args,
        )
        os.makedirs(aux_dir, exist_ok=True)

        if not parallelize:
            start_time = time()
            if type(self) is not Paml:
                os.chdir(
                    aux_dir
                )  # move to aux dir as rate4site generates extra files in current running directory
            for cmd in command:
                if "cd " in cmd:
                    os.chdir(cmd.replace("cd ", ""))
                else:
                    res = os.system(
                        f"{cmd} > /dev/null 2>&1"
                    )  # for some reason, rate4 prints some logs into the stderr,
                    # making the typical test (raise error i=f stderr > 0) invalid in this case
                    if res != 0:
                        raise RuntimeError(f"command {cmd} failed to execute.")
            end_time = time()
            return (end_time - start_time) / 60
        else:
            commands = (
                [
                    f"cd {aux_dir.replace(os.environ['container_data_dir'], cluster_data_dir)}",
                    """timestamp() {
                      date +"%T" # current time
                    }
                    timestamp""",
                ]
                + command
                + ["timestamp"]
            )

            job = Job(
                name=self.name,
                sh_dir=aux_dir,
                output_dir=aux_dir,
                commands=commands,
                priority=priority,
                queue=queue,
            )
            completion_validator = job.submit(
                wait_until_complete=wait_until_complete,
                get_completion_validator=get_completion_validator,
            )
            return completion_validator

    @staticmethod
    def parse_output(
        output_path: str, job_output_dir: t.Optional[str] = None
    ) -> t.Dict[str, t.Any]:
        """
        :param output_path: path holding the output of the program
        :param job_output_dir: directory holding the output of the job, in case parallelization was chosen
        :return: a dictionary holding the parsed result of the program execution
        """
        if not os.path.exists(output_path):
            raise ValueError(f"output path {output_path} does not exist")
        with open(output_path, "r") as out:
            output_content = out.read()

        result = {"raw_output": output_content}

        if job_output_dir:
            paths_by_time = sorted(
                [
                    f"{job_output_dir}/{path}"
                    for path in os.listdir(job_output_dir)
                    if ".OU" in path
                ],
                key=os.path.getmtime,
            )
            if len(paths_by_time) > 0:
                job_output_path = paths_by_time[0]
                with open(job_output_path, "r") as job_output_file:
                    content = job_output_file.readlines()
                try:
                    start_time = datetime.strptime(content[0].rstrip(), "%H:%M:%S")
                    end_time = datetime.strptime(content[-1].rstrip(), "%H:%M:%S")
                    duration = (end_time - start_time).total_seconds() / 60
                    result["duration(minutes)"] = duration
                except:
                    pass

        return result

    @staticmethod
    def parse_reference_data(input_paths: t.Dict[str, str]) -> t.Dict[str, t.Any]:
        """
        :param input_paths: paths to the reference data
        :return: a dictionary with the parsed reference data
        """
        result = dict()
        result["raw_output"] = dict()
        if "parameters_reference" in input_paths:
            input_path = input_paths["parameters_reference"]
            if os.path.exists(input_path):
                with open(input_path, "r") as input_file:
                    result["raw_output"]["parameters"] = input_file.read()
        if "per_position_reference" in input_paths:
            input_path = input_paths["per_position_reference"]
            if os.path.exists(input_path):
                with open(input_path, "r") as input_file:
                    result["raw_output"]["position_wise"] = input_file.read()
        return result

    @staticmethod
    def write_result_json(
        input_path: str, job_output_dir: t.Optional[str], output_path: str
    ):
        """
        :param input_path: path to file with the raw program results
        :param job_output_dir: path to the job's output which holds duration measures in case of parallelization
        :param output_path: path ot write the json file to
        :return:
        """
        output_dir = os.path.dirname(output_path)
        os.makedirs(output_dir, exist_ok=True)
        result = Program.parse_output(
            input_path=input_path, job_output_dir=job_output_dir
        )
        with open(output_path, "w") as output:
            json.dump(result, output)

    @staticmethod
    def get_error(
        reference_data: t.Dict[str, t.Any], test_data: t.Dict[str, t.Any], **kwargs
    ) -> pd.Series:
        """
        :param reference_data: reference data to compute results by reference to
        :param test_data: test data to compare to the reference data
        :return: the output of pd series with indices as the members for which error it assessed (be it positions in a sequence of sequences) and the values are the error values computed for them
        """
        pass  # is overloaded by implementations n the inherited classes

    @staticmethod
    def get_result(data: t.Dict[str, t.Any], **kwargs) -> pd.Series:
        """
        :param data: dictionary mapping results
        :return: the relevant data to compute error for
        """
        pass  # will be overloaded by inheriting classes

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
        pass  # is overloaded by inheriting classes
