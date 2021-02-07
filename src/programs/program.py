import os
import typing as t
from dataclasses import dataclass
from time import sleep

from pydantic import BaseModel
import json
from enum import Enum

import subprocess

from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv())

import logging

logger = logging.getLogger(__name__)


class Queue(Enum):
    ITAYM = "itaym"
    ITAYMAA = "itaymaa"
    ITAYM1 = "itaym1"
    ITAYM2 = "itaym2"
    ITAYM3 = "itaym3"
    ITAYM4 = "itaym4"


class Job(BaseModel):
    name: str
    sh_dir: str
    output_dir: str
    commands: t.List[str]
    cpus_number: int = 1
    mem_alloc: int = 2  # memory allocation in gb
    queue: Queue = Queue.ITAYM
    priority: int = 0

    def create_sh(self):
        """
        :return:
        """
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.sh_dir, exist_ok=True)

        with open(f"{self.sh_dir}{self.name}.sh", "w") as handle:
            handle.write(
                f"# !/bin/bash\n\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -r y\n#PBS -q {self.queue.value}\n"
            )
            handle.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n")
            handle.write(
                f"#PBS -N {self.name}\n#PBS -e {self.output_dir}\n#PBS -o {self.output_dir}\n"
            )
            handle.write(
                f"#PBS -l select=ncpus={self.cpus_number}:mem={self.mem_alloc}gb\n"
            )
            handle.write("\n".join(self.commands))
            handle.write(f"touch {self.output_dir}{self.name}.touch")

    def submit(
        self, wait_until_complete: bool = False, get_completion_validator: bool = True
    ) -> t.Optional[str]:
        """
        :param wait_until_complete:
        :param get_completion_validator:
        :return:
        """
        self.create_sh()
        logger.info(f"submitting job {self.name}")
        res = os.system(f"qsub -p {self.priority} {self.sh_dir}{self.name}.sh")
        if wait_until_complete:
            while not os.path.exists(f"{self.output_dir}{self.name}.touch"):
                sleep(5)
        if get_completion_validator:
            return f"{self.output_dir}{self.name}.touch"


@dataclass
class Program:
    name: str
    program_exe: str
    cluster_program_exe: str
    input_param_name: str
    output_param_name: str  # maps parameter name to parameter value
    module_to_load: t.Optional[str] = None

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
    ) -> t.Optional[str]:
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
        :return:
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

        additional_params_str = ""
        if additional_params:
            for field in additional_params:
                if "/" in additional_params[field]:
                    additional_params[field] = (
                        f"{os.environ['container_data_dir']}{additional_params[field]}"
                        if not parallelize
                        else f"{cluster_data_dir}{additional_params[field]}"
                    )
            additional_params_str = " ".join(
                [
                    f"{param_name} {additional_params[param_name]}"
                    for param_name in additional_params
                ]
            )
        command = f"{self.program_exe if not parallelize else self.cluster_program_exe} {input_str} {output_str} {additional_params_str} "

        os.makedirs(aux_dir, exist_ok=True)

        if not parallelize:
            process = subprocess.call(
                command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE
            )
            if process != 0:
                raise ValueError(
                    f"command {command} failed to execute due to error {subprocess.PIPE}"
                )
        else:
            commands = [
                f"cd {aux_dir.replace(os.environ['container_data_dir'], cluster_data_dir)}"
            ]
            if self.module_to_load:
                commands.append(f"module load {self.module_to_load}")
            commands.append(command)

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
    def parse_output(output_path: str) -> t.Dict[str, t.Any]:
        if not os.path.exists(output_path):
            raise ValueError(f"output path {output_path} does not exist")
        with open(output_path, "r") as out:
            output_content = out.read()

        result = {"row_output": output_content}
        return result

    @staticmethod
    def write_result_json(input_path: str, output_path: str):
        """
        :param input_path: path to file with the raw program results
        :param output_path: path ot write the json file to
        :return:
        """
        output_dir = os.path.dirname(output_path)
        os.makedirs(output_dir, exist_ok=True)
        result = Program.parse_output(input_path)
        with open(output_path, "w") as output:
            json.dump(result, output)
