import os
import typing as t
from dataclasses import dataclass, field
from time import sleep

from pydantic import BaseModel
import json
from enum import Enum

import subprocess


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
        subprocess.run(f"qsub -p {self.priority} {self.sh_dir}{self.name}.sh")
        if wait_until_complete:
            while not os.path.exists(f"{self.output_dir}{self.name}.touch"):
                sleep(5)
        if get_completion_validator:
            return f"{self.output_dir}{self.name}.touch"


@dataclass
class Program:
    name: str
    program_exe: str
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
        priority: int = 0,
        queue: str = "itaym",
        wait_until_complete: bool = False,
        get_completion_validator: bool = True,
    ) -> t.Optional[str]:
        """
        :param input_path:
        :param output_path:
        :param additional_params:
        :param parallelize:
        :param aux_dir:
        :param priority:
        :param queue:
        :param wait_until_complete:
        :param get_completion_validator:
        :return:
        """
        workdir = os.getcwd()
        input_str = f"{self.input_param_name} {input_path}"
        output_str = f"{self.output_param_name} {output_path}"

        additional_params_str = ""
        if additional_params:
            additional_params_str = " ".join(
                [
                    f"{param_name} {additional_params[param_name]}"
                    for param_name in additional_params
                ]
            )
        command = f"{self.program_exe} {input_str} {output_str} {additional_params_str}"

        os.makedirs(aux_dir, exist_ok=True)
        os.chdir(aux_dir)

        if not parallelize:
            process = subprocess.run(command, shell=True, capture_output=True)
            if process.returncode != 0:
                raise ValueError(
                    f"command {command} failed to execute due to error {process.stderr}"
                )
        else:
            commands = []
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
        os.chdir(workdir)

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
