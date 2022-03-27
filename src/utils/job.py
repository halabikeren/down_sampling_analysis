import socket
import typing as t

import os
from time import sleep

from pydantic import BaseModel
from .types import Queue

import logging

logger = logging.getLogger(__name__)


class Job(BaseModel):
    name: str
    sh_dir: str
    output_dir: str
    commands: t.List[str]
    cpus_number: int = 1
    mem_alloc: int = 8  # memory allocation in gb
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
            if 'tau' in socket.gethostname() or 'power' in socket.gethostname():
                handle.write(f"{os.environ['conda_act_cmd']}\n")
            handle.write("\n".join(self.commands))
            handle.write(f"\ntouch {self.output_dir}{self.name}.touch")

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
