import os
import shutil
import typing as t

import subprocess
from Bio import SeqIO
from dataclasses import dataclass

from ete3 import Tree


@dataclass
class Sampler:
    sequences_path: str
    tree: Tree
    sequences: t.List[SeqIO.SeqRecord] = None

    def get_sample(
        self, k: int, aux_dir: str, **kwargs
    ) -> t.Union[str, t.List[SeqIO.SeqRecord]]:
        """
        :param k: number of sequences to sample
        :param aux_dir directory to generate auxiliary files in
        :return: either a path to the generated sample or a list of samples sequence names
        """
        self.sequences = list(SeqIO.parse(self.sequences_path, "fasta"))
        if k < 0 or k > len(self.sequences):
            raise ValueError(
                f"sample size {k} is invalid for data of size {len(self.sequences)}"
            )
        return self.sequences_path

    def write_sample(self, k: int, aux_dir: str, output_path: str, **kwargs):
        """
        writes the sampled data to an output file
        :param k required sample size
        :param aux_dir directory to generate auxiliary files in
        :param output_path path to write the sample output to
        :return: None
        """
        sample = self.get_sample(k, aux_dir, **kwargs)
        if type(sample) is str:
            shutil.copyfile(sample, output_path)
        else:
            SeqIO.write(sample, output_path, "fasta")
