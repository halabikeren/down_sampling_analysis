from Bio import SeqIO
import typing as t
from dataclasses import dataclass
import random
from .sampler import Sampler


@dataclass
class Random(Sampler):
    def get_sample(
        self, k: int, aux_dir: str, **kwargs
    ) -> t.Union[str, t.List[SeqIO.SeqRecord]]:
        """
        :param k: number of sequences to sample
        :param aux_dir directory to generate auxiliary files in
        :return: either a path to the generated sample or a list of samples sequence names
        """
        sample = super(Random, self).get_sample(k, aux_dir)
        if k < len(self.sequences):
            sample = random.sample(self.sequences, k)
        return sample