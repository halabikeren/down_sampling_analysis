from Bio import SeqIO
import typing as t
from dataclasses import dataclass
import random

from ete3 import Tree

from .sampler import Sampler


@dataclass
class Random(Sampler):

    def __init__(self, sequence_data_path: str, tree_path: str, exclude_a_ref_sequence: bool = False, sequences: t.Optional[t.List[SeqIO.SeqRecord]] = None):
        super(Random, self).__init__(sequence_data_path=sequence_data_path, tree_path=tree_path, sequences=sequences)

    def get_sample(
        self, k: int, aux_dir: str, **kwargs
    ) -> t.List[SeqIO.SeqRecord]:
        """
        :param k: number of sequences to sample
        :param aux_dir directory to generate auxiliary files in
        :return: either a path to the generated sample or a list of samples sequence names
        """
        sample = super(Random, self).get_sample(k, aux_dir)
        if k < len(self.sequences):
            sample_size = k-1 if self.exclude_a_ref_sequence else k
            sample = random.sample(self.sequences, sample_size)
            if self.exclude_a_ref_sequence:
                sample.append(self.saved_sequence)  # saved sequence must be in every sample
        return sample
