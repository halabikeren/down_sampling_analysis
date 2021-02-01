import os
import unittest

from Bio import SeqIO
from ete3 import Tree
from samplers import Random


class TestRandom(unittest.TestCase):

    sequences_data_path = f"{os.path.dirname(os.path.realpath(__file__))}/data/seq_data.fas"
    aux_dir = f"{os.path.dirname(os.path.realpath(__file__))}/aux/random/"
    tree = Tree("((A:2,B:1):1,(C:1,(D:8,E:2):2):2);")


    def test_sample_one(self):
        random_sampler = Random(sequences_path=self.sequences_data_path,
                                tree=self.tree,
                                aux_dir=self.aux_dir)
        accepted_sample = random_sampler.get_sample(1)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(len(accepted_sample), 1)

    def test_sample_no_repeats(self):
        random_sampler = Random(sequences_path=self.sequences_data_path,
                                tree=self.tree,
                                aux_dir=self.aux_dir)
        accepted_sample = random_sampler.get_sample(3)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(len(accepted_sample), 3)

    def test_sample_all(self):
        random_sampler = Random(sequences_path=self.sequences_data_path,
                                tree=self.tree,
                                aux_dir=self.aux_dir)
        expected_sample = {"A", "B", "C", "D", "E"}
        accepted_sample = random_sampler.get_sample(5)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(accepted_sample, expected_sample)


if __name__ == "__main__":
    unittest.main()
