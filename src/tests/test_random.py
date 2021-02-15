import unittest

from Bio import SeqIO
from ete3 import Tree
from samplers import Random


class TestRandom(unittest.TestCase):
    sequences_data_path = "/data/test/seq_data.fas"
    tree_path = "/data/test/tree.nwk"
    aux_dir = "/data/test/aux_random/"

    def test_sample_one(self):
        random_sampler = Random(sequence_data_path=self.sequences_data_path,
                                tree_path=self.tree_path)
        accepted_sample = random_sampler.get_sample(1, self.aux_dir)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(len(accepted_sample), 1)

    def test_sample_no_repeats(self):
        random_sampler = Random(sequence_data_path=self.sequences_data_path,
                                tree_path=self.tree_path)
        accepted_sample = random_sampler.get_sample(3, self.aux_dir)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(len(accepted_sample), 3)

    def test_sample_all(self):
        random_sampler = Random(sequence_data_path=self.sequences_data_path,
                                tree_path=self.tree_path)
        accepted_sample = random_sampler.get_sample(8, self.aux_dir)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(len(accepted_sample), 8)


if __name__ == "__main__":
    unittest.main()
