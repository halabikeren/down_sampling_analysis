import unittest
from Bio import SeqIO
from samplers import CdHit
from ete3 import Tree


class TestCDHIT(unittest.TestCase):

    sequences_data_path = "/data/test/seq_data.fas"
    aux_dir = "/data/test/aux_cdhit/"
    tree = Tree("((A:2,B:1):1,(C:1,(D:8,E:2):2):2);")

    def test_sample_one(self):
        cdhit = CdHit(
            sequences_path=self.sequences_data_path,
            tree=self.tree,
        )
        expected_sample = {"A"}
        accepted_sample = cdhit.get_sample(1, self.aux_dir)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(accepted_sample, expected_sample)

    def test_sample_all(self):
        cdhit = CdHit(
            sequences_path=self.sequences_data_path,
            tree=self.tree,
        )
        accepted_sample = cdhit.get_sample(8, self.aux_dir)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(len(accepted_sample), 8)

    def test_sample_attainable(self):
        cdhit = CdHit(
            sequences_path=self.sequences_data_path,
            tree=self.tree,
        )
        expected_sample = {"A", "C"}
        accepted_sample = cdhit.get_sample(2, self.aux_dir)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(accepted_sample, expected_sample)

    def test_sample_unattainable(self):
        cdhit = CdHit(
            sequences_path=self.sequences_data_path,
            tree=self.tree,
        )
        accepted_sample = cdhit.get_sample(6, self.aux_dir)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(len(accepted_sample), 4)


if __name__ == "__main__":
    unittest.main()