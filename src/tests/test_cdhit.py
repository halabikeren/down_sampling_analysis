import unittest
from Bio import SeqIO
from samplers import CdHit
from ete3 import Tree


class TestCDHIT(unittest.TestCase):

    sequences_data_path = "/data/test/seq_data.fas"
    tree_path = "/data/test/tree.nwk"
    aux_dir = "/data/test/aux_cdhit/"

    def test_sample_one(self):
        cdhit = CdHit(
            sequence_data_path=self.sequences_data_path,
            tree_path=self.tree_path,
        )
        expected_sample = {"A"}
        accepted_sample = cdhit.get_sample(1, self.aux_dir)
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(accepted_sample, expected_sample)

    def test_sample_all(self):
        cdhit = CdHit(
            sequence_data_path=self.sequences_data_path,
            tree_path=self.tree_path,
        )
        accepted_sample = cdhit.get_sample(8, self.aux_dir)
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(len(accepted_sample), 8)

    def test_sample_attainable(self):
        cdhit = CdHit(
            sequence_data_path=self.sequences_data_path,
            tree_path=self.tree_path,
        )
        expected_sample = {"A", "C"}
        accepted_sample = cdhit.get_sample(2, self.aux_dir)
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(accepted_sample, expected_sample)

    def test_sample_unattainable(self):
        cdhit = CdHit(
            sequence_data_path=self.sequences_data_path,
            tree_path=self.tree_path,
        )
        accepted_sample = cdhit.get_sample(5, self.aux_dir)
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(len(accepted_sample), 4)


if __name__ == "__main__":
    unittest.main()
