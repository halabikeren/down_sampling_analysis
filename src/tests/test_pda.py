import unittest

from Bio import SeqIO

from samplers import Pda
from ete3 import Tree


class TestPDA(unittest.TestCase):
    sequences_data_path = "/data/test/seq_data.fas"
    tree_path = "/data/test/tree.nwk"
    aux_dir = "/data/test/aux_pda/"
    tree = Tree("((A:2,B:1):1,(C:1,(D:8,E:2):2):2);")

    def test_unweighted_sample(self):
        pda = Pda(sequence_data_path=self.sequences_data_path,
                  tree_path=self.tree_path)
        expected_sample = {"A", "D", "E"}
        accepted_sample = pda.get_sample(3, self.aux_dir)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(accepted_sample, expected_sample)

    def test_unweighted_sample_external(self):
        pda = Pda(sequence_data_path=self.sequences_data_path,
                  tree_path=self.tree_path)
        expected_sample = {"A", "D", "E"}
        accepted_sample = pda.get_sample(3, self.aux_dir, use_external=True)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(accepted_sample, expected_sample)

    def test_null_weighted_sample(self):
        pda = Pda(sequence_data_path=self.sequences_data_path,
                  tree_path=self.tree_path,
                  taxon_to_weight={"A": 1, "B": 1, "C": 1, "D": 1, "E": 1})
        expected_sample = {"A", "D", "E"}
        accepted_sample = pda.get_sample(3, self.aux_dir)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(accepted_sample, expected_sample)

    def test_null_weighted_sample_external(self):
        pda = Pda(sequence_data_path=self.sequences_data_path,
                  tree_path=self.tree_path,
                  taxon_to_weight={"A": 1, "B": 1, "C": 1, "D": 1, "E": 1})
        expected_sample = {"A", "D", "E"}
        accepted_sample = pda.get_sample(3, self.aux_dir)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(accepted_sample, expected_sample)

    def test_weighted_sample_no_normalization(self):
        pda = Pda(sequence_data_path=self.sequences_data_path,
                  tree_path=self.tree_path,
                  taxon_to_weight={"A": 0.5, "B": 0.5, "C": 1, "D": 0.01, "E": 0.5})
        expected_sample = {"A", "D", "E"}
        accepted_sample = pda.get_sample(3, self.aux_dir, is_weighted=True)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(accepted_sample, expected_sample)

    def test_weighted_sample_no_normalization_external(self):
        pda = Pda(sequence_data_path=self.sequences_data_path,
                  tree_path=self.tree_path,
                  taxon_to_weight={"A": 0.5, "B": 0.5, "C": 1, "D": 0.01, "E": 0.5})
        expected_sample = {"A", "D", "E"}
        accepted_sample = pda.get_sample(3, self.aux_dir, use_external=True)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(accepted_sample, expected_sample)

    def test_weighted_sample_with_normalization(self):
        pda = Pda(sequence_data_path=self.sequences_data_path,
                  tree_path=self.tree_path,
                  taxon_to_weight={"A": 0.5, "B": 0.5, "C": 1, "D": 0.01, "E": 0.5})
        pda.norm_factor = 0.01
        expected_sample = {"A", "C", "E"}
        accepted_sample = pda.get_sample(3, self.aux_dir, is_weighted=True)
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(accepted_sample, expected_sample)

    def test_weighted_sample_with_normalization_external(self):
        pda = Pda(
            sequence_data_path=self.sequences_data_path,
            tree_path=self.tree_path,
            taxon_to_weight={"A": 0.5, "B": 0.5, "C": 1, "D": 0.01, "E": 0.5},
        )
        pda.norm_factor = 0.01
        expected_sample = {"A", "C", "E"}
        accepted_sample = pda.get_sample(
            3, self.aux_dir, is_weighted=True, use_external=True
        )
        if type(accepted_sample) is str:
            accepted_sample = list(SeqIO.parse(accepted_sample, "fasta"))
        accepted_sample = set([record.id for record in accepted_sample])
        self.assertEqual(accepted_sample, expected_sample)

    def test_computed_weights(self):
        pda = Pda(sequence_data_path=self.sequences_data_path,
                  tree_path=self.tree_path)
        pda.compute_taxon_weights(f"/data/test/aligned_seq_data.fas")
        self.assertEqual(
            pda.taxon_to_weight,
            {
                "A": 0.68359375,
                "B": 0.68359375,
                "C": 0.3645833333333333,
                "D": 0.5911458333333334,
                "E": 0.7356770833333334,
                "F": 0.7356770833333334,
                "G": 0.7356770833333334,
                "H": 0.7356770833333334,
            },
        )


if __name__ == "__main__":
    unittest.main()
