import unittest
from pydantic import ValidationError
from samplers import PDA
from ete3 import Tree


class TestPDA(unittest.TestCase):

    def test_class_validation(self):
        tree = Tree("((A:2,B:1):1,(C:1,(D:8,E:2):2):2);")
        with self.assertRaises(ValidationError):
            pda = PDA(tree=tree, taxon_to_weight=[0, 1, 2, 3])

    def test_unweighted_sample(self):
        pda = PDA(tree=Tree("((A:2,B:1):1,(C:1,(D:8,E:2):2):2);"))
        pda.compute_sample(3)
        expected_sample = {"A", "D", "E"}
        accepted_sample = set(pda.sample_subtree.get_leaf_names())
        self.assertEqual(accepted_sample, expected_sample)

    def test_null_weighted_sample(self):
        pda = PDA(tree=Tree("((A:2,B:1):1,(C:1,(D:8,E:2):2):2);"), taxon_to_weight={"A": 1, "B": 1, "C": 1, "D": 1, "E": 1})
        pda.compute_sample(3)
        expected_sample = {"A", "D", "E"}
        accepted_sample = set(pda.sample_subtree.get_leaf_names())
        self.assertEqual(accepted_sample, expected_sample)

    def test_weighted_sample(self):
        pda = PDA(tree=Tree("((A:2,B:1):1,(C:1,(D:8,E:2):2):2);"), taxon_to_weight={"A": 0.5, "B": 0.5, "C": 1, "D": 0.01, "E": 0.5})
        pda.compute_sample(3)
        expected_sample = {"A", "D", "E"}
        accepted_sample = set(pda.sample_subtree.get_leaf_names())
        self.assertEqual(accepted_sample, expected_sample)
        pda.norm_factor = 0.01
        expected_sample = {"A", "C", "E"}
        pda.compute_sample(3, is_weighted=True)
        accepted_sample = set(pda.sample_subtree.get_leaf_names())
        self.assertEqual(accepted_sample, expected_sample)


if __name__ == '__main__':
    unittest.main()
