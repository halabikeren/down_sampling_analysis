import os
import unittest
from pydantic import ValidationError
from samplers import CDHIT


class TestCDHIT(unittest.TestCase):

    def test_class_validation(self):
        with self.assertRaises(ValidationError):
            cdhit = CDHIT(sequences_path=None, aux_dir=f"{os.path.dirname(os.path.realpath(__file__))}/aux/cdhit/")

    def test_sample_one(self):
        sequence_data_path = f"{os.path.dirname(os.path.realpath(__file__))}/data/seq_data.fas"
        cdhit = CDHIT(sequences_path=sequence_data_path, aux_dir=f"{os.path.dirname(os.path.realpath(__file__))}/aux/cdhit/")
        cdhit.compute_sample(1)
        expected_sample = ["A"]
        accepted_sample = cdhit.sample_members
        self.assertEqual(accepted_sample, expected_sample)

    def test_sample_all(self):
        sequence_data_path = f"{os.path.dirname(os.path.realpath(__file__))}/data/seq_data.fas"
        cdhit = CDHIT(sequences_path=sequence_data_path, aux_dir=f"{os.path.dirname(os.path.realpath(__file__))}/aux/cdhit/")
        cdhit.compute_sample(4)
        expected_sample = ["A", "C", "D", "E"]
        accepted_sample = cdhit.sample_members
        self.assertEqual(accepted_sample, expected_sample)

    def test_sample_attainable(self):
        sequence_data_path = f"{os.path.dirname(os.path.realpath(__file__))}/data/seq_data.fas"
        cdhit = CDHIT(sequences_path=sequence_data_path, aux_dir=f"{os.path.dirname(os.path.realpath(__file__))}/aux/cdhit/")
        cdhit.compute_sample(2)
        expected_sample = ["A", "C"]
        accepted_sample = cdhit.sample_members
        self.assertEqual(accepted_sample, expected_sample)

    def test_sample_unattainable(self):
        sequence_data_path = f"{os.path.dirname(os.path.realpath(__file__))}/data/seq_data.fas"
        cdhit = CDHIT(sequences_path=sequence_data_path, aux_dir=f"{os.path.dirname(os.path.realpath(__file__))}/aux/cdhit/")
        cdhit.compute_sample(5)
        expected_sample = ["A", "C", "D", "E"]
        accepted_sample = cdhit.sample_members
        self.assertEqual(accepted_sample, expected_sample)


if __name__ == '__main__':
    res = os.system(f"mkdir -p {os.path.dirname(os.path.realpath(__file__))}/aux/cdhit/")
    unittest.main()
