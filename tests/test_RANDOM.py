import os
import unittest
from pydantic import ValidationError
from samplers import RANDOM


class TestRANDOM(unittest.TestCase):

    def test_class_validation(self):
        with self.assertRaises(ValidationError):
            random_sampler = RANDOM(all_sequences_path=None)

    def test_sample_one(self):
        sequence_data_path = f"{os.path.dirname(os.path.realpath(__file__))}/data/seq_data.fas"
        random_sampler = RANDOM(all_sequences_path=sequence_data_path)
        random_sampler.compute_sample(1)
        accepted_sample = random_sampler.sample_members
        self.assertEqual(len(accepted_sample), 1)

    def test_sample_no_repeats(self):
        sequence_data_path = f"{os.path.dirname(os.path.realpath(__file__))}/data/seq_data.fas"
        random_sampler = RANDOM(all_sequences_path=sequence_data_path)
        random_sampler.compute_sample(3)
        accepted_sample = set(random_sampler.sample_members)
        self.assertEqual(len(accepted_sample), 3)


    def test_sample_all(self):
        sequence_data_path = f"{os.path.dirname(os.path.realpath(__file__))}/data/seq_data.fas"
        random_sampler = RANDOM(all_sequences_path=sequence_data_path)
        random_sampler.compute_sample(5)
        expected_sample = {"A", "B", "C", "D", "E"}
        accepted_sample = set(random_sampler.sample_members)
        self.assertEqual(accepted_sample, expected_sample)


if __name__ == '__main__':
    unittest.main()
