import shutil
import typing as t
import os
import unittest
import json
from Bio import SeqIO
from pipeline_utils import PipelineInput, Pipeline


class TestPipeline(unittest.TestCase):
    json_path = "/data/test/input.json"

    def generate_pipeline_input(self) -> PipelineInput:
        os.chdir(os.path.dirname(self.json_path))
        with open(self.json_path, "r") as json_file:
            pipeline_json_input = json.load(json_file)
        os.makedirs(pipeline_json_input["pipeline_dir"], exist_ok=True)
        pipeline_input = PipelineInput(**pipeline_json_input)
        return pipeline_input

    def test_pipeline_input_parsing(self):
        pipeline_input = self.generate_pipeline_input()

    def test_pipeline_parsing(self):
        pipeline_input = self.generate_pipeline_input()
        pipeline = Pipeline(pipeline_input)

    def call_samples_generation(self) -> t.Union[PipelineInput, Pipeline]:
        pipeline_input = self.generate_pipeline_input()
        pipeline = Pipeline(pipeline_input)
        pipeline.generate_samples(pipeline_input)
        return pipeline_input, pipeline

    def test_samples_generation(self):
        pipeline_input, pipeline = self.call_samples_generation()
        full_data_size = len(
            list(SeqIO.parse(pipeline_input.unaligned_sequence_data_path, "fasta"))
        )
        for fraction in pipeline_input.sampling_fractions:
            for method in pipeline_input.sampling_methods:
                if method.value[0] != "cdhit" and fraction != 0.75:
                    path = pipeline.samples_info[fraction][method.value][
                        "unaligned_sequence_data_path"
                    ]
                    self.assertTrue(os.path.exists(path))
                    records = list(SeqIO.parse(path, "fasta"))
                    self.assertEqual(len(records), int(full_data_size * fraction))

    def test_program_exec(self):
        pipeline_input, pipeline = self.call_samples_generation()
        pipeline.execute_programs(pipeline_input)
        for program_name in pipeline_input.programs:
            for fraction in pipeline_input.sampling_fractions:
                for method in pipeline_input.sampling_methods:
                    program_info = pipeline.samples_info[fraction][method.value][
                        "programs_performance"
                    ][program_name.value]
                    self.assertIsInstance(program_info["result"], dict)

    def tearDown(self):
        with open(TestPipeline.json_path, "r") as json_file:
            pipeline_json_input = json.load(json_file)
        if os.path.exists(pipeline_json_input["pipeline_dir"]):
            os.system(f"rm -rf {pipeline_json_input['pipeline_dir']}")


if __name__ == "__main__":
    unittest.main()
