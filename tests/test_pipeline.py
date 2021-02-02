import os
import unittest
import json
import subprocess
from Bio import SeqIO
from pipeline_utils import PipelineInput, Pipeline


class TestPipeline(unittest.TestCase):
    json_path = f"{os.path.dirname(os.path.realpath(__file__))}/data/input.json"

    def test_pipeline_input_parsing(self):
        with open(self.json_path, "r") as json_file:
            pipeline_json_input = json.load(json_file)
        if not os.path.exists(pipeline_json_input["pipeline_dir"]):
            subprocess.run(
                f"mkdir -p {pipeline_json_input['pipeline_dir']}",
                shell=True,
                capture_output=True,
            )
        pipeline_input = PipelineInput(**pipeline_json_input)

    def test_pipeline_parsing(self):
        with open(self.json_path, "r") as json_file:
            pipeline_json_input = json.load(json_file)
        if not os.path.exists(pipeline_json_input["pipeline_dir"]):
            subprocess.run(
                f"mkdir -p {pipeline_json_input['pipeline_dir']}",
                shell=True,
                capture_output=True,
            )
        pipeline_input = PipelineInput(**pipeline_json_input)
        pipeline = Pipeline(pipeline_input)

    def test_samples_generation(self):
        with open(self.json_path, "r") as json_file:
            pipeline_json_input = json.load(json_file)
        if not os.path.exists(pipeline_json_input["pipeline_dir"]):
            subprocess.run(
                f"mkdir -p {pipeline_json_input['pipeline_dir']}",
                shell=True,
                capture_output=True,
            )
        pipeline_input = PipelineInput(**pipeline_json_input)
        pipeline = Pipeline(pipeline_input)
        pipeline.generate_samples(
            pipeline_input.sampling_fractions, pipeline_input.sampling_methods
        )
        full_data_size = len(
            list(SeqIO.parse(pipeline_input.sequence_data_path, "fasta"))
        )
        for fraction in pipeline_input.sampling_fractions:
            for method in pipeline_input.sampling_methods:
                if method.value[0] != "cdhit" and fraction != 0.75:
                    path = pipeline.samples_info[fraction][method.value[0]]["path"]
                    self.assertTrue(os.path.exists(path))
                    records = list(SeqIO.parse(path, "fasta"))
                    self.assertEqual(len(records), int(full_data_size * fraction))

    def tearDown(self):
        with open(TestPipeline.json_path, "r") as json_file:
            pipeline_json_input = json.load(json_file)
        subprocess.run(
            f"rm -rf {pipeline_json_input['pipeline_dir']}",
            shell=True,
            capture_output=True,
        )


if __name__ == "__main__":
    unittest.main()
