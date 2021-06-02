import typing as t
import os
import unittest
import json
from pipeline_utils import PipelineInput, Pipeline
from utils import SimulationInput, SimulationTools


class TestSimulationPipeline(unittest.TestCase):
    simulation_params_path = "/data/test/simulations.json"
    pipeline_input_json_paths = []

    def generate_simulations(self):
        with open(self.simulation_params_path, "r") as input_file:
            simulation_params = json.load(input_file)
        os.makedirs(simulation_params["simulations_output_dir"], exist_ok=True)
        simulation_input = SimulationInput(**simulation_params)
        self.pipeline_input_json_paths = SimulationTools.simulate(
            simulation_input=simulation_input
        )

    def generate_pipeline_inputs(self) -> t.List[PipelineInput]:
        if len(self.pipeline_input_json_paths) == 0:
            self.generate_simulations()
        pipeline_inputs = []
        for json_path in self.pipeline_input_json_paths:
            self.assertTrue(os.path.exists(json_path))
            os.chdir(os.path.dirname(json_path))
            with open(json_path, "r") as json_file:
                pipeline_json_input = json.load(json_file)
            os.makedirs(pipeline_json_input["pipeline_dir"], exist_ok=True)
            pipeline_input = PipelineInput(**pipeline_json_input)
            pipeline_inputs.append(pipeline_input)
        return pipeline_inputs

    def test_pipeline_parsing(self):
        pipeline_inputs = self.generate_pipeline_inputs()
        for pipeline_input in pipeline_inputs:
            pipeline = Pipeline(pipeline_input)

    def tearDown(self):
        for json_path in self.pipeline_input_json_paths:
            with open(json_path, "r") as json_file:
                pipeline_json_input = json.load(json_file)
            if os.path.exists(pipeline_json_input["pipeline_dir"]):
                os.system(f"rm -rf {pipeline_json_input['pipeline_dir']}")
        self.pipeline_input_json_paths = []


if __name__ == "__main__":
    unittest.main()
