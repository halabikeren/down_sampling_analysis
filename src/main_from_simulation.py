import json
import logging
import os
import sys

import click
from utils import BaseTools, Job

from dotenv import load_dotenv

load_dotenv()


@click.command()
@click.option(
    "--input_path",
    help="path to a json file with input parameters for simulations for pipeline execution",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
)
def exec_pipeline_on_simulations(input_path: click.Path):
    """Program to simulate multiple datasets and then submit pipeline jobs for each one
       For example of the json format parameters, see data/test/simulation.json"""

    # process input json file
    with open(input_path, "r") as input_file:
        simulation_params = json.load(input_file)

    # intialize the logger
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s module: %(module)s function: %(funcName)s line: %(lineno)d %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(
                f"{os.path.dirname(input_path)}/simulations.log"
            ),
        ],
    )
    logger = logging.getLogger(__name__)

    logger.info(f"Simulating data in {simulation_params['simulations_output_dir']}")
    pipeline_input_json_paths = BaseTools.simulate(simulation_params=simulation_params)
    logger.info(f"Simulation is complete.")

    logger.info(f"submitting pipeline jobs for the simulated data")
    for json_path in pipeline_input_json_paths:
        aux_dir = f"{os.path.dirname(json_path)}/job_aux/"
        job = Job(
            name="pipeline_on_simulated_data",
            sh_dir=aux_dir,
            output_dir=aux_dir,
            commands=[f"python /groups/itay_mayrose/halabikeren/down_sampling_analysis/src/main.py --input_path={json_path}"],
            priority=simulation_params["priority"],
            queue=simulation_params["queue"],
        )
        job.submit(
            wait_until_complete=False,
            get_completion_validator=False,
        )
    logger.info(f"Job submission is complete")


if __name__ == "__main__":
    exec_pipeline_on_simulations()
