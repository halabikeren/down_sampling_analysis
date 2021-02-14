import json
import logging
import os
import sys

import click
from pipeline_utils import PipelineInput, Pipeline

from dotenv import load_dotenv

load_dotenv()


@click.command()
@click.option(
    "--input_path",
    help="path to a json file with input parameters for down sampling pipeline",
    type=click.Path(exists=True, file_okay=True, readable=True),
    required=True,
)
def exec_pipeline(input_path: click.Path):
    """Program to activate down sampling analysis pipeline given parameters input file in a json format.
        For example of the json format parameters, see data/input_json_parameters.txt"""

    # process input json file
    json_dir = os.path.dirname(input_path)
    os.chdir(json_dir)
    with open(input_path, "r") as input_file:
        pipeline_json_input = json.load(input_file)
    os.makedirs(
        os.path.join(os.path.dirname(input_path), pipeline_json_input["pipeline_dir"]),
        exist_ok=True,
    )

    # intialize the logger
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s module: %(module)s function: %(funcName)s line: %(lineno)d %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(
                f"{os.path.join(os.path.dirname(input_path), pipeline_json_input['pipeline_dir'])}pipeline.log"
            ),
        ],
    )
    logger = logging.getLogger(__name__)
    logger.info("Json input has been successfully processed")

    pipeline_input = PipelineInput(**pipeline_json_input)
    logger.info("Json input has been successfully parsed as pipeline input")

    pipeline = Pipeline(pipeline_input)
    logger.info("Pipeline has been successfully initialized")

    pipeline.generate_samples(pipeline_input)
    logger.info("Generated samples successfully.")

    program_names = [prog.value for prog in pipeline_input.programs]
    logger.info(f"Executing programs {program_names}")
    pipeline.execute_programs(pipeline_input)
    logger.info(f"Executed program {program_names} on the generated samples successfully.")

    pipeline_output_path = f"{pipeline_input.pipeline_dir}/pipeline_output.json"
    logger.info(f"Writing pipeline output to {pipeline_output_path}")
    pipeline.write_results(pipeline_output_path)
    logger.info(f"Written pipeline output successfully.")

    logger.info(f"Plotting results")
    pipeline.analyze_results(pipeline_input)
    logger.info(f"Plotted pipeline output successfully.\nPipeline is complete.")



if __name__ == "__main__":
    exec_pipeline()
