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

    json_dir = os.path.dirname(input_path)
    os.chdir(json_dir)
    with open(input_path, "r") as input_file:
        pipeline_json_input = json.load(input_file)
    os.makedirs(
        os.path.join(os.path.dirname(input_path), pipeline_json_input["pipeline_dir"]),
        exist_ok=True,
    )

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(module)s %(funcName)s %(lineno)d %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(
                f"{os.path.join(os.path.dirname(input_path), pipeline_json_input['pipeline_dir'])}pipeline.log"
            ),
        ],
    )
    logger = logging.getLogger(__name__)
    logger.info(f"Json input has been successfully processed")

    pipeline_input = PipelineInput(**pipeline_json_input)
    logger.info(f"Json input has been successfully parsed as pipeline input")

    pipeline = Pipeline(pipeline_input)
    logger.info(f"Pipeline has been successfully initialized")

    pipeline.generate_samples(pipeline_input)
    samples_paths_description = "\n".join(
        [
            pipeline.samples_info[fraction][method.value][
                "unaligned_sequence_data_path"
            ]
            for fraction in pipeline_input.sampling_fractions
            for method in pipeline_input.sampling_methods
        ]
    )
    logger.info(
        f"Generated samples successfully. Sample files are available at {samples_paths_description}"
    )

    pipeline.execute_programs(pipeline_input)
    logger.info(
        f"Executed program {[program_name.value for program_name in pipeline_input.programs]} on the generated "
        f"samples successfully. "
    )


if __name__ == "__main__":
    exec_pipeline()