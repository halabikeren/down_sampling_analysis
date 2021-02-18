import json
import logging
import os
import sys
from time import sleep
import pandas as pd
import click
from utils import BaseTools, Job
import matplotlib.pyplot as plt
from dotenv import load_dotenv
import seaborn as sns

load_dotenv()


def plot_large_scale_results(df: pd.DataFrame, output_path: str):
    """
    :param df: dataframe with a column "replicate" and other, program specific, columns
    :param output_path: path to plot the data in
    :return: none
    """

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.grid(False)
    ncols = 2 if not "relative_accuracy_to_ref" in df.columns and "relative_accuracy_to_full" in column else 1
    fig, axis = plt.subplots(
        nrows=1,
        ncols=ncols,
        sharex="none",
        sharey="none",
        figsize=[ncols * 8.5 + 2 + 2, 7.58 + 2],
        frameon=True,
    )
    if "relative_accuracy_to_ref" in df.columns:
        sns.boxplot(ax=axis[0], y="relative_accuracy_to_ref", x="sampling_fraction", data=df.groupby(["replicate"]).mean().reset_index(),
                           palette="colorblind",
                           hue="sampling_method")
        axis[0].set_ylabel(f"mean relative accuracy ({len(df['replicate'].unique())} replicates)")
        axis[0].set_xlabel("sampling fraction")
        axis[0].set_ylim(bottom=0, top=1)
        axis[0].set_title("reference: simulated rates")
    if "relative_accuracy_to_full" in df.columns:
        sns.boxplot(ax=axis[1], y="relative_accuracy_to_full", x="sampling_fraction", data=df.groupby(["replicate"]).mean().reset_index(),
                           palette="colorblind",
                           hue="sampling_method")
        axis[0].set_ylabel(f"mean relative accuracy ({len(df['replicate'].unique())} replicates)")
        axis[0].set_xlabel("sampling fraction")
        axis[0].set_title("reference: inferred rates based on full dataset")
    fig.subplots_adjust()
    fig.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", transparent=True)
    plt.clf()


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
    simulations_exist = False
    simulations_exec_complete = False
    repetitions_num = int(simulation_params["nrep"])
    if os.path.exists(simulation_params['simulations_output_dir']) and os.listdir(
            simulation_params['simulations_output_dir']) == repetitions_num:
        simulations_exist = True
        all_exist = True
        for path in os.listdir(simulation_params['simulations_output_dir']):
            completion_validator = f"{simulation_params['simulations_output_dir']}/{path}/job_aux/pipeline_on_simulated_data.touch"
            if not os.path.exists(completion_validator):
                all_exist = False
                break
        if all_exist:
            simulations_exec_complete = True

    if not simulations_exist:
        pipeline_input_json_paths = BaseTools.simulate(simulation_params=simulation_params)
        aux_dirs = [f"{os.path.dirname(json_path)}/job_aux/" for json_path in pipeline_input_json_paths]
        logger.info(f"Simulation is complete.")

    else:
        aux_dirs = [f"{simulation_params['simulations_output_dir']}/{path}/job_aux/" for path in
                    os.listdir(simulation_params['simulations_output_dir'])]

    if not simulations_exec_complete:
        logger.info(f"submitting pipeline jobs for the simulated data")
        completion_validators = []
        for aux_dir in aux_dirs:
            if not os.path.exists(f"{aux_dir}/pipeline_on_simulated_data.touch"):
                job = Job(
                    name="pipeline_on_simulated_data",
                    sh_dir=aux_dir,
                    output_dir=aux_dir,
                    commands=[
                        f"python /groups/itay_mayrose/halabikeren/down_sampling_analysis/src/main.py --input_path={json_path}"],
                    priority=simulation_params["priority"],
                    queue=simulation_params["queue"],
                )
                completion_validators.append(job.submit(
                    wait_until_complete=False,
                    get_completion_validator=True,
                ))
        logger.info(f"Job submission is complete")

        # wait for jobs to complete
        for validator in completion_validators:
            while not os.path.exists(validator):
                sleep(60)

    # analyze large scale results
    for program in simulation_params["programs"]:
        data = []
        output_path = f"{simulation_params['simulations_output_dir']}/{program}.svg"
        for path in os.listdir(simulation_params['simulations_output_dir']):
            df_path = f"{simulation_params['simulations_output_dir']}/{path}/pipeline_dir/tables/{program}_summary.csv"
            try:
                rep_data = pd.read_csv(df_path)
                rep_data["replicate"] = path
                data.append(rep_data)
            except Exception as e:
                logger.error(f"Failed to load dataframe from {df_path} due to error {e}")
        full_df = pd.concat(data)
        # plot large scale data
        plot_large_scale_results(df=full_df, output_path=output_path)


if __name__ == "__main__":
    exec_pipeline_on_simulations()
