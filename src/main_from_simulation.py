import json
import logging
import os
import sys
from time import sleep
import pandas as pd
import click
from utils import SimulationInput, BaseTools, Job
import matplotlib.pyplot as plt
from dotenv import load_dotenv
import seaborn as sns

load_dotenv()


def plot_large_scale_error(df: pd.DataFrame, output_path: str, use_relative_error: bool = False):
    """
    :param df: dataframe with a column "replicate" and other, program specific, columns
    :param output_path: path to plot the data in
    :param use_relative_error: boolean indicating weather the error should be relative or absolute
    :return: none
    """

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.grid(False)
    prefix = "absolute"
    if use_relative_error:
        prefix = "relative"
    ncols = 2 if f"{prefix}_error_to_ref" in df.columns and f"{prefix}_error_to_full" in df.columns else 1
    fig, axis = plt.subplots(
        nrows=1,
        ncols=ncols,
        sharex="none",
        sharey="none",
        figsize=[ncols * 8.5 + 2 + 2, 7.58 + 2],
        frameon=True,
    )
    if f"{prefix}_error_to_ref" in df.columns:
        sns.boxplot(ax=axis[0], y=f"{prefix}_error_to_ref", x="sampling_fraction", data=df,
                    palette="colorblind",
                    hue="sampling_method")
        axis[0].set_ylabel(f"mean {prefix} error ({len(df['replicate'].unique())} replicates)")
        axis[0].set_xlabel("sampling fraction")
        axis[0].set_title("reference: simulated rates")
    if f"{prefix}_error_to_full" in df.columns:
        sns.boxplot(ax=axis[1], y=f"{prefix}_error_to_full", x="sampling_fraction", data=df,
                    palette="colorblind",
                    hue="sampling_method")
        axis[1].set_ylabel(f"mean {prefix} error ({len(df['replicate'].unique())} replicates)")
        axis[1].set_xlabel("sampling fraction")
        axis[1].set_title("reference: inferred rates based on full dataset")
    fig.subplots_adjust()
    fig.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", transparent=True)
    plt.clf()


def plot_large_scale_bias(df: pd.DataFrame, output_path: str):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.grid(False)
    fig, axis = plt.subplots(
        nrows=1,
        ncols=2,
        sharex="none",
        sharey="none",
        figsize=[2 * 8.5 + 2 + 2, 7.58 + 2],
        frameon=True,
    )
    # plot rate next to simulated
    sns.boxplot(ax=axis[0], y="full_bias", x="sampling_fraction", data=df,
                palette="colorblind",
                hue="sampling_method")
    # plot rate next to full
    axis[0].set_ylabel(f"rates diff ({len(df['replicate'].unique())} replicates)")
    axis[0].set_xlabel("sampling fraction")
    axis[0].set_title("reference: simulated rates")

    sns.boxplot(ax=axis[1], y="simulated_bias", x="sampling_fraction", data=df,
                palette="colorblind",
                hue="sampling_method")
    axis[1].set_ylabel(f"rates diff ({len(df['replicate'].unique())} replicates)")
    axis[1].set_xlabel("sampling fraction")
    axis[1].set_title("reference: inferred rates based on full dataset")

    fig.subplots_adjust()
    fig.tight_layout()
    plt.savefig(output_path, bbox_inches="tight", transparent=True)
    plt.clf()


def plot_large_scale_samples_overlap(df: pd.DataFrame, output_path: str):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.grid(False)
    fig, axis = plt.subplots(
        nrows=1,
        ncols=1,
        sharex="none",
        sharey="none",
        figsize=[1 * 8.5 + 2 + 2, 7.58 + 2],
        frameon=True,
    )
    sns.boxplot(ax=axis, y="overlap_fraction", x="sampling_fraction", data=df,
                palette="colorblind",
                hue="compared_methods")
    axis.set_ylabel(f"% overlap ({len(df['replicate'].unique())} replicates)")
    axis.set_xlabel("sampling fraction")
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
    os.makedirs(
        simulation_params["simulations_output_dir"],
        exist_ok=True,
    )

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
    logger.info("Json input has been successfully processed")

    logger.info(f"Processing simulation input from {input_path}")
    simulation_input = SimulationInput(**simulation_params)
    logger.info("Json input has been successfully parsed as simulation input")

    logger.info(f"Simulating data in {simulation_input.simulations_output_dir}")
    simulations_exist = False
    simulations_exec_complete = False
    repetitions_num = simulation_input.nrep
    if os.path.exists(simulation_input.simulations_output_dir) and os.listdir(
            simulation_input.simulations_output_dir) == repetitions_num:
        simulations_exist = True
        all_exist = True
        for path in os.listdir(simulation_params.simulations_output_dir):
            completion_validator = f"{simulation_params.simulations_output_dir}/{path}/job_aux/pipeline_on_simulated_data.touch"
            if not os.path.exists(completion_validator):
                all_exist = False
                break
        if all_exist:
            simulations_exec_complete = True

    if not simulations_exist:
        pipeline_input_json_paths = BaseTools.simulate(simulation_input=simulation_input)
        simulations_dirs = [f"{os.path.dirname(json_path)}/" for json_path in pipeline_input_json_paths]
        logger.info(f"Simulation is complete.")

    else:
        simulations_dirs = [f"{simulation_input.simulations_output_dir}/{path}/" for path in
                            os.listdir(simulation_input.simulations_output_dir)]

    if not simulations_exec_complete:
        logger.info(f"submitting pipeline jobs for the simulated data")
        completion_validators = []
        for simulations_dir in simulations_dirs:
            aux_dir = f"{simulations_dir}/job_aux/"
            json_path = f"{simulations_dir}/input.json"
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
    paths = [path for path in os.listdir(simulation_input.simulations_output_dir) if "rep" in path]
    overlap_dfs = []
    for path in paths:
        overlap_df_path = f"{simulation_input.simulations_output_dir}/{path}/pipeline_dir/samples/samples_overlap.csv"
        overlap_df = pd.read_csv(overlap_df_path)
        overlap_df["replicate"] = path
        overlap_df["compared_methods"] = (overlap_df["method_1"], overlap_df["method_2"])

        overlap_dfs.append(overlap_df)
    full_overlap_df = pd.concat(overlap_dfs)
    plot_large_scale_samples_overlap(df=full_overlap_df, output_path=f"{simulation_input.simulations_output_dir}/samples_overlap.svg")



    for program in simulation_input.programs:
        data = []
        paths = [path for path in os.listdir(simulation_input.simulations_output_dir) if "rep" in path]
        for path in paths:
            df_path = f"{simulation_input.simulations_output_dir}/{path}/pipeline_dir/tables/{program}_summary.csv"
            try:
                rep_data = pd.read_csv(df_path)
                rep_data["replicate"] = path
                data.append(rep_data)
            except Exception as e:
                logger.error(f"Failed to load dataframe from {df_path} due to error {e}")
        full_df = pd.concat(data)
        full_df["full_bias"] = full_df["full_result"] - full_df["result"]
        full_df["simulated_bias"] = full_df["simulated"] - full_df["result"]
        full_df_grouped = full_df.groupby(["replicate", "sampling_fraction", "sampling_method"]).mean().reset_index()
        full_df_grouped.to_csv(f"{simulation_params['simulations_output_dir']}/{program}_aggregated_data.csv")

        # plot large scale data
        plot_large_scale_error(df=full_df_grouped,
                               output_path=f"{simulation_input.simulations_output_dir}/{program}_absolute_error.svg",
                               use_relative_error=False)
        plot_large_scale_error(df=full_df_grouped,
                               output_path=f"{simulation_input.simulations_output_dir}/{program}_relative_error.svg",
                               use_relative_error=True)
        plot_large_scale_bias(df=full_df_grouped,
                              output_path=f"{simulation_input.simulations_output_dir}/{program}_bias.svg")


if __name__ == "__main__":
    exec_pipeline_on_simulations()
