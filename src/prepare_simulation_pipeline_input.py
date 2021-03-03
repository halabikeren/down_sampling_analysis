import logging
import os
import pickle
import re
import json
import sys
from time import sleep
import click
from Bio import SeqIO
from random import sample
from utils import SequenceDataType, AlignmentMethod, BaseTools
from programs import program_to_callable
import typing as t


def remove_duplicates(sequence_records):
    i = 0
    gb_regex = re.compile("gb\:(.*?)\|")
    gbs = set()
    while i < len(sequence_records):
        record = sequence_records[i]
        gb = gb_regex.search(record.description).group(1)
        if gb not in gbs:
            gbs.add(gb)
            i += 1
        else:
            sequence_records.pop(i)


def convert_names(translator_path: str, records: t.List[SeqIO.SeqRecord]):
    new_to_old = dict()
    i = 1
    for record in records:
        old = record.description
        new = f"S{i}"
        i += 1
        record.id = record.name = record.description = new
        new_to_old[new] = old
    with open(translator_path, "wb") as outfile:
        pickle.dump(new_to_old, outfile)


def sample_data(full_data: t.List[SeqIO.SeqRecord], output_dir: str, required_data_size: int, num_of_repeats: int) -> t.List[str]:
    if required_data_size > len(full_data):
        raise ValueError(f"full data size {len(full_data)} cannot be down-sampled to {required_data_size}")
    os.makedirs(output_dir, exist_ok=True)
    samples_paths = []
    for i in range(1, num_of_repeats+1):
        os.makedirs(f"{output_dir}/sample_{i}/", exist_ok=True)
        sampled_data = sample(full_data, required_data_size)
        output_path = f"{output_dir}/sample_{i}/sequences_old_names.fasta"
        SeqIO.write(sampled_data, output_path, "fasta")
        convert_names(translator_path=f"{output_dir}/sample_{i}/new_to_old_names_map.pickle", records=sampled_data)
        output_path = f"{output_dir}/sample_{i}/sequences_new_names.fasta"
        SeqIO.write(sampled_data, output_path, "fasta")
        samples_paths.append(output_path)
    return samples_paths

def run_program(sequence_data_path: click.Path, sequence_data_type: SequenceDataType, additional_params: t.Optional[t.Dict] = None) -> t.Union[str, str, str, str]:
    """
    :param sequence_data_path: unaligned sequence data
    :param sequence_data_type: sequence data type
    :param additional_params: additional program parameters, if needed
    :return: path to job completion validator file
    """
    program_name = "paml" if sequence_data_type == SequenceDataType.CODON else "phyml" # program will provide the parameters to simulate with
    # align the data
    working_dir = os.path.dirname(sequence_data_path)
    alignment_path = f"{working_dir}/aligned_sequences.fasta"
    BaseTools.align(input_path=sequence_data_path, output_path=alignment_path, sequence_data_type=sequence_data_type, alignment_method=AlignmentMethod.MAFFT)
    output_path = working_dir if program_name == "phyml" else f"{working_dir}/paml.out"

    # create a program instance
    program_to_exec = program_to_callable[program_name]()

    # run the inference program
    completion_validator_path = program_to_exec.exec(
        input_path=alignment_path,
        output_path=output_path,
        aux_dir=working_dir,
        additional_params=additional_params,
        parallelize=True,
        cluster_data_dir=os.path.dirname(alignment_path),
        priority=0,
        queue="itaym",
        wait_until_complete=True,
        get_completion_validator=True,
        control_file_path=f"{working_dir}/input.ctl",
        input_tree_path=f"{working_dir}/tree.nwk"
    )
    return alignment_path, output_path, working_dir, completion_validator_path


@click.command()
@click.option("--sequence_data_path",
              help="path to the full sequence data in a fasta format",
              type=click.Path(exists=True, file_okay=True, readable=True),
              required=True)
@click.option("--sequence_data_type",
              help="type of sequence data. options are: nucleotide, codon, amino_acid",
              type=str,
              required=True)
@click.option("--required_data_size",
              help="integer of the required data size",
              type=int,
              required=False,
              default=500)
@click.option("--num_of_repeats",
              help="number of repeats of down sampling the data",
              type=int,
              required=False,
              default=1)
@click.option("--output_dir",
              help="directory to write the sampled data to",
              type=click.Path(exists=False, dir_okay=True),
              required=True)
@click.option("--additional_program_parameters",
              help="a path ot a json file with extra program parameters",
              default=None,
              required=False)
@click.option("--additional_simulation_parameters",
              help="path to json file with additional simulation parameters",
              default=None,
              required=False)
@click.option("--log_path",
              help="path ot log file",
              type=click.Path(exists=False, dir_okay=True),
              required=True)
def prepare_data(sequence_data_path: click.Path,
                 sequence_data_type: str, required_data_size: int,
                 num_of_repeats: int,
                 output_dir: click.Path,
                 additional_program_parameters: t.Optional[click.Path],
                 additional_simulation_parameters: t.Optional[click.Path],
                 log_path: click.Path):
    """reduced the given data to a required size by randomly sampling sequences without repeats.
       can repeat the procedure multiple times to augment data"""

    # intialize the logger
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s module: %(module)s function: %(funcName)s line: %(lineno)d %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_path),
        ],
    )
    logger = logging.getLogger(__name__)
    logger.info("Initiating simulations preparation from real data")

    sequence_data_type = SequenceDataType(sequence_data_type)
    os.makedirs(output_dir, exist_ok=True)
    full_data = list(SeqIO.parse(sequence_data_path, "fasta"))
    logger.info("Data loaded successfully")
    remove_duplicates(sequence_records=full_data)
    logger.info("Accession duplicates removed from data")
    sampled_data_paths = sample_data(full_data=full_data, output_dir=output_dir, required_data_size=required_data_size, num_of_repeats=num_of_repeats)
    logger.info(f"{num_of_repeats} samples of size {required_data_size} generated successfully")
    sample_to_output = dict()
    for path in sampled_data_paths:
        if additional_program_parameters and os.path.exists(additional_program_parameters):
            with open(additional_program_parameters, "r") as input_file:
                additional_program_parameters = json.load(input_file)
        logger.info(f"executing program on sample {path}")
        alignment_path, program_output_path, job_output_dir, completion_validator_path = run_program(sequence_data_path=path, sequence_data_type=sequence_data_type, additional_params=additional_program_parameters)
        sample_to_output[path] = {"alignment_path": alignment_path, "program_output_path": program_output_path, "job_output_dir": job_output_dir, "completion_validator_path": completion_validator_path}

    # wait for the program to finish
    for path in sample_to_output:
        while not os.path.exists(sample_to_output[path]["completion_validator_path"]):
            sleep(10)
    logger.info("execution of program on all samples is complete")

    # write simulations pipeline input according to the output of the program
    logger.info(f"parsing program output ro simulation inputs")
    program_name = "paml" if sequence_data_type == SequenceDataType.CODON else "phyml"
    program_to_exec = program_to_callable[program_name]()
    for path in sample_to_output:
        output = program_to_exec.parse_output(output_path=sample_to_output[path]["program_output_path"], job_output_dir=sample_to_output[path]["job_output_dir"])
        if not additional_simulation_parameters:
            additional_simulation_parameters = dict()
        if not "simulations_output_dir" in additional_simulation_parameters:
            additional_simulation_parameters["simulations_output_dir"] = f"{output_dir}/simulations/"
        if not "sequence_data_type" in additional_simulation_parameters:
            additional_simulation_parameters["sequence_data_type"] = sequence_data_type
        if not "seq_len" in additional_simulation_parameters:
            additional_simulation_parameters["seq_len"] = len(list(SeqIO.parse(sample_to_output[path]["alignment_path"])[0]))
        if not "ntaxa" in additional_simulation_parameters:
            additional_simulation_parameters["ntaxa"] = len(full_data)
        program_to_exec.write_output_to_simulation_pipeline_json(program_output=output, output_path=f"{os.path.dirname(path)}/simulations.json", additional_simulation_parameters=additional_simulation_parameters)
        logger.info(f"parsing complete")

if __name__ == '__main__':
    prepare_data()
