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
from programs import ProgramName, program_to_callable
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


def sample_data(full_data: t.List[SeqIO.SeqRecord], output_dir: str, required_data_size: int, num_of_repeats: int) -> \
        t.List[str]:
    if required_data_size > len(full_data):
        raise ValueError(f"full data size {len(full_data)} cannot be down-sampled to {required_data_size}")
    os.makedirs(output_dir, exist_ok=True)
    samples_paths = []
    for i in range(1, num_of_repeats + 1):
        os.makedirs(f"{output_dir}/sample_{i}/", exist_ok=True)
        output_path = f"{output_dir}/sample_{i}/sequences_old_names.fasta"
        if not os.path.exists(output_path):
            sampled_data = sample(full_data, required_data_size)
            SeqIO.write(sampled_data, output_path, "fasta")
        output_path = f"{output_dir}/sample_{i}/sequences_new_names.fasta"
        if not os.path.exists(output_path):
            convert_names(translator_path=f"{output_dir}/sample_{i}/new_to_old_names_map.pickle", records=sampled_data)
            output_path = f"{output_dir}/sample_{i}/sequences_new_names.fasta"
            SeqIO.write(sampled_data, output_path, "fasta")
        samples_paths.append(output_path)
    return samples_paths


def run_program(program_name: ProgramName, sequence_data_path: click.Path, alignment_path: str, sequence_data_type: SequenceDataType, program_output_path: str,
                additional_params: t.Optional[t.Dict] = None) -> t.Union[str, str, str, str]:
    """
    :param program_name: name of program to execute
    :param sequence_data_path: unaligned sequence data
    :param alignment_path: path in which an alignment will be created
    :param program_output_path: path to which the program output will be written
    :param sequence_data_type: sequence data type
    :param additional_params: additional program parameters, if needed
    :return: path to job completion validator file
    """

    # align the data
    BaseTools.align(input_path=sequence_data_path, output_path=alignment_path,
                    sequence_data_type=sequence_data_type,
                    alignment_method=AlignmentMethod.MAFFT)

    # create a program instance
    program_to_exec = program_to_callable[program_name]()

    # run the inference program (in the case of paml, the control file will be generated in the default directory
    completion_validator_path = program_to_exec.exec(
        input_path=alignment_path,
        output_path=program_output_path,
        aux_dir=f"{os.path.dirname(sequence_data_path)}/",
        additional_params=additional_params,
        parallelize=True,
        cluster_data_dir=os.path.dirname(alignment_path),
        priority=0,
        queue="itaym",
        wait_until_complete=False,
        get_completion_validator=True
    )

    return completion_validator_path


def output_exists(program_name: str, output_dir: str) -> bool:
    """
    :param program_name: name of program
    :param output_dir: directory that should bol its output
    :return: boolean indicating weather output already exists or not
    """
    if program_name == "paml" and "paml.out" in output_dir:
        if os.path.exists(output_dir):
            return True
        return False

    for path in os.listdir(output_dir):
        if program_name == "busted" and "BUSTED" in path:
            return True
        elif program_name == "phyml" and "phyml_stats" in path:
            return True
    return False


@click.command()
@click.option("--sequence_data_path",
              help="path to the full sequence data in a fasta format",
              type=click.Path(exists=True, file_okay=True, readable=True),
              required=True)
@click.option("--sequence_data_type",
              help="type of sequence data",
              type=click.Choice(['nucleotide', 'codon', 'amino_acid']),
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
@click.option("--program_name",
              help="name of the program to infect simulation parameters from the provided data with",
              type=click.Choice(['paml', 'phyml', 'busted']),
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
                 sequence_data_type: str,
                 required_data_size: int,
                 num_of_repeats: int,
                 output_dir: click.Path,
                 program_name: str,
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
    sampled_data_paths = sample_data(full_data=full_data, output_dir=output_dir,
                                     required_data_size=required_data_size,
                                     num_of_repeats=num_of_repeats)
    logger.info(f"{num_of_repeats} samples of size {required_data_size} generated successfully")
    sample_to_output = dict()
    for input_path in sampled_data_paths:
        if additional_program_parameters and os.path.exists(additional_program_parameters):
            with open(additional_program_parameters, "r") as input_file:
                additional_program_parameters = json.load(input_file)
        working_dir = f"{os.path.dirname(input_path)}/"
        program_output_path = working_dir if program_name == "phyml" or program_name == "busted" else f"{working_dir}/paml.out"
        alignment_path = str(sequence_data_path).replace(".fas", "_aligned.fas")
        completion_validator_path = None
        if not output_exists(program_name=program_name, output_dir=program_output_path):
            logger.info(f"executing program on sample {input_path}")
            completion_validator_path = run_program(
                program_name=program_name, sequence_data_path=input_path, alignment_path=alignment_path, sequence_data_type=sequence_data_type,
                program_output_path=program_output_path, additional_params=additional_program_parameters)
        sample_to_output[input_path] = {"alignment_path": alignment_path, "program_output_path": program_output_path,
                                  "job_output_dir": os.path.dirname(sequence_data_path)}
        if completion_validator_path:
            sample_to_output[input_path]["completion_validator_path"] = completion_validator_path

    # wait for the program to finish
    for input_path in sample_to_output:
        if "completion_validator_path" in sample_to_output[input_path]:
            while not os.path.exists(sample_to_output[input_path]["completion_validator_path"]):
                sleep(10)
    logger.info("execution of program on all samples is complete")

    # write simulations pipeline input according to the output of the program
    logger.info(f"parsing program output ro simulation inputs")
    if additional_simulation_parameters and os.path.exists(additional_simulation_parameters):
        with open(additional_simulation_parameters, "r") as input_file:
            additional_simulation_parameters = json.load(input_file)
    else:
        additional_simulation_parameters = dict()
    program_to_exec = program_to_callable[program_name]()
    for input_path in sample_to_output:
        output = program_to_exec.parse_output(output_path=sample_to_output[input_path]["program_output_path"],
                                              job_output_dir=sample_to_output[input_path]["job_output_dir"])
        additional_simulation_parameters["simulations_output_dir"] = f"{os.path.dirname(input_path)}/simulations/"
        os.makedirs(additional_simulation_parameters["simulations_output_dir"], exist_ok=True)
        if not "sequence_data_type" in additional_simulation_parameters:
            additional_simulation_parameters["sequence_data_type"] = sequence_data_type
        if not "seq_len" in additional_simulation_parameters:
            additional_simulation_parameters["seq_len"] = len(
                list(SeqIO.parse(sample_to_output[input_path]["alignment_path"], "fasta"))[0].seq) // (3 if sequence_data_type == SequenceDataType.CODON else 1)
        if not "ntaxa" in additional_simulation_parameters:
            additional_simulation_parameters["ntaxa"] = len(full_data)
        program_to_exec.write_output_to_simulation_pipeline_json(program_output=output,
                                                                 output_path=f"{os.path.dirname(input_path)}/simulations.json",
                                                                 additional_simulation_parameters=additional_simulation_parameters)
        logger.info(f"parsing complete")


if __name__ == '__main__':
    prepare_data()
