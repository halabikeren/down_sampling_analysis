import subprocess
import typing as t
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timedelta, date, time
from decimal import Decimal
from enum import Enum
from pathlib import PosixPath
from types import GeneratorType

from ete3 import Tree
import os
from Bio import SeqIO
from Bio.Seq import Seq
from dotenv import load_dotenv, find_dotenv
from pydantic import BaseModel
from pydantic.json import isoformat

from .types import SequenceDataType, AlignmentMethod, TreeReconstructionMethod

load_dotenv(find_dotenv())

import logging

logger = logging.getLogger(__name__)


@dataclass
class BaseTools:
    JSONABLE_OK = (str, int, float, type(None))
    SEQUENCES = (list, set, frozenset, GeneratorType, tuple)
    ENCODERS_BY_TYPE = {
        datetime: isoformat,
        date: isoformat,
        time: isoformat,
        timedelta: lambda td: td.total_seconds(),
        set: list,
        frozenset: list,
        GeneratorType: list,
        bytes: lambda o: o.decode(),
        Decimal: float,
    }

    @staticmethod
    def scale_tree(tree: Tree, required_size: float):
        """
        :param tree: tree to rescale
        :param required_size: required tree size
        :return: none. rescaled the tree
        """
        orig_size = 0
        for node in tree.traverse():
            orig_size += node.dist
        rescaling_factor = required_size / orig_size
        for node in tree.traverse():
            node.dist *= rescaling_factor

    @staticmethod
    def simplify_names(
        input_path: str,
        output_path: str,
        names_translator: t.Optional[t.Dict[str, str]] = None,
    ) -> t.Optional[t.Dict[str, str]]:
        """
        :param input_path: path with the original sequence names
        :param output_path:  path to which the sequences with the new names will be written
        :param names_translator: translator of new to old names. if not provided, simple names will be generated and returned
        :return:
        """
        input_is_tree = False
        if ".nwk" in str(input_path):
            input_is_tree = True
        if not input_is_tree:
            seq_records = list(SeqIO.parse(input_path, "fasta"))
            if not names_translator:
                s = 1
                new_to_orig_name = dict()
                for record in seq_records:
                    new_to_orig_name[f"S{s}"] = record.description
                    record.description = record.id = record.name = f"S{s}"
                    s += 1
                SeqIO.write(seq_records, output_path, "fasta")
                return new_to_orig_name
            else:
                reversed_names_translator = {
                    names_translator[key]: key for key in names_translator
                }
                for record in seq_records:
                    record.description = (
                        record.name
                    ) = record.id = reversed_names_translator[record.description]
                SeqIO.write(seq_records, output_path, "fasta")
        else:
            with open(input_path, "r") as infile:
                tree_str = infile.read()
            tree = Tree(tree_str, format=1)
            tree_leaves = tree.get_leaves()
            if not names_translator:
                s = 1
                new_to_orig_name = dict()
                for leaf in tree_leaves:
                    new_to_orig_name[f"S{s}"] = leaf.name
                    leaf.name = f"S{s}"
                    s += 1
                tree.write(outfile=output_path, format=5)
                return new_to_orig_name
            else:
                reversed_names_translator = {
                    names_translator[key]: key for key in names_translator
                }
                for leaf in tree_leaves:
                    leaf.name = reversed_names_translator[leaf.name]
                tree.write(outfile=output_path, format=5)

    @staticmethod
    def translate(sequence_records: t.List[SeqIO.SeqRecord], output_path: str):
        """
        :param sequence_records: list of coding sequence records
        :param output_path: the output path to which the translated sequences corresponding to the coding sequence
        records will be written to in fasta format
        :return: None
        """
        aa_records = deepcopy(sequence_records)
        for aa_record in aa_records:
            aa_record.seq = aa_record.seq.translate(table=1)
        SeqIO.write(aa_records, output_path, "fasta")

    @staticmethod
    def reverse_translate(
        unaligned_codon_path: str, aligned_aa_path: str, aligned_codon_path: str
    ):
        """
        :param unaligned_codon_path: path to non-aligned coding sequences
        :param aligned_aa_path: path to the aligned corresponding translated sequences
        :param aligned_codon_path: path to write the aligned coding sequences to
        :return: None
        """
        unaligned_codon_records = list(SeqIO.parse(unaligned_codon_path, "fasta"))
        aligned_aa_records = list(SeqIO.parse(aligned_aa_path, "fasta"))
        aligned_codon_records = []
        for aligned_aa_record in aligned_aa_records:
            unaligned_codon_record = [
                sequence_record
                for sequence_record in unaligned_codon_records
                if sequence_record.description == aligned_aa_record.description
            ][0]
            aligned_codon_record = deepcopy(unaligned_codon_record)
            aligned_codon_record_seq = ""
            aa_index = 0
            codon_index = 0
            while aa_index < len(aligned_aa_record.seq):
                if aligned_aa_record.seq[aa_index] != "-":
                    aligned_codon_record_seq += str(
                        unaligned_codon_record.seq[
                            codon_index * 3 : codon_index * 3 + 3
                        ]
                    )
                    codon_index += 1
                else:
                    aligned_codon_record_seq += "---"
                aa_index += 1
            aligned_codon_record.seq = Seq(aligned_codon_record_seq)
            aligned_codon_records.append(aligned_codon_record)
        SeqIO.write(aligned_codon_records, aligned_codon_path, "fasta")

    @staticmethod
    def align(
        input_path: str,
        output_path: str,
        sequence_data_type: SequenceDataType,
        alignment_method: AlignmentMethod,
        alignment_params: t.Optional[t.Dict[str, t.Any]] = None,
    ):
        """
        Aligns sequence data according to given method
        input_path: unaligned sequence data path
        output_path: the path to write the alignment to
        alignment_method: the method that the sequence data should be aligned with
        alignment_params: the parameters that the method should be run with
        :return:
        """
        if os.path.exists(output_path) and (
            os.path.getsize(output_path) > os.path.getsize(input_path)
        ):
            logger.info(
                f"{output_path} already exists. The program will assume it is complete"
            )
            return

        sequence_records = list(SeqIO.parse(input_path, "fasta"))
        alignment_input_path = input_path
        alignment_output_path = output_path
        if (
            sequence_data_type == SequenceDataType.CODON
            and alignment_method == AlignmentMethod.MAFFT
        ):
            alignment_input_path = input_path.replace(".fasta", "_translated.fasta")
            BaseTools.translate(
                sequence_records,
                input_path.replace(".fasta", "_translated.fasta"),
            )
            alignment_output_path = output_path.replace(".fasta", "_translated.fasta")

        cmd = ""
        if not alignment_method or alignment_method == AlignmentMethod.MAFFT:
            cmd = f"(mafft --localpair --maxiterate 1000 {alignment_input_path} > {alignment_output_path}) > /dev/null 2>&1"
        elif alignment_method == AlignmentMethod.PRANK:
            cmd = f"(prank -d={alignment_input_path} -o={alignment_output_path} -f=fasta -support {'-codon' if sequence_data_type == SequenceDataType.CODON else ''} -iterate=100 -showtree) > /dev/null 2>&1"
        if alignment_params:
            cmd += " ".join(
                [
                    f"{param_name} {alignment_params[param_name]}"
                    for param_name in alignment_params
                ]
            )
        if os.path.exists(alignment_output_path) and (
            os.path.getsize(alignment_output_path)
            > os.path.getsize(alignment_input_path)
        ):
            logger.info(
                f"Temporary alignment {alignment_output_path} already exists and will not be recreated."
            )
        else:
            process = os.system(cmd)
            if process != 0:
                raise IOError(
                    f"failed to align {alignment_output_path} with {alignment_method.value}"
                )
        if (
            alignment_method == AlignmentMethod.MAFFT
            and sequence_data_type == SequenceDataType.CODON
        ):
            BaseTools.reverse_translate(
                input_path,
                alignment_output_path,
                output_path,
            )
            os.remove(alignment_input_path)
            os.remove(alignment_output_path)

    @staticmethod
    def build_tree(
        input_path: str,
        output_path: str,
        sequence_data_type: SequenceDataType,
        tree_reconstruction_method: TreeReconstructionMethod,
        tree_reconstruction_params: t.Optional[t.Dict[str, t.Any]] = None,
    ):
        """
        :param input_path path to aligned sequence data in a fasta format
        :param output_path path in which the tree should be written in newick format
        :param sequence_data_type either nucleotide, amino_acid or codon
        :param tree_reconstruction_method: enum representing the tree reconstruction method
        :param tree_reconstruction_params: map of parameter names to parameter values
        :return: None
        """
        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            logger.info(
                f"{output_path} exists. The program will assume it is complete."
            )
            return

        if tree_reconstruction_method in [
            TreeReconstructionMethod.UPGMA,
            TreeReconstructionMethod.NJ,
        ]:
            sequence_data_type_hyphy = (
                "1"
                if sequence_data_type in [SequenceDataType.NUC, SequenceDataType.AA]
                else "2\\n1"
            )
            dist_method_hyphy = (
                "4"
                if sequence_data_type in [SequenceDataType.NUC, SequenceDataType.AA]
                else "3\\n3"
            )
            cmd = f"printf '12\\n1\\n1\\n{sequence_data_type_hyphy}\\n{input_path}\\n1\\n{dist_method_hyphy}\\ny\\n{output_path}\\n' | hyphy"
            process = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            if len(process.stderr.read()) > 0:
                raise IOError(
                    f"failed to reconstruct {tree_reconstruction_method.value} tree with hyphy due to error {process.stderr.read()}. Execution output is {process.stdout.read()}"
                )

        elif tree_reconstruction_method == TreeReconstructionMethod.FASTTREE:
            cmd = f"(fasttree {input_path} > {output_path}) > /dev/null 2>&1"
            res = os.system(cmd)
            if res:
                raise IOError(
                    f"failed to reconstruct {tree_reconstruction_method.value} tree with fasttree program."
                )
            tree = Tree(output_path)
            tree.write(outfile=output_path, format=5)

        elif tree_reconstruction_method == TreeReconstructionMethod.ML:
            output_dir, output_filename = os.path.split(output_path)
            aux_dir = f"{output_dir}/raxml_aux/"
            os.mkdir(aux_dir)
            os.chdir(aux_dir)
            model = (
                tree_reconstruction_params["-m"]
                if tree_reconstruction_params and "-m" in tree_reconstruction_params
                else "GTRGAMMA"
            )
            num_of_categories = (
                tree_reconstruction_params["-n"]
                if tree_reconstruction_params and "-n" in tree_reconstruction_params
                else 4
            )
            cmd = f"raxmlHPC -s {input_path} -n out -m {model} -c {num_of_categories}"
            process = subprocess.Popen(
                cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            if len(process.stderr.read()) > 0:
                raise IOError(
                    f"failed to reconstruct ML tree with raxml due to error {process.stderr.read()}. Execution output is {process.stdout.read()}"
                )
            os.rename(f"{aux_dir}RAxML_bestTree.out", output_path)
            os.remove(aux_dir)

