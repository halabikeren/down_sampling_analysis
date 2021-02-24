import shutil
import typing as t
from Bio import SeqIO
from dataclasses import dataclass

from ete3 import Tree


@dataclass
class Sampler:
    sequences_path: str
    tree: Tree
    sequences: t.List[SeqIO.SeqRecord]
    exclude_a_ref_sequence: bool = False
    saved_sequence: t.Optional[SeqIO.SeqRecord] = None  # this sequence will appear in every sample

    def __init__(self, sequence_data_path: str, tree_path: str, exclude_a_ref_sequence: bool = False, sequences: t.Optional[t.List[SeqIO.SeqRecord]] = None):
        if sequences:
            self.sequences = sequences
        else:
            self.sequences = list(SeqIO.parse(sequence_data_path, "fasta"))
        self.exclude_a_ref_sequence = exclude_a_ref_sequence
        if self.exclude_a_ref_sequence:
            self.saved_sequence = self.sequences[0]
            self.sequences.remove(self.saved_sequence)
            self.sequences_path = sequence_data_path.replace(".fas", f"_without_reference_{self.saved_sequence.name}.fas")
        else:
            self.sequences_path = sequence_data_path
        SeqIO.write(self.sequences, self.sequences_path, "fasta")
        self.tree = Tree(tree_path, format=1)

    def get_sample(
        self, k: int, aux_dir: str, **kwargs
    ) -> t.List[SeqIO.SeqRecord]:
        """
        :param k: number of sequences to sample
        :param aux_dir directory to generate auxiliary files in
        :return: either a path to the generated sample or a list of samples sequence names
        """
        if k < 0 or k > len(self.sequences)+1:
            raise ValueError(
                f"sample size {k} is invalid for data of size {len(self.sequences)+1}"
            )
        if k == 1 and self.exclude_a_ref_sequence:
            return [self.saved_sequence]
        sample = self.sequences
        if self.exclude_a_ref_sequence:
            sample += [self.saved_sequence]
        return sample

    def write_sample(self, k: int, aux_dir: str, output_path: str, **kwargs) -> t.Union[str, t.List[SeqIO.SeqRecord]]:
        """
        writes the sampled data to an output file
        :param k required sample size
        :param aux_dir directory to generate auxiliary files in
        :param output_path path to write the sample output to
        :return: None
        """
        sample = self.get_sample(k, aux_dir, **kwargs)
        SeqIO.write(sample, output_path, "fasta")
        return sample
