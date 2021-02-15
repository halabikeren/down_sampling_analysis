import subprocess
from Bio import SeqIO
import typing as t
from dataclasses import dataclass, field
import os
import logging

from ete3 import Tree

from .sampler import Sampler

log = logging.getLogger(__name__)


@dataclass
class CdHit(Sampler):
    thr_to_sample: t.Dict[float, t.List[SeqIO.SeqRecord]] = field(default_factory=dict)

    def __init__(self, sequence_data_path: str, tree_path: str, sequences: t.Optional[t.List[SeqIO.SeqRecord]] = None):
        super(CdHit, self).__init__(sequence_data_path=sequence_data_path, tree_path=tree_path, sequences=sequences)
        self.thr_to_sample = dict()

    def run_cd_hit(self, threshold: float, aux_dir: str) -> str:
        """
        executes cd-hit on the sequence data of the instance given a threshold
        :param threshold: similarity threshold for ch-hit
        :param aux_dir directory to generate auxiliary files in
        :return: path of the output file provided by cd-hit
        """
        input_name = os.path.splitext(os.path.basename(self.sequences_path))[0]
        records = list((SeqIO.parse(self.sequences_path, "fasta")))
        for record in records:
            if "-" in record.seq:
                raise ValueError(
                    f"sequence data at {self.sequences_path}is aligned, so cd-hit cannot run on it"
                )
        output_file = f"{aux_dir}/{input_name}_threshold_{threshold}"
        word_len = (
            (5 if threshold > 0.7 else 4)
            if threshold > 0.6
            else (3 if threshold > 0.5 else 2)
        )
        cmd = f"cd-hit -i {self.sequences_path} -o {output_file} -c {threshold} -n {word_len}"
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if len(process.stderr.read()) > 0:
            raise RuntimeError(
                f"CD-HIT failed to properly execute and provide an output file with error {process.stderr.read()} and output is {process.stdout.read()}"
            )
        return output_file

    def get_cdhit_output(
        self, threshold: float, aux_dir: str
    ) -> t.List[SeqIO.SeqRecord]:
        """
        :param threshold: similarity threshold for cd-hit
        :param aux_dir directory to put auxiliary files in
        :return: a list of members in the sample of cd-hit according to the given threshold
        """
        output_path = self.run_cd_hit(threshold, aux_dir)
        chosen_subset = list(SeqIO.parse(output_path, "fasta"))
        return chosen_subset

    def get_similarity_threshold(
        self, k: int, left_thr: float, right_thr: float, aux_dir: str
    ) -> float:
        """
        recursive function for searching the similarity threshold that, when given to cd-hit, would produce a sample of
        the required size that is held in self.sample_size
        The function fills the threshold matching the sample size,
        or the one that it as close as can be obtained, and changes the sample size accordingly
        :param k required sample size
        :param left_thr: smallest threshold to search from
        :param right_thr: largest threshold to search up to
        :param aux_dir directory to put auxiliary files in
        """
        if left_thr not in self.thr_to_sample:
            self.thr_to_sample[left_thr] = self.get_cdhit_output(left_thr, aux_dir)
            self.thr_to_sample[left_thr] = self.get_cdhit_output(left_thr, aux_dir)
        left_sample_size = len(self.thr_to_sample[left_thr])

        if right_thr not in self.thr_to_sample:
            self.thr_to_sample[right_thr] = self.get_cdhit_output(right_thr, aux_dir)
        right_sample_size = len(self.thr_to_sample[right_thr])

        # stop conditions
        if k < left_sample_size:
            log.info(
                f"there is no convergence of thresholds for the clusters number {k}"
            )
            log.info(
                f"the closest clusters' number to the required one is {left_sample_size}"
            )
            return left_thr
        if (
            (right_thr - left_thr) < 0.01
            or left_thr >= right_thr
            or k > right_sample_size
        ):
            log.info(
                f"there is no convergence of thresholds for the clusters number {k}"
            )
            log.info(
                f"the closest clusters' number to the required one is {right_sample_size}"
            )
            return right_thr
        if k == right_sample_size:
            return right_thr
        if k == left_sample_size:
            return left_thr

        # recursive search
        middle_thr = (
            (right_thr - left_thr)
            / ((right_sample_size - left_sample_size) / (k - left_sample_size))
        ) % 1 + left_thr

        if middle_thr not in self.thr_to_sample:
            self.thr_to_sample[middle_thr] = self.get_cdhit_output(middle_thr, aux_dir)
        middle_sample_size = len(self.thr_to_sample[middle_thr])

        if k == middle_sample_size:
            return middle_thr

        elif k > middle_sample_size:
            return self.get_similarity_threshold(k, middle_thr, right_thr, aux_dir)

        else:  # self.sample_size < middle_sample_size:
            return self.get_similarity_threshold(k, left_thr, middle_thr, aux_dir)

    def get_sample(
        self, k: int, aux_dir: str, **kwargs
    ) -> t.List[SeqIO.SeqRecord]:
        """
        :param k: number of sequences to sample
        :param aux_dir directory to generate auxiliary files in
        :return: either a path to the generated sample or a list of samples sequence names
        """
        sample = super(CdHit, self).get_sample(k, aux_dir)
        if 1 < k < len(self.sequences):
            os.makedirs(aux_dir, exist_ok=True)
            thr = self.get_similarity_threshold(
                k-1, 0.4, 1, aux_dir
            )  # the minimum threshold 0.4 is set based on experience
            sample = self.thr_to_sample[thr]
            sample.append(self.saved_sequence)  # saved sequence must be in every sample
        return sample
