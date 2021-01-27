from Bio import SeqIO
import typing as t
from pydantic import BaseModel
import os
import logging
log = logging.getLogger(__name__)

class CDHIT(BaseModel):
    sequences_path: str
    aux_dir: str = os.getcwd()
    sample_size: int = 0
    similarity_threshold: float = 0
    sample_members: t.List[str] = []
    thr_to_sample: t.Dict[float, t.List[str]] = dict()

    class Config:
        arbitrary_types_allowed = True

    def run_cd_hit(self, threshold: float) -> str:
        """
        executes cd-hit on the sequence data of the instance given a threshold
        :param threshold: similarity threshold for ch-hit
        :return: path of the output file provided by cd-hit
        """
        input_name = os.path.splitext(os.path.basename(self.sequences_path))[0]
        output_file = f"{self.aux_dir}/{input_name}_threshold_{threshold}"
        word_len = (5 if threshold > 0.7 else 4) if threshold > 0.6 else (3 if threshold > 0.5 else 2)
        res = os.system(f"cd-hit -i {self.sequences_path} -o {output_file} -c {threshold} -n {word_len}")
        if res != 0:
            raise IOError("CD-HIT failed to properly execute and provide an output file")
        return output_file

    def get_cdhit_output(self, threshold: float) -> t.List[str]:
        """
        :param threshold: similarity threshold for cd-hit
        :return: a list of members in the sample of cd-hit according to the given threshold
        """
        output_path = self.run_cd_hit(threshold)
        chosen_subset = [record.description for record in SeqIO.parse(output_path, "fasta")]
        return chosen_subset

    def compute_similarity_threshold(self, left_thr: float, right_thr: float):
        """
        recursive function for searching the similarity threshold that, when given to cd-hit, would produce a sample of
        the required size that is held in self.sample_size
        The function fills the threshold matching the sample size,
        or the one that it as close as can be obtained, and changes the sample size accordingly
        :param left_thr: smallest threshold to search from
        :param right_thr: largest threshold to search up to
        """
        if left_thr not in self.thr_to_sample:
            self.thr_to_sample[left_thr] = self.get_cdhit_output(left_thr)
        left_sample_size = len(self.thr_to_sample[left_thr])

        if right_thr not in self.thr_to_sample:
            self.thr_to_sample[right_thr] = self.get_cdhit_output(right_thr)
        right_sample_size = len(self.thr_to_sample[right_thr])

        # stop conditions
        if self.sample_size < left_sample_size:
            log.info(f"there is no convergence of thresholds for the clusters number {self.sample_size}")
            log.info(f"the closest clusters' number to the required one is {left_sample_size}")
            self.sample_size = left_sample_size
            self.similarity_threshold = left_thr
            self.sample_members = self.thr_to_sample[left_thr]
            return
        if (right_thr - left_thr) < 0.01 or left_thr >= right_thr or self.sample_size > right_sample_size:
            log.info(f"there is no convergence of thresholds for the clusters number {self.sample_size}")
            log.info(f"the closest clusters' number to the required one is {right_sample_size}")
            self.sample_size = right_sample_size
            self.similarity_threshold = right_thr
            self.sample_members = self.thr_to_sample[right_thr]
            return
        if self.sample_size == right_sample_size:
            self.sample_size = right_sample_size
            self.similarity_threshold = right_thr
            self.sample_members = self.thr_to_sample[right_thr]
            return
        if self.sample_size == left_sample_size:
            self.sample_size = left_sample_size
            self.similarity_threshold = left_thr
            self.sample_members = self.thr_to_sample[left_thr]
            return

        # recursive search
        middle_thr = ((right_thr - left_thr) / ((right_sample_size - left_sample_size) / (
                self.sample_size - left_sample_size))) % 1 + left_thr

        if middle_thr not in self.thr_to_sample:
            self.thr_to_sample[middle_thr] = self.get_cdhit_output(middle_thr)
        middle_sample_size = len(self.thr_to_sample[middle_thr])

        if self.sample_size == middle_sample_size:
            self.similarity_threshold = middle_thr
            self.sample_members = self.thr_to_sample[middle_thr]
            return

        elif self.sample_size > middle_sample_size:
            self.compute_similarity_threshold(middle_thr, right_thr)

        else:  # self.sample_size < middle_sample_size:
            self.compute_similarity_threshold(left_thr, middle_thr)

    def compute_sample(self, k: int) -> t.List[str]:
        """
        :param k: the required sample size
        :return: the sampled sequence names
        """
        res = os.system(f"mkdir -p {self.aux_dir}")
        self.sample_size = k
        self.compute_similarity_threshold(0.4, 1)  # the minimum threshold is set based on experience
        self.sample_members = self.thr_to_sample[self.similarity_threshold]
        return self.sample_members
