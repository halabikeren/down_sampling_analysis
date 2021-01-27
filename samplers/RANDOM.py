from Bio import SeqIO
import typing as t
from pydantic import BaseModel
import random
import os
import logging

log = logging.getLogger(__name__)


class RANDOM(BaseModel):
    all_sequences_path: str
    sample_size: int = 0
    sample_members: t.List[str] = []
    sampled_sequences_path: str = None

    def compute_sample(self, k: int, write: bool = False):
        sample_candidates = list(SeqIO.parse(self.all_sequences_path, "fasta"))
        if k < 0 or k > len(sample_candidates):
            raise ValueError(f"sample size {k} is invalid for data of size {len(sample_candidates)}")
        self.sample_size = k
        if k == len(sample_candidates):
            self.sample_members = [record.description for record in sample_candidates]
            if write:
                res = os.system(f"cp -r {self.all_sequences_path} {self.sampled_sequences_path}")
        else:
            sampled_records = random.sample(sample_candidates, k)
            self.sample_members = [record.description for record in sampled_records]
            if write:
                SeqIO.write(sampled_records, self.sampled_sequences_path, "fasta")
