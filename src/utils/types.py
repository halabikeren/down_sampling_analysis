from enum import Enum


class SequenceDataType(Enum):
    NUC = "nucleotide"
    AA = "amino_acid"
    CODON = "codon"


class AlignmentMethod(Enum):
    MAFFT = "mafft"
    PRANK = "prank"


class TreeReconstructionMethod(Enum):
    UPGMA = "upgma"
    FASTTREE = "fasttree"
    NJ = "nj"
    ML = "ml"


class Queue(Enum):
    ITAYM = "itaym"
    ITAYMAA = "itaymaa"
    ITAYM1 = "itaym1"
    ITAYM2 = "itaym2"
    ITAYM3 = "itaym3"
    ITAYM4 = "itaym4"
