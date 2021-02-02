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
    NJ = "nj"
    ML = "ml"
