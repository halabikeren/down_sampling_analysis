import typing as t

import os

from utils import (
    SequenceDataType,
    AlignmentMethod,
    TreeReconstructionMethod,
)
from programs import ProgramName
from pydantic import BaseModel, FilePath, DirectoryPath, Field, validator
from samplers import SamplingMethod

import logging

logger = logging.getLogger(__name__)


class PipelineInput(BaseModel):
    pipeline_dir: DirectoryPath  # path in which the entire pipeline input and output will be placed
    unaligned_sequence_data_path: FilePath  # full path to a fasta file with sequence data that is either aligned or
    # not aligned
    aligned_sequence_data_path: t.Optional[
        FilePath
    ] = None  # full path to a fasta file with sequence data that is either aligned or not aligned
    sequence_data_type: SequenceDataType  # the type of provided sequence data in the form of SequenceDataType
    alignment_method: AlignmentMethod = (
        AlignmentMethod.MAFFT
    )  # the method that the alignment should be build with, in case that it is not provided
    alignment_params: t.Optional[
        t.Dict[str, t.Any]
    ] = None  # the parameters that the alignment method should be executed with. The default ones are viewable in #
    # the method Pipeline.align()
    tree_path: t.Optional[
        FilePath
    ] = None  # a full path to a tree file in newick format. Should be provided only if the user wishes PDA to use it
    # during sampling
    tree_reconstruction_method: TreeReconstructionMethod = (
        TreeReconstructionMethod.FASTTREE
    )  # the method that the tree should be build with. Can be either UPGMA, NJ or ML (which uses RaxML)
    tree_reconstruction_params: t.Optional[
        t.Dict[str, str]
    ] = None  # the parameters that the tree reconstruction method should be executed with. The default ones are
    # available in the method Pipeline.build_tree()
    sampling_fractions: t.List[float] = Field(
        default_factory=lambda: [0.25, 0.5, 0.75]
    )  # the fractions of sampling that should be generated in the scope of the pipeline (must be > 0 and < 1)
    sampling_methods: t.List[SamplingMethod] = Field(
        default_factory=lambda: [m for m in SamplingMethod]
    )  # the methods of sampling that should be used in the scope of the pipeline
    weight_pda: bool = (
        False  # indicator weather when using PDA, weighting should be used or not
    )
    use_external_pda: bool = True  # indicator weather external pda should be used
    samples_alignment_method: AlignmentMethod = (
        AlignmentMethod.MAFFT
    )  # the method that the alignment should be build with, in case that it is not provided
    samples_alignment_params: t.Optional[t.Dict[str, t.Any]] = None
    programs: t.List[ProgramName] = Field(
        default_factory=lambda: [n for n in ProgramName]
    )  # the programs that should be executed on the sampled data in the scope of the pipeline
    programs_params: t.Optional[
        t.Dict[str, t.Any]
    ] = None  # a map of programs to parameters it should be executed with
    exclude_ref_seq: bool = False
    use_full_alignment_in_sample: bool = False  # indicates weather the full alignment should be trimmed ot create an alignment for the sampled, or if alignment should be reconstructed for the sample from scratch
    use_full_tree_in_sample: bool = False  # indicates weather the full alignment should be trimmed ot create an alignment for the sampled, or if alignment should be reconstructed for the sample from scratch
    exec_on_full_data: bool = (
        True  # indicates weather the program should be executed on the full dataset
    )
    reference_data_paths: t.Optional[
        t.Dict[str, t.Dict[str, FilePath]]
    ] = None  # maps a program to its relevant simulated reference data
    parallelize: bool = True  # indicator weather execution of programs on the samples should be parallelized or not
    cluster_data_dir: t.Optional[str] = None
    priority: int = (
        0  # in case of parallelization, this parameter sets the priority of the jobs
    )
    queue: str = "itaym"  # in case of parallelization, this parameter sets the queue that jobs should be submitted to

    def __init__(self, **data: t.Any):
        for field in [
            "pipeline_dir",
            "unaligned_sequence_data_path",
            "aligned_sequence_data_path",
            "tree_path",
        ]:
            if field in data and data[field]:
                data[field] = os.path.join(os.getcwd(), data[field])
        super().__init__(**data)

    @validator("cluster_data_dir")
    def given_if_parallelize(cls, v, values):
        if not v and values["parallelize"]:
            raise ValueError(
                "Cannot set parallelization without providing cluster data dir"
            )
        return v

    @validator("sampling_fractions")
    def between_zero_and_one(cls, v):
        for item in v:
            if item <= 0 or item >= 1:
                raise ValueError(
                    f"Sampling fraction {item} is invalid. A value must be between 0 and 1, excluded"
                )
        return v
