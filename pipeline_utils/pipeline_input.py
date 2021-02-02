import typing as t
from pipeline_utils import SequenceDataType, AlignmentMethod, TreeReconstructionMethod
from pydantic import BaseModel, FilePath, DirectoryPath, Field
from samplers import SamplingMethod

import logging

logger = logging.getLogger(__name__)


class PipelineInput(BaseModel):
    pipeline_dir: DirectoryPath  # path in which the entire pipeline input and output will be placed
    sequence_data_path: FilePath  # full path to a fasta file with sequence data that is either aligned or not aligned
    sequence_data_type: SequenceDataType  # the type of provided sequence data in the form of SequenceDataType
    alignment_method: AlignmentMethod = (
        AlignmentMethod.MAFFT
    )  # the method that the alignment should be build with, in case that it is not provided
    alignment_params: t.Optional[
        t.Dict[str, str]
    ] = None  # the parameters that the alignment method should be executed with. The default ones are viewable in
    # the method Pipeline.align()
    tree_path: t.Optional[
        str
    ] = None  # a full path to a tree file in newick format. Should be provided only if the user wishes PDA to use it
    # during sampling
    tree_reconstruction_method: TreeReconstructionMethod = (
        TreeReconstructionMethod.UPGMA
    )  # the method that the tree should be build with. Can be either UPGMA, NJ or ML (which uses RaxML)
    tree_reconstruction_params: t.Optional[
        t.Dict[str, str]
    ] = None  # the parameters that the tree reconstruction method should be executed with. The default ones are
    # available in the method Pipeline.build_tree()
    sampling_fractions: t.List[float] = Field(
        default_factory=lambda: [0.25, 0.5, 0.75, 1]
    )  # the fractions of sampling that should be generated in the scope of the pipeline
    sampling_methods: t.List[SamplingMethod] = Field(
        default_factory=lambda: [m for m in SamplingMethod]
    )  # the methods of sampling that should be used in the scope of the pipeline
    programs: t.List[str] = Field(
        default_factory=lambda: ["r4s"]
    )  # the programs that should be executed on the sampled data in the scope of the pipeline
    program_params: t.Optional[
        t.Dict[str, t.Any]
    ] = None  # a map of programs to parameters it should be executed with
    weight_pda: bool = (
        False  # indicator weather when using PDA, weighting should be used or not
    )
    parallelize: bool = True  # indicator weather execution of programs on the samples should be parallelized or not
    priority: int = (
        0  # in case of parallelization, this parameter sets the priority of the jobs
    )
    queue: str = "itaym"  # in case of parallelization, this parameter sets the queue that jobs should be submitted to
