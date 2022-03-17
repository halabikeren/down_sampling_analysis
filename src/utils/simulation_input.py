import typing as t
import os
from pydantic import BaseModel, FilePath, Field, validator
from .types import SequenceDataType, AlignmentMethod, TreeReconstructionMethod
from samplers import SamplingMethod

import logging

logger = logging.getLogger(__name__)


class SimulationInput(BaseModel):
    simulations_output_dir: str  # path to write simulations to
    sequence_data_type: SequenceDataType  # the type of provided sequence data in the form of SequenceDataType
    substitution_model: str  # substitution model. will be build as a enum later
    substitution_model_params: t.Dict[str , t.Union[float,str]]  # maps tuples of two states to the rate of substitution between them
    states_frequencies: t.Dict[str, float]  # maps state (character / triplet of characters in case of codons) to their frequencies
    tree_rooted: bool = True
    tree_random: bool = True
    tree_length: t.Optional[float] = None
    simulation_tree_path: t.Optional[str] = None
    ntaxa: int
    seq_len: int = 1000
    birth_rate: float = 0.3
    death_rate: float = 0.1
    sample_rate: float = 0
    mutation_rate: float = 0.5
    pinv: float = 0.1,
    alpha: float = 0.5
    ngamcat: int = 16
    nrep: int = 30
    use_simulated_alignment: bool = True
    use_simulated_tree: bool = True
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
        default_factory=lambda: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
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
    programs: t.List[str] = Field(
        default_factory=lambda: ["rate4site"]
    )  # the programs that should be executed on the sampled data in the scope of the pipeline
    programs_params: t.Optional[
        t.Dict[str, t.Any]
    ] = None  # a map of programs to parameters it should be executed with
    use_full_alignment_in_sample: bool = False  # indicates weather the full alignment should be trimmed ot create an alignment for the sampled, or if alignment should be reconstructed for the sample from scratch
    use_full_tree_in_sample: bool = False  # indicates weather the full alignment should be trimmed ot create an alignment for the sampled, or if alignment should be reconstructed for the sample from scratch
    exec_on_full_data: bool = (
        True  # indicates weather the program should be executed on the full dataset
    )
    reference_data_paths: t.Optional[
        t.Dict[str, FilePath]] = None  # maps a program to its relevant simulated reference data
    parallelize: bool = True  # indicator weather execution of programs on the samples should be parallelized or not
    cluster_data_dir: t.Optional[str] = None
    priority: int = (
        0  # in case of parallelization, this parameter sets the priority of the jobs
    )
    queue: str = "itaym"  # in case of parallelization, this parameter sets the queue that jobs should be submitted to

    def __init__(self, **data: t.Any):
        for field in [
            "simulations_output_dir"
            "pipeline_dir",
            "unaligned_sequence_data_path",
            "aligned_sequence_data_path",
            "tree_path",
        ]:
            if field in data:
                data[field] = os.path.join(os.getcwd(), data[field])
        super().__init__(**data)

    @validator("simulation_tree_path")
    def given_if_not_random(cls, v, values):
        if not v and not values["tree_random"]:
            raise ValueError(
                "If tree for simulations is not set to be random, it must be provided in simulation_tree_path"
            )
        return v

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
