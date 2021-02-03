from enum import Enum
from samplers.sampler import Sampler
from samplers.pda import Pda
from samplers.cdhit import CdHit
from samplers.random import Random


class SamplingMethod(Enum):
    PDA = "pda"
    CDHIT = "cdhit"
    RANDOM = "random"


method_to_callable = {
    SamplingMethod.PDA.value: Pda,
    SamplingMethod.CDHIT.value: CdHit,
    SamplingMethod.RANDOM.value: Random,
}

__all__ = ["method_to_callable", "SamplingMethod"] + [
    c.__name__ for c in method_to_callable.values()
]
