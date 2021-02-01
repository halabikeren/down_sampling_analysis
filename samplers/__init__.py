from enum import Enum
from samplers.Sampler import Sampler
from samplers.Pda import Pda
from samplers.CdHit import CdHit
from samplers.Random import Random


class SamplingMethod(Enum):
    Pda = "pda"
    CdHit = "cdhit"
    Random = "random"


__all__ = [en.value for en in SamplingMethod] + ["SamplingMethod"]
