from enum import Enum
from functools import partial

from samplers.sampler import Sampler
from samplers.pda import Pda
from samplers.cdhit import CdHit
from samplers.random import Random


class SamplingMethod(Enum):
    PDA = ("pda", partial(Pda))
    CDHIT = ("cdhit", partial(CdHit))
    RANDOM = ("random", partial(Random))

    def __call__(self):
        return self._call()


SamplingMethod.PDA._call = Pda
SamplingMethod.CDHIT._call = CdHit
SamplingMethod.RANDOM._call = Random


__all__ = [en.value for en in SamplingMethod] + ["SamplingMethod"]
