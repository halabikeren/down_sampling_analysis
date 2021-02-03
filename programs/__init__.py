from programs.program import Program
from programs.rate4site import Rate4Site
from enum import Enum


class ProgramName(Enum):
    RATE4SITE = "rate4site"


program_to_callable = {
    ProgramName.RATE4SITE.value: Rate4Site,
}

__all__ = ["program_to_callable", "ProgramName"] + [
    c.__name__ for c in program_to_callable.values()
]
