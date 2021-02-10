from .rate4site import Rate4Site
from .paml import PAML
from .program import Program
from enum import Enum


class ProgramName(Enum):
    RATE4SITE = "rate4site"


program_to_callable = {
    ProgramName.RATE4SITE.value: Rate4Site,
}
