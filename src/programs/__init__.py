from .rate4site import Rate4Site
from .paml import Paml
from .program import Program
from enum import Enum


class ProgramName(Enum):
    RATE4SITE = "rate4site"
    PAML = "paml"


program_to_callable = {
    ProgramName.RATE4SITE.value: Rate4Site,
    ProgramName.PAML.value: Paml,
}
