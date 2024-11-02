from dataclasses import dataclass
from dataclasses import field as dataclass_field
from enum import IntEnum
from pathlib import Path

import torch


class Fields(IntEnum):
    NDIM = 1
    PRIM = 2
    SUPER = 3
    ORB = 4
    HAMILT = 5
    SYMM = 6
    PHASE = 7
    BONDS = 8
    PAIRS = 9
    DILUT = 10


def get_pi() -> float:
    """Get value of PI using torch"""
    return torch.acos(torch.tensor(-1.0))


@dataclass
class GeomParams:
    RDIM: int = 3
    TOLL: float = 1e-6
    PI: float = dataclass_field(default_factory=get_pi)
    IM: complex = 1j

    # Input fields
    INPUT_FIELDS: list[str] = (
        "#NDIM",
        "#PRIM",
        "#SUPER",
        "#ORB",
        "#HAMILT",
        "#SYMM",
        "#PHASE",
        "#BONDS",
        "#PAIR",
        "#DILUT",
    )

    def __post_init__(self):
        self.found_fields = {field_name: False for field_name in self.INPUT_FIELDS}

    def analyze_input(self, input_file: Path) -> None:
        """Analyze which fields are specified in the input file"""
        with open(input_file) as f:
            content = f.read()

        # Check which fields are present
        for field_name in self.INPUT_FIELDS:
            self.found_fields[field_name] = field_name in content

        # Check compulsory fields using list comprehension
        missing = [field_name for field_name in self.INPUT_FIELDS[:5] if not self.found_fields[field_name]]

        if missing:
            raise ValueError(f"Missing compulsory fields: {', '.join(missing)}")

        # Check PAIRS dependency on BONDS
        if self.found_fields["#PAIR"] and not self.found_fields["#BONDS"]:
            raise ValueError("#PAIR requires #BONDS to be specified in input")
