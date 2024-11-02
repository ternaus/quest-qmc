# src/symmetry.py
from dataclasses import dataclass

import numpy as np
import torch

DEVICE = torch.get_default_device()


@dataclass
class SymmetryOperations:
    """Symmetry operations for lattice"""

    ntotsymm: int  # total number of symmetries from input
    nsymm: int  # number of symmetries compatible with supercell
    ntransl: int  # number of translations
    device: torch.device = DEVICE

    def __post_init__(self):
        # Symmetry mappings
        self.map_symm = None  # site mapping under symmetry
        self.map_symm_k = None  # k-point mapping under symmetry
        self.map_symm_g = None  # k-point mapping with time reversal
        self.map_symm_b = None  # bond mapping under symmetry
        self.map_symm_p = None  # pair mapping under symmetry

        # Translation information
        self.translback = None  # translation mapping to primitive cell
        self.translate = None  # translation action on sites

        # Symmetry operation details
        self.valid_symm = None  # valid symmetry operations
        self.symmangle = None  # rotation angles
        self.symmpoint = None  # symmetry operation points
        self.symmaxis = None  # symmetry operation axes
        self.symmlabel = None  # operation labels (C,D,I,S)

        # Status flags
        self.initialized = False
        self.lattice_mapped = False
        self.recip_lattice_mapped = False
        self.bonds_mapped = False
        self.add_time_rev = False

    @torch.no_grad()
    def read_symmetries(self, input_file: str) -> None:
        """Read symmetry operations from input file"""
        # Initialize arrays
        self.symmangle = torch.zeros(self.ntotsymm, device=self.device)
        self.symmpoint = torch.zeros((3, self.ntotsymm), device=self.device)
        self.symmaxis = torch.zeros((3, self.ntotsymm), device=self.device)
        self.symmlabel = [""] * self.ntotsymm

        # Read symmetry operations
        with open(input_file) as f:
            for i in range(self.ntotsymm):
                line = f.readline().split()
                label = line[0]
                self.symmlabel[i] = label[0]

                if label[0] in ["C", "S"]:  # Rotation/rotoreflection
                    axis_order = int(label[1:])
                    self.symmangle[i] = 2 * np.pi / axis_order
                    self.symmpoint[:, i] = torch.tensor([float(x) for x in line[1:4]], device=self.device)
                    self.symmaxis[:, i] = torch.tensor([float(x) for x in line[4:7]], device=self.device)

                elif label[0] == "D":  # Mirror plane
                    self.symmpoint[:, i] = torch.tensor([float(x) for x in line[1:4]], device=self.device)
                    self.symmaxis[:, i] = torch.tensor([float(x) for x in line[4:7]], device=self.device)

                elif label[0] == "I":  # Inversion
                    self.symmpoint[:, i] = torch.tensor([float(x) for x in line[1:4]], device=self.device)
