# src/structure.py
from dataclasses import dataclass
from enum import IntEnum
from typing import list

import torch

DEVICE = torch.get_default_device()


class StructureFlag(IntEnum):
    INIT = 0
    DIM = 1
    ADJ = 2
    CLASS = 3
    WAVE = 4
    BOND = 5
    PHASE = 6
    FT = 7


@dataclass
class CompressedColumnStorage:
    """Compressed Column Storage for sparse matrices"""

    n: int  # dimension of matrix
    nnz: int  # number of nonzeros
    values: torch.Tensor  # nonzero elements
    row_indices: torch.Tensor  # row indices
    col_ptrs: torch.Tensor  # column pointers
    device: torch.device = DEVICE

    def to_dense(self) -> torch.Tensor:
        """Convert to dense tensor"""
        dense = torch.zeros((self.n, self.n), device=self.device)
        for j in range(self.n):
            start, end = self.col_ptrs[j], self.col_ptrs[j + 1]
            dense[self.row_indices[start:end], j] = self.values[start:end]
        return dense


@dataclass
class LatticeStructure:
    """Lattice structure for DQMC"""

    nsites: int  # number of sites
    ncell: int  # number of unit cells
    name: str  # structure name
    device: torch.device = DEVICE

    def __post_init__(self):
        # Initialize structure components
        self.dim: torch.Tensor | None = None
        self.n_t: int = 0  # number of hopping types
        self.T = None  # hopping matrix (CCS format)

        # Distance classification
        self.nclass: int = 0
        self.D = torch.zeros((self.nsites, self.nsites), dtype=torch.int64, device=self.device)
        self.F = None  # frequency of distances
        self.clabel: list[str] = []

        # Phase and wave information
        self.P = torch.zeros((self.nsites, self.nsites), device=self.device)
        self.nwave: int = 0
        self.W = None
        self.wlabel: list[str] = []

        # Fourier transform
        self.FT = None

        # Status flags
        self.checklist = torch.zeros(len(StructureFlag), dtype=torch.bool, device=self.device)
