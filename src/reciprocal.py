# src/reciprocal.py
from dataclasses import dataclass

import numpy as np
import torch

from .lattice import Lattice

DEVICE = torch.get_default_device()


@dataclass
class ReciprocalLattice:
    """Reciprocal lattice implementation"""

    ndim: int
    nkpts: int  # number of k-points (equal to ncell)
    device: torch.device = DEVICE

    def __post_init__(self):
        # Initialize k-points and vectors
        self.klist = torch.zeros((self.nkpts, self.ndim), device=self.device)
        self.kpoint = torch.zeros(self.ndim, device=self.device)
        self.ktwist = torch.zeros(self.ndim, device=self.device)

        # Reciprocal lattice vectors
        self.kcs = torch.zeros((self.ndim, self.ndim), device=self.device)
        self.ks = torch.zeros((self.ndim, self.ndim), device=self.device)
        self.kc = torch.zeros((self.ndim, self.ndim), device=self.device)

        # K-point classification
        self.nclass_k = 0
        self.myclass_k = None
        self.class_size_k = None
        self.class_repr_k = None

        # Fourier transform
        self.fourier_c = None

        # Status flags
        self.initialized = False
        self.constructed = False
        self.analyzed = False

    def init_from_lattice(self, lattice: Lattice, apply_twist: bool = True) -> None:
        """Initialize reciprocal lattice from real space lattice"""
        if not lattice.initialized:
            raise ValueError("Lattice must be initialized first")

        # Get reciprocal lattice vectors
        self.kc = 2 * np.pi * torch.linalg.inv(lattice.ac.t())
        self.kcs = 2 * np.pi * torch.linalg.inv(lattice.scc.t())

        # Apply twist if requested
        if apply_twist:
            self.ktwist = self.kpoint * np.pi

        self.initialized = True
