# src/geometry.py
from dataclasses import dataclass
from pathlib import Path

import torch

DEVICE = torch.get_default_device()


@dataclass
class LatticeConfig:
    ndim: int
    natom: int
    nsites: int
    primitive_vectors: torch.Tensor  # (ndim, ndim)
    atom_positions: torch.Tensor  # (natom, ndim)
    device: torch.device = DEVICE


@dataclass
class ReciprocalLattice:
    nkpts: int
    nclass_k: int
    klist: torch.Tensor  # (nkpts, ndim)
    myclass_k: torch.Tensor  # (nkpts)
    include_gamma: bool = False
    device: torch.device = DEVICE


@dataclass
class Hamiltonian:
    t_up: torch.Tensor  # Hopping parameters up spin
    t_dn: torch.Tensor  # Hopping parameters down spin
    U: torch.Tensor  # Hubbard U values
    mu_up: torch.Tensor  # Chemical potential up
    mu_dn: torch.Tensor  # Chemical potential down
    device: torch.device = DEVICE


class GeometryWrapper:
    def __init__(self, config: dict):
        self.config = config
        self.device = config.get("device", torch.device("cuda" if torch.cuda.is_available() else "cpu"))

        # Initialize components
        self.lattice: LatticeConfig | None = None
        self.recip_lattice: ReciprocalLattice | None = None
        self.gamma_lattice: ReciprocalLattice | None = None
        self.hamiltonian: Hamiltonian | None = None

    def init_from_file(self, geom_file: Path) -> None:
        """Initialize geometry from input file"""
        if not geom_file.exists():
            raise FileNotFoundError(f"Cannot open geometry definition file {geom_file}")

        # Parse input file
        with open(geom_file) as f:
            self.input_content = f.read()

        # Initialize lattice
        self._init_lattice()

        # Construct full lattice
        self._construct_lattice()

        # Initialize and construct reciprocal lattices
        self._init_recip_lattice(include_gamma=False)  # Regular k-points
        self._init_recip_lattice(include_gamma=True)  # Including Gamma point

        # Construct Hamiltonian
        self._construct_hamiltonian()

    def _init_lattice(self) -> None:
        """Initialize basic lattice information"""
        # Implementation will depend on input file format

    def _construct_lattice(self) -> None:
        """Construct full lattice including neighbors"""

    def _init_recip_lattice(self, include_gamma: bool = False) -> None:
        """Initialize reciprocal lattice"""

    def print_header_ft(self, apply_twist: bool = True) -> None:
        """Print Fourier transform grid information"""
        lattice = self.recip_lattice if apply_twist else self.gamma_lattice

        print(" Grid for " + ("Green's function" if apply_twist else "spin/charge correlations"))
        print("  K-points")
        print("  Class")

        for i in range(lattice.nclass_k):
            first = True
            for k in range(lattice.nkpts):
                if lattice.myclass_k[k] == i:
                    if first:
                        print(f"  {i+1:3d}      {' '.join(f'{x:10.5f}' for x in lattice.klist[k])}")
                        first = False
                    else:
                        print(f"           {' '.join(f'{x:10.5f}' for x in lattice.klist[k])}")
            print()
