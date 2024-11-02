# src/lattice.py
from dataclasses import dataclass

import torch

DEVICE = torch.get_default_device()


@dataclass
class Lattice:
    """Real space lattice implementation"""

    nsites: int  # number of total sites
    natom: int  # number of sites in primitive cell
    ncell: int  # number of cells in supercell
    ndim: int  # number of extended dimensions
    device: torch.device = DEVICE

    def __post_init__(self):
        # Initialize lattice vectors and positions
        self.sc = torch.zeros((3, 3), device=self.device)  # supercell in fractional coords
        self.ac = torch.zeros((3, 3), device=self.device)  # primitive cell in cartesian
        self.scc = torch.zeros((3, 3), device=self.device)  # supercell in cartesian

        # Site positions
        self.pos = torch.zeros((3, self.nsites), device=self.device)  # fractional
        self.cartpos = torch.zeros((3, self.nsites), device=self.device)  # cartesian
        self.xat = torch.zeros((3, self.natom), device=self.device)  # primitive cell

        # Phases and translations
        self.phase = torch.zeros(self.nsites, device=self.device)
        self.translation = torch.zeros((3, self.ncell), device=self.device)

        # Classification
        self.nclass = 0
        self.myclass = None
        self.class_size = None
        self.class_label = None
        self.gf_phase = None

        # Labels
        self.olabel = [""] * self.natom

        # Status flags
        self.initialized = False
        self.constructed = False
        self.analyzed = False

    def assign_gf_phase(self, twist: torch.Tensor) -> None:
        """Assign phases for Green's function calculation"""
        self.gf_phase = torch.zeros((self.nsites, self.nsites), device=self.device)

        for ic in range(self.nclass):
            csize = 0
            for j in range(self.nsites):
                for i in range(self.nsites):
                    if self.myclass[i, j] != ic:
                        continue

                    d = self.translation[:, j // self.natom] - self.translation[:, i // self.natom]

                    if csize == 0:
                        d0 = d.clone()
                        self.gf_phase[j, i] = 1
                    else:
                        # Calculate phases
                        rphase = torch.cos(torch.sum(twist * (d0 - d)))
                        iphase = torch.sin(torch.sum(twist * (d0 - d)))

                        if torch.abs(iphase) < 1e-6:
                            rounded = torch.round(rphase)
                            if torch.abs(rphase - rounded) < 1e-6:
                                self.gf_phase[j, i] = rounded
                            else:
                                # Try alternative phase
                                rphase = torch.cos(torch.sum(twist * (d0 + d)))
                                iphase = torch.sin(torch.sum(twist * (d0 + d)))
                                if torch.abs(iphase) < 1e-6:
                                    rounded = torch.round(rphase)
                                    if torch.abs(rphase - rounded) < 1e-6:
                                        self.gf_phase[j, i] = rounded
                                    else:
                                        raise ValueError("Problem with phase")

                    csize += 1

            if csize != self.class_size[ic]:
                raise ValueError("Problem with classes")
