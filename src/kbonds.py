# src/kbonds.py
from dataclasses import dataclass

import torch

from .lattice import Lattice

DEVICE = torch.get_default_device()


@dataclass
class KBonds:
    """K-space bonds and symmetries"""

    nak: int  # number of (atom,k) pairs
    nbonds: int  # number of bonds
    nmomenta: int  # number of momenta
    device: torch.device = DEVICE

    def __post_init__(self):
        # Initialize symmetry maps
        self.map_symm_ak = None  # maps (atom,k) under symmetry
        self.map_symm_bak = None  # maps bonds under symmetry
        self.class_size = None  # size of each class
        self.myclass = None  # class assignments
        self.nclass = None  # number of classes

        # Bond information
        self.bond_origin = None  # origin sites of bonds
        self.bond_target = None  # target sites of bonds
        self.bmap = None  # bond mappings
        self.ksum = None  # k-space sums

    def construct_kbond_classes(self, lattice: Lattice) -> None:
        """Construct classes of k-space bonds"""
        ntotbond = self.nbonds
        maxclass = 0

        # Initialize arrays
        self.myclass = torch.zeros((ntotbond, ntotbond, self.nmomenta), dtype=torch.int64, device=self.device)

        for im in range(self.nmomenta):
            # Process each momentum
            bond1 = []  # Will store bond origins
            bond2 = []  # Will store bond targets
            csizev = []  # Will store class sizes

            mclass = 0
            for i in range(ntotbond):
                for j in range(ntotbond):
                    # Process bond pairs
                    found = False
                    for ic in range(len(csizev)):
                        if self._check_bond_equivalence(i, j, bond1[ic], bond2[ic]):
                            found = True
                            break

                    if not found:
                        # New class
                        bond1.append(i)
                        bond2.append(j)
                        csizev.append(1)
                        mclass += 1

            # Store class information
            self.nclass[im] = mclass
            maxclass = max(maxclass, mclass)

        # Allocate and fill class sizes
        self.class_size = torch.zeros((maxclass, self.nmomenta), device=self.device)
        for im in range(self.nmomenta):
            for ib in range(ntotbond):
                for jb in range(ntotbond):
                    class_idx = self.myclass[ib, jb, im]
                    self.class_size[class_idx, im] += 1
