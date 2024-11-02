# src/gtau.py
from dataclasses import dataclass

import torch

DEVICE = torch.get_default_device()


@dataclass
class GTau:
    """Time-dependent Green's function"""

    n: int  # number of sites
    L: int  # number of time slices
    nb: int  # number of time slices in A_up, A_dn
    nnb: int  # size of A_up, A_dn
    it_up: int  # time indices
    i0_up: int
    it_dn: int
    i0_dn: int
    north: int  # Number of orthogonalizations
    device: torch.device = DEVICE

    def __post_init__(self):
        # Initialize matrices on GPU
        self.U0up = torch.zeros((self.n, self.n), device=self.device)
        self.U0dn = torch.zeros((self.n, self.n), device=self.device)
        self.e0up = torch.zeros(self.n, device=self.device)
        self.e0dn = torch.zeros(self.n, device=self.device)
        self.g0_stored = False

    @torch.no_grad()
    def get_g0(self, tau_idx: int, spin: int = 0) -> torch.Tensor:
        """Get non-interacting Green's function

        Args:
            tau_idx: Time slice index
            spin: Spin index (0=both, 1=up, -1=down)
        """
        if not self.g0_stored:
            # Compute eigendecomposition of B matrices
            self.e0up, self.U0up = torch.linalg.eigh(self.B_up)
            self.e0dn, self.U0dn = torch.linalg.eigh(self.B_dn)
            self.g0_stored = True

        g0tau = torch.zeros((self.n, self.n), device=self.device)

        if spin in [0, 1]:  # UP or BOTH
            g = self.e0up.pow(tau_idx) / (1.0 + self.e0up.pow(self.L))
            W2 = self.U0up.clone()
            W2 = W2 * g.view(1, -1)
            g0tau += torch.mm(W2, self.U0up.t())

        if spin in [0, -1]:  # DOWN or BOTH
            g = self.e0dn.pow(tau_idx) / (1.0 + self.e0dn.pow(self.L))
            W2 = self.U0dn.clone()
            W2 = W2 * g.view(1, -1)
            g0tau += torch.mm(W2, self.U0dn.t())

        if spin == 0:  # BOTH
            g0tau *= 0.5

        return g0tau
