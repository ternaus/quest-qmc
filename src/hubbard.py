# src/hubbard.py
from dataclasses import dataclass
from typing import tuple

import torch

from .green import GreenFunction

DEVICE = torch.get_default_device()


@dataclass
class HubbardModel:
    """Hubbard model implementation"""

    n: int  # Number of sites
    L: int  # Number of time slices
    beta: float  # Inverse temperature
    U: torch.Tensor  # Hubbard U
    mu: tuple[float, float]  # Chemical potential (up, down)
    t: tuple[torch.Tensor, torch.Tensor]  # Hopping (up, down)
    device: torch.device = DEVICE

    def __post_init__(self):
        # Initialize HSF fields
        self.HSF = torch.randint(2, (self.L, self.n), device=self.device) * 2 - 1
        self.lambda_vec = torch.sqrt(torch.abs(self.U) * self.beta / self.L)

        # Initialize V matrices
        self.V_up = torch.zeros((self.n, self.L), device=self.device)
        self.V_dn = torch.zeros((self.n, self.L), device=self.device)
        self.update_V()

    @torch.no_grad()
    def update_V(self):
        """Update V matrices based on HSF configuration"""
        if torch.all(self.U >= 0):  # Attractive U
            temp = self.HSF * self.lambda_vec.view(-1, 1)
            self.V_up = torch.exp(temp)
            self.V_dn = torch.exp(-temp)
        else:  # Repulsive U
            temp = self.HSF * self.lambda_vec.view(-1, 1)
            self.V_up = torch.exp(temp)
            self.V_dn = self.V_up  # Point to same tensor

    @torch.no_grad()
    def compute_greens(self) -> tuple[GreenFunction, GreenFunction]:
        """Compute Green's functions for up and down spins"""
        G_up = GreenFunction(
            n=self.n,
            L=self.L,
            G=torch.zeros((self.n, self.n), device=self.device),
            V=self.V_up,
            device=self.device,
        )

        G_dn = GreenFunction(
            n=self.n,
            L=self.L,
            G=torch.zeros((self.n, self.n), device=self.device),
            V=self.V_dn,
            device=self.device,
        )

        # Initialize Green's functions
        G_up.compute_G()
        G_dn.compute_G()

        return G_up, G_dn
