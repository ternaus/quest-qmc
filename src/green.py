# src/green.py
from dataclasses import dataclass

import torch

DEVICE = torch.get_default_device()


@dataclass
class GreenFunction:
    n: int  # Number of sites
    L: int  # Number of time slices
    G: torch.Tensor  # Green's function matrix
    V: torch.Tensor  # HSF field matrices
    ilb: int  # Index of leftmost B
    det: float  # Determinant
    sign: float  # Sign of det(G)

    # Numerical stability parameters
    n_wrap: int = 5  # Wrapping frequency
    wps: int = 0  # Wraps performed
    last_wrap: int = 0  # Last wrap position
    max_wrap: int = 100  # Maximum wraps
    fix_wrap: int = -1  # Fixed wrap point
    diff_lim: float = 0.1  # Difference limit
    err_rate: float = 0.0  # Error rate

    # Block update parameters
    n_blk: int = 10  # Block size for updates
    blk_sz: int = 0  # Current block size
    U: torch.Tensor | None = None  # Update matrices
    W: torch.Tensor | None = None
    n_modify: int = 0

    device: torch.device = DEVICE

    def __post_init__(self):
        if self.U is None:
            self.U = torch.zeros((self.n, self.n_blk), device=self.device)
        if self.W is None:
            self.W = torch.zeros((self.n, self.n_blk), device=self.device)

    @torch.no_grad()
    def update_G(self, j: int, gamma: float) -> None:
        """Update Green's function with rank-1 update"""
        blk_sz = self.blk_sz

        # Compute x and y vectors
        x = self.G[:, j].clone()
        y = self.G[j, :].clone()
        y[j] -= 1.0

        if blk_sz > 0:
            # Add effects of previous updates
            xx = self.W[j, :blk_sz]
            x += torch.mv(self.U[:, :blk_sz], xx)
            xx = self.U[j, :blk_sz]
            y += torch.mv(self.W[:, :blk_sz], xx)

        y *= gamma

        # Store update vectors
        self.U[:, blk_sz] = x
        self.W[:, blk_sz] = y
        self.blk_sz += 1

        # Apply updates if block is full
        self.apply_update()

    @torch.no_grad()
    def apply_update(self, forced: bool = False) -> None:
        """Apply blocked updates to G"""
        if forced or self.blk_sz == self.n_blk:
            blk_sz = self.blk_sz
            self.G += torch.mm(self.U[:, :blk_sz], self.W[:, :blk_sz].t())
            self.blk_sz = 0
