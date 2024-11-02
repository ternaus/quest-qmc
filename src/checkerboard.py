# src/checkerboard.py
from dataclasses import dataclass

import torch

DEVICE = torch.get_default_device()


@dataclass
class MatB:
    """Checkerboard decomposition matrix B"""

    n: int  # dim of B
    m: int  # number of neighbors of lattice
    A: torch.Tensor  # Adjacency info (3, m)
    sinht: torch.Tensor  # sinh(t)
    cosht: torch.Tensor  # cosh(t)
    exptaumu: torch.Tensor  # parameters for checkerboard method
    work: torch.Tensor  # work array
    name: str = "Checkerboard"
    device: torch.device = DEVICE


class CheckerBoard:
    def __init__(self, config: dict):
        self.config = config
        self.device = config.get("device", torch.device("cuda" if torch.cuda.is_available() else "cpu"))

    def init_B(
        self,
        n: int,
        adj: torch.Tensor,
        ckb: torch.Tensor,
        t: torch.Tensor,
        mu: torch.Tensor,
        dtau: float,
    ) -> MatB:
        """Initialize MatB structure"""
        if not torch.all(ckb.row == adj.row) or not torch.all(ckb.cstart == adj.cstart):
            raise ValueError("ckb and adj do not conform")

        m = ckb.nnz // 2
        nckb = ckb.max()
        nt = adj.max()

        # Initialize parameters
        sinht = torch.sinh(dtau * 0.5 * t[:nt])
        cosht = torch.cosh(dtau * 0.5 * t[:nt])
        exptaumu = torch.exp(dtau * mu)

        # Build adjacency info
        A = self._build_adjacency(n, ckb, adj, nckb)

        return MatB(
            n=n,
            m=m,
            A=A.to(self.device),
            sinht=sinht.to(self.device),
            cosht=cosht.to(self.device),
            exptaumu=exptaumu.to(self.device),
            work=torch.zeros(n, device=self.device),
        )
