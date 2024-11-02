# src/workspace.py
from dataclasses import dataclass
from enum import IntEnum

import torch

DEVICE = torch.get_default_device()


class LapackOp(IntEnum):
    SYEV = 0  # Symmetric eigenvalue decomposition
    GEQRF = 1  # QR factorization
    ORGQR = 2  # Generate Q matrix
    ORMQR = 3  # Multiply by Q matrix
    GETRI = 4  # Matrix inversion


@dataclass
class Workspace:
    """Memory workspace for DQMC calculations"""

    n: int  # matrix dimension
    device: torch.device = DEVICE

    def __post_init__(self):
        # Real matrices
        self.R1 = torch.zeros((self.n, self.n), device=self.device)
        self.R2 = torch.zeros((self.n, self.n), device=self.device)
        self.R3 = torch.zeros((self.n, self.n), device=self.device)
        self.R4 = torch.zeros((self.n, self.n), device=self.device)
        self.R8 = torch.zeros((self.n, self.n), device=self.device)

        # Real vectors
        self.R5 = torch.zeros(self.n, device=self.device)
        self.R6 = torch.zeros(self.n, device=self.device)

        # Integer arrays
        self.I1 = torch.zeros(self.n, dtype=torch.int64, device=self.device)
        self.I2 = torch.zeros(self.n, dtype=torch.int64, device=self.device)

        # LAPACK workspace sizes
        self.lwork = torch.zeros(len(LapackOp), dtype=torch.int64, device=self.device)
        self._compute_workspace_sizes()

        # Allocate workspace
        max_work = self.lwork.max().item()
        self.R7 = torch.zeros(max_work, device=self.device)

    def _compute_workspace_sizes(self) -> None:
        """Compute optimal workspace sizes for LAPACK operations"""
        n = self.n

        # SYEV workspace
        self.lwork[LapackOp.SYEV] = max(1, 3 * n - 1)

        # GEQRF workspace
        self.lwork[LapackOp.GEQRF] = max(1, n)

        # ORGQR workspace
        self.lwork[LapackOp.ORGQR] = max(1, n)

        # ORMQR workspace
        self.lwork[LapackOp.ORMQR] = max(1, n)

        # GETRI workspace
        self.lwork[LapackOp.GETRI] = max(1, n)
