# src/tdm.py
from dataclasses import dataclass

import torch

DEVICE = torch.get_default_device()


@dataclass
class TimeDependentMeasurements:
    """Time-dependent measurements for DQMC"""

    nsites: int  # number of sites
    ntau: int  # number of time slices
    nbins: int  # number of measurement bins
    nclass: int  # number of correlation classes
    device: torch.device = DEVICE

    def __post_init__(self):
        # Initialize measurement arrays
        self.gfun = torch.zeros((self.nclass, self.ntau + 1, self.nbins + 2), device=self.device)
        self.gfun_up = torch.zeros_like(self.gfun)
        self.gfun_dn = torch.zeros_like(self.gfun)
        self.spxx = torch.zeros_like(self.gfun)
        self.spzz = torch.zeros_like(self.gfun)
        self.dens = torch.zeros_like(self.gfun)
        self.pair = torch.zeros_like(self.gfun)
        self.cond = torch.zeros_like(self.gfun)

        # Working arrays
        self.work = torch.zeros((self.nclass, self.ntau + 1), device=self.device)
        self.temp = torch.zeros((self.ntau + 1, 2), device=self.device)

        # Labels
        self.labels = [""] * self.nclass

    @torch.no_grad()
    def measure_tdm(self, gtau_up: torch.Tensor, gtau_dn: torch.Tensor, sign: float, bin_idx: int) -> None:
        """Perform time-dependent measurements"""
        # Single-particle Green's functions
        self.gfun_up[..., bin_idx] += sign * gtau_up
        self.gfun_dn[..., bin_idx] += sign * gtau_dn
        self.gfun[..., bin_idx] += sign * (gtau_up + gtau_dn) / 2

        # Spin correlations
        self.spxx[..., bin_idx] += sign * (gtau_up * gtau_dn)
        self.spzz[..., bin_idx] += sign * (gtau_up - gtau_dn) * (gtau_up - gtau_dn) / 4

        # Density correlations
        self.dens[..., bin_idx] += sign * (gtau_up + gtau_dn) * (gtau_up + gtau_dn) / 4
