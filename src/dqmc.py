# src/dqmc.py
import torch

from config import DQMCConfig
from lattice import SquareLattice


class DQMC:
    def __init__(self, config: DQMCConfig):
        self.config = config
        self.device = config.device
        self.lattice = SquareLattice(config)

        # Initialize HSF fields
        self.hsf = torch.randint(2, (self.config.n_slices, self.config.L_sites), device=self.device) * 2 - 1

        # Initialize Green's functions
        self.G_up = torch.zeros((self.config.L_sites, self.config.L_sites), device=self.device)
        self.G_dn = torch.zeros((self.config.L_sites, self.config.L_sites), device=self.device)

    @torch.no_grad()
    def warmup(self):
        """Perform warmup sweeps"""
        for _ in range(self.config.n_warm):
            self.sweep(measure=False)

    @torch.no_grad()
    def sweep(self, measure: bool = True):
        """Perform one sweep through all HSF fields"""
        # Implementation coming soon
