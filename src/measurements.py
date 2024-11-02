# src/measurements.py
from dataclasses import dataclass
from enum import IntEnum

import torch

DEVICE = torch.get_default_device()


class Observable(IntEnum):
    UP_OCC = 0  # Up spin occupancy
    DN_OCC = 1  # Down spin occupancy
    POT_ENERGY = 2  # Potential energy
    KIN_ENERGY = 3  # Kinetic energy
    TOT_ENERGY = 4  # Total energy
    DENSITY = 5  # Density
    XX_FM_SF = 6  # XX ferromagnetic structure factor
    ZZ_FM_SF = 7  # ZZ ferromagnetic structure factor
    XX_AF_SF = 8  # XX antiferromagnetic structure factor
    ZZ_AF_SF = 9  # ZZ antiferromagnetic structure factor
    RMS_XX_AF = 10  # RMS XX AF SF
    RMS_ZZ_AF = 11  # RMS ZZ AF SF
    AVG_SIGN = 12  # Average sign
    UP_SIGN = 13  # Average up sign
    DN_SIGN = 14  # Average down sign


@dataclass
class PhysicalMeasurements:
    """Physical measurements during DQMC simulation"""

    n: int  # Number of sites
    n_bins: int  # Number of measurement bins
    n_class: int  # Number of correlation classes
    device: torch.device = DEVICE

    def __post_init__(self):
        # Initialize scalar measurements
        self.scalars = torch.zeros((len(Observable), self.n_bins), device=self.device)

        # Initialize correlation functions
        self.green_fn = torch.zeros((self.n_class, self.n_bins), device=self.device)
        self.dens_uu = torch.zeros((self.n_class, self.n_bins), device=self.device)
        self.dens_ud = torch.zeros((self.n_class, self.n_bins), device=self.device)
        self.spin_xx = torch.zeros((self.n_class, self.n_bins), device=self.device)
        self.spin_zz = torch.zeros((self.n_class, self.n_bins), device=self.device)

        # Statistics
        self.curr_bin = 0
        self.measurements = 0
        self.avg_bin = self.n_bins
        self.err_bin = self.n_bins + 1

    @torch.no_grad()
    def measure(
        self,
        G_up: torch.Tensor,
        G_dn: torch.Tensor,
        mu: tuple[float, float],
        t: torch.Tensor,
        sign_up: float,
        sign_dn: float,
    ) -> None:
        """Perform measurements using Green's functions"""
        # Occupancies
        n_up = torch.diagonal(G_up).sum() / self.n
        n_dn = torch.diagonal(G_dn).sum() / self.n

        # Energies
        pot = self.compute_potential(G_up, G_dn)
        kin = self.compute_kinetic(G_up, G_dn, t)

        # Store scalar measurements
        self.scalars[Observable.UP_OCC, self.curr_bin] += n_up
        self.scalars[Observable.DN_OCC, self.curr_bin] += n_dn
        self.scalars[Observable.POT_ENERGY, self.curr_bin] += pot
        self.scalars[Observable.KIN_ENERGY, self.curr_bin] += kin
        self.scalars[Observable.TOT_ENERGY, self.curr_bin] += pot + kin
        self.scalars[Observable.DENSITY, self.curr_bin] += (n_up + n_dn) / 2

        # Signs
        self.scalars[Observable.UP_SIGN, self.curr_bin] += sign_up
        self.scalars[Observable.DN_SIGN, self.curr_bin] += sign_dn
        self.scalars[Observable.AVG_SIGN, self.curr_bin] += sign_up * sign_dn

        # Correlation functions
        self.measure_correlations(G_up, G_dn)

        self.measurements += 1
