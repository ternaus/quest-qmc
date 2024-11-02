# src/dqmc.py

import torch

from config import DQMCConfig
from lattice import SquareLattice
from measurements import TimeDependentMeasurements
from workspace import Workspace


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

        # Initialize workspace for matrix operations
        self.workspace = Workspace(config.L_sites)

        # Initialize measurements
        self.measurements = TimeDependentMeasurements(
            nsites=config.L_sites,
            ntau=config.n_slices,
            nbins=config.n_bins,
            nclass=self.lattice.nclass,
        )

        # Initialize B matrices for propagation
        self.B_up = torch.zeros((self.config.L_sites, self.config.L_sites), device=self.device)
        self.B_dn = torch.zeros((self.config.L_sites, self.config.L_sites), device=self.device)
        self._init_B_matrices()

    def _init_B_matrices(self) -> None:
        """Initialize B matrices for propagation"""
        # Kinetic energy part
        self.B_up = -self.config.dt * self.lattice.hopping
        self.B_dn = self.B_up.clone()

        # Add chemical potential
        mu_mat = torch.eye(self.config.L_sites, device=self.device) * self.config.mu
        self.B_up -= mu_mat
        self.B_dn -= mu_mat

        # Exponentiate
        self.B_up = torch.matrix_exp(self.B_up)
        self.B_dn = torch.matrix_exp(self.B_dn)

    @torch.no_grad()
    def warmup(self) -> None:
        """Perform warmup sweeps"""
        for _ in range(self.config.n_warm):
            self.sweep(measure=False)

    @torch.no_grad()
    def sweep(self, measure: bool = True) -> None:
        """Perform one sweep through all HSF fields"""
        # Get workspace arrays
        U = self.workspace.R1
        D = self.workspace.R5
        T = self.workspace.R2

        # Loop over all time slices
        for slice_idx in range(self.config.n_slices):
            # Update HSF fields
            self._update_slice(slice_idx, U, D, T)

            # Compute Green's functions if needed
            if slice_idx % self.config.n_delay == 0:
                self._compute_greens(slice_idx, U, D, T)

                # Perform measurements if requested
                if measure and slice_idx % self.config.n_measure == 0:
                    self._measure(slice_idx)

    def _update_slice(self, slice_idx: int, U: torch.Tensor, D: torch.Tensor, T: torch.Tensor) -> None:
        """Update HSF fields for a given time slice"""
        # Local updates using Metropolis algorithm
        for site in range(self.config.L_sites):
            # Compute acceptance ratio
            delta_s = self._compute_delta_s(slice_idx, site)

            # Accept/reject update
            if torch.rand(1, device=self.device) < torch.min(torch.exp(-delta_s), torch.ones_like(delta_s)):
                self.hsf[slice_idx, site] *= -1
                self._update_greens(slice_idx, site, U, D, T)

    def _compute_delta_s(self, slice_idx: int, site: int) -> torch.Tensor:
        """Compute change in action for proposed update"""
        # Get current HSF value and its proposed flip
        hsf_old = self.hsf[slice_idx, site]
        hsf_new = -hsf_old

        # Compute change in interaction term
        delta_v = self.config.U * (hsf_new - hsf_old)

        # Compute determinant ratio using Green's functions
        det_ratio_up = 1 + (1 - self.G_up[site, site]) * (torch.exp(delta_v) - 1)
        det_ratio_dn = 1 + (1 - self.G_dn[site, site]) * (torch.exp(-delta_v) - 1)

        # Return change in action
        return -torch.log(det_ratio_up * det_ratio_dn)

    def _update_greens(self, slice_idx: int, site: int, U: torch.Tensor, D: torch.Tensor, T: torch.Tensor) -> None:
        """Update Green's functions after accepted update using Sherman-Morrison formula"""
        # Get interaction change
        delta_v = 2 * self.config.U * self.hsf[slice_idx, site]

        # Update up spin Green's function
        x = (torch.exp(delta_v) - 1) * (1 - self.G_up[site, site])
        self.G_up -= torch.ger(self.G_up[:, site], self.G_up[site, :]) * x / (1 + x)

        # Update down spin Green's function
        x = (torch.exp(-delta_v) - 1) * (1 - self.G_dn[site, site])
        self.G_dn -= torch.ger(self.G_dn[:, site], self.G_dn[site, :]) * x / (1 + x)

    def _compute_greens(self, slice_idx: int, U: torch.Tensor, D: torch.Tensor, T: torch.Tensor) -> None:
        """Compute Green's functions using UDV decomposition"""
        # Initialize UDV decomposition
        U.copy_(torch.eye(self.config.L_sites, device=self.device))
        D.copy_(torch.ones(self.config.L_sites, device=self.device))
        T.copy_(torch.eye(self.config.L_sites, device=self.device))

        # Multiply B matrices and interaction terms
        for i in range(slice_idx, -1, -1):
            # Multiply by B matrix
            self._multiply_B(U, D, T, up_spin=True)
            self._multiply_B(U, D, T, up_spin=False)

            # Multiply by interaction term
            self._multiply_V(i, U, D, T)

        # Compute Green's functions
        self.G_up = torch.inverse(U + T)
        self.G_dn = self.G_up.clone()

    def _multiply_B(self, U: torch.Tensor, D: torch.Tensor, T: torch.Tensor, up_spin: bool) -> None:
        """Multiply UDT decomposition by B matrix"""
        B = self.B_up if up_spin else self.B_dn
        U.copy_(torch.mm(B, U))
        T.copy_(torch.mm(T, B.t()))

    def _multiply_V(self, slice_idx: int, U: torch.Tensor, D: torch.Tensor, T: torch.Tensor) -> None:
        """Multiply UDT decomposition by interaction term"""
        V = torch.diag(torch.exp(self.config.U * self.hsf[slice_idx]))
        U.copy_(torch.mm(V, U))
        T.copy_(torch.mm(T, V))

    def _measure(self, slice_idx: int) -> None:
        """Perform measurements"""
        self.measurements.measure_tdm(
            gtau_up=self.G_up,
            gtau_dn=self.G_dn,
            sign=self._compute_sign(),
            bin_idx=self.config.current_bin,
        )

    def _compute_sign(self) -> float:
        """Compute Monte Carlo sign"""
        # For Hubbard model, sign is product of determinants
        sign_up = torch.det(torch.eye(self.config.L_sites, device=self.device) + self.G_up).sign()
        sign_dn = torch.det(torch.eye(self.config.L_sites, device=self.device) + self.G_dn).sign()
        return float(sign_up * sign_dn)
