# src/kernels.py
import torch


@torch.jit.script
def diag(A: torch.Tensor) -> torch.Tensor:
    """Extract diagonal of matrix A"""
    return torch.diagonal(A, 0)


@torch.jit.script
def normcol(A: torch.Tensor, D: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    """Normalize columns of A by D and compute norms"""
    # Scale columns
    A = A * D.unsqueeze(-2)  # Broadcasting for column scaling
    # Compute column norms
    c = torch.sum(A * A, dim=-2)
    return A, c


@torch.jit.script
def permute(A: torch.Tensor, ipiv: torch.Tensor) -> torch.Tensor:
    """Permute rows of A according to pivot indices"""
    return A[ipiv]


@torch.jit.script
def scalerowperm(D: torch.Tensor, Q: torch.Tensor, ipiv: torch.Tensor) -> torch.Tensor:
    """Scale rows and permute: T = D^-1 * R * P"""
    n = Q.shape[0]
    T = torch.zeros_like(Q)

    # Create mask for upper triangular part
    mask = torch.triu(torch.ones(n, n, dtype=torch.bool, device=Q.device))

    # Scale and permute in one operation
    T[ipiv] = torch.where(mask, Q / D.unsqueeze(-2), torch.zeros_like(Q))
    return T


@torch.jit.script
def scalerow(h: torch.Tensor, B: torch.Tensor) -> torch.Tensor:
    """Scale rows of B by h"""
    return B * h.unsqueeze(-2)  # Broadcasting for row scaling


@torch.jit.script
def scalerowcol(h: torch.Tensor, G: torch.Tensor) -> torch.Tensor:
    """Scale rows and columns of G by h"""
    return (G * h.unsqueeze(-2)) / h.unsqueeze(-1)  # Broadcasting for row and column scaling


@torch.jit.script
def scalerowadd(
    Db: torch.Tensor,
    U: torch.Tensor,
    D: torch.Tensor,
    T: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Scale rows and add: G = U/Db, T = G + D*T"""
    G = U.t() / Db.unsqueeze(-1)  # Transpose and scale
    T = G + D.unsqueeze(-1) * T
    return G, T


@torch.jit.script
def sort_pivot(Db: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    """Sort Db and return permutation indices"""
    values, indices = torch.sort(Db, descending=True)
    return values, indices
