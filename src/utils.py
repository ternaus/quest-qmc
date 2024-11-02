import torch


def get_default_device() -> torch.device:
    """Get default device (CUDA if available, else CPU)"""
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")
