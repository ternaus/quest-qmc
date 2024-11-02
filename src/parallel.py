# src/parallel.py
from dataclasses import dataclass
from enum import IntEnum

import torch
import torch.distributed as dist

DEVICE = torch.get_default_device()


class Channel(IntEnum):
    AGGR = 1
    MEAS = 2
    GFUN = 3


class ParallelLevel(IntEnum):
    LEVEL_1 = 1
    LEVEL_2 = 2
    LEVEL_3 = 3
    LEVEL_4 = 4


@dataclass
class ParallelQMC:
    """Parallel processing configuration"""

    level: int  # parallelization level
    rank: int = 0  # global rank
    size: int = 1  # global size
    device: torch.device = DEVICE

    def __post_init__(self):
        # Initialize ranks
        self.aggr_rank = self.rank
        self.meas_rank = self.rank
        self.gfun_rank = self.rank

        # Initialize sizes
        self.aggr_size = self.size
        self.meas_size = self.size
        self.gfun_size = self.size

        # Initialize roots
        self.aggr_root = 0
        self.meas_root = 0
        self.gfun_root = 0

        # Initialize communicators
        if dist.is_initialized():
            self.aggr_comm = dist.group.WORLD
            self.meas_comm = dist.group.WORLD
            self.gfun_comm = dist.group.WORLD

    @staticmethod
    def init(level: int = ParallelLevel.LEVEL_1) -> "ParallelQMC":
        """Initialize parallel processing"""
        if not dist.is_initialized():
            dist.init_process_group("nccl" if torch.cuda.is_available() else "gloo")

        return ParallelQMC(level=level, rank=dist.get_rank(), size=dist.get_world_size())

    def is_root(self, channel: Channel) -> bool:
        """Check if current process is root for given channel"""
        if channel == Channel.AGGR:
            return self.aggr_rank == self.aggr_root
        return False
