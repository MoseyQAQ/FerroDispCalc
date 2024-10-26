from abc import ABC, abstractmethod
from ase import Atoms
import numpy as np

class ComputeBackend(ABC):
    def __init__(self, input, type_map: list[str]=None, prefix: str=None):
        self.traj = input
        self.type_map = type_map
        self.prefix = prefix
    
    @abstractmethod
    def get_averaged_structure(self, select: list[int]):
        pass 

    @abstractmethod
    def get_displacement(self, select: list[int], nl: np.ndarray):
        pass

    @abstractmethod
    def get_polarization(self, select: list[int], nl: np.ndarray):
        pass