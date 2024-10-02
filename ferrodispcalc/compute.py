'''
This a wrapper calss for the computation of the ferroelectric displacement and polarization.
It has two backends: c++ or python.
When the input is a list of pymatgen.Structure or ase.Atoms, the python backend is used.
When the input is a string, the c++ backend is used. (lmp-dump file), for better performance.

The class has the following methods:
- get_avgeraged_structure: 
    Get the averaged structure of the input structures.
- polarization:
    Compute the polarization of the input structures. (only support ABO3 perovskite)
- displacement:
    Compute the displacement of the input structures.
- lattice:
    Compute the lattice parameters of the input structures. (only support cubic perovskite)
- rotation:
    Compute the rotation of the input structures. (only support cubic perovskite)
'''

from ase import Atoms
from pymatgen.core import Structure
import numpy as np
from fdc import get_averaged_structure
class Compute:
    def __init__(self, input: str | list[Structure] | list[Atoms], 
                 format: str=None,
                 type_map: list[str]=None) -> None:
        
        self.input = input
        self.format = format
        self.type_map = type_map
        
        # determine the backend
        if format == 'lmp-dump':
            self.backend = 'cpp'
        else:
            self.backend = 'py'
        
        # determine the type of the input
        if isinstance(input, list):
            if isinstance(input[0], Structure):
                self.input_type = 'pymatgen'
            elif isinstance(input[0], Atoms):
                self.input_type = 'ase'
        elif isinstance(input, str):
            self.input_type = 'lmp-dump'
        else:
            raise ValueError('Invalid input type')

    def get_avgeraged_structure(self, select: slice=None) -> Structure:
        pass 

    def get_polarization(self, select: slice, nl_ba: np.ndarray, nl_bo: np.ndarray, born_effective_charge: dict):
        pass
    
    def get_displacement(self, select: slice, nl: np.ndarray):
        pass 
    
    def get_local_lattice(self, select: slice, nl_ba: np.ndarray, rotate:dict):
        pass 
    
    def get_octahedron_rotation(self, select: slice, nl_bo: np.ndarray, rotate:dict):
        pass 

# ------------------- Python backend ------------------- #
def calculate_avgeraged_structure(input: list[Structure] | list[Atoms]) -> Structure:
    pass

def calculate_polarization(input: list[Structure] | list[Atoms]) -> np.ndarray:
    pass

def calculate_displacement(input: list[Structure] | list[Atoms]) -> np.ndarray:
    pass

def calculate_local_lattice(input: list[Structure] | list[Atoms]) -> np.ndarray:
    pass

def calculate_octahedron_rotation(input: list[Structure] | list[Atoms]) -> np.ndarray:
    pass

# ------------------- C++ backend ------------------- #
pass 