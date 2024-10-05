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
from pymatgen.core import Structure, Lattice
from ferrodispcalc.io import LAMMPSdump
import numpy as np
import os
from fdc import get_averaged_structure, get_displacement, get_polarization

# ------------------- Main class ------------------- #
class Compute:
    def __init__(self, input: str | list[Structure] | list[Atoms], 
                 format: str=None,
                 type_map: list[str]=None,
                 prefix: str=None) -> None:
        
        self.input = input
        self.format = format
        self.type_map = type_map
        self.prefix = prefix
        
        self.backend = self.__set_backend()
        self.input_type = self.__checkinput()

    def get_averaged_structure(self, select: slice=None) -> Structure:
        
        # Workflow:
        # 1. convert slice to list
        # 2. perform the computation, depending on the backend
        # 3. store the result in self.stru, and return it.
        select = convert_slice_to_list(self.input, select)
        if self.backend == 'py':
            pass 
        elif self.backend == 'cpp':
            stru = calculate_averaged_structure_cpp(self.input, self.type_map, select)
        else:
            raise ValueError('Invalid backend') 
        
        self.stru = stru
        return stru

    def get_polarization(self, select: slice, nl_ba: np.ndarray, nl_bx: np.ndarray, born_effective_charge: dict) -> np.ndarray:
        select = convert_slice_to_list(self.input, select)
        if self.backend == 'py':
            pass
        elif self.backend == 'cpp':
            polar = calculate_polarization_cpp(self.input, select, nl_ba, nl_bx, born_effective_charge, self.type_map)
        self.polar = polar
        return polar
    
    def get_displacement(self, nl: np.ndarray | str, select: slice=None) -> np.ndarray:
        select = convert_slice_to_list(self.input, select)
        if self.backend == 'py':
            pass
        elif self.backend == 'cpp':
            disp = calculate_displacement_cpp(self.input, select, nl)
        self.disp = disp
        return disp
    
    def get_local_lattice(self, select: slice, nl_ba: np.ndarray, rotate:dict):
        raise NotImplementedError('This method is not implemented yet')
    
    def get_octahedron_rotation(self, select: slice, nl_bx: np.ndarray, rotate:dict):
        raise NotImplementedError('This method is not implemented yet')

    def __set_backend(self) -> str:
        '''
        set the backend of the computation.
        The backend can be 'cpp' or 'py'.

        Returns:
        -------
        str:
            The backend of the computation.
        '''
        if self.format == 'lmp-dump':
            backend = 'cpp'
        else:
            backend = 'py'
        return backend
    
    def __checkinput(self) -> str:
        '''
        do some check for the input type.

        Returns:
        -------
        str:
            The type of the input. Can be 'pymatgen', 'ase' or 'lmp-dump'.
        '''
        if isinstance(self.input, list):
            if isinstance(self.input[0], Structure):
                input_type = 'pymatgen'
            elif isinstance(self.input[0], Atoms):
                input_type = 'ase'
        elif isinstance(self.input, str):
            input_type = 'lmp-dump'
            if self.type_map is None:
                raise ValueError('type_map is required for lmp-dump file')
            if os.path.exists(self.input) is False:
                raise FileNotFoundError('The input file does not exist')
        else:
            raise ValueError('Invalid input type')

        return input_type

# ------------------- Python backend ------------------- #
def calculate_avgeraged_structure_py(input: list[Structure] | list[Atoms]) -> Structure:
    NotImplementedError('This method is not implemented yet')

def calculate_polarization_py(input: list[Structure] | list[Atoms]) -> np.ndarray:
    NotImplementedError('This method is not implemented yet')

def calculate_displacement_py(input: list[Structure] | list[Atoms]) -> np.ndarray:
    NotImplementedError('This method is not implemented yet')

def calculate_local_lattice_py(input: list[Structure] | list[Atoms]) -> np.ndarray:
    NotImplementedError('This method is not implemented yet')

def calculate_octahedron_rotation_py(input: list[Structure] | list[Atoms]) -> np.ndarray:
    NotImplementedError('This method is not implemented yet')

# ------------------- C++ backend ------------------- #
def calculate_averaged_structure_cpp(input: str, type_map: list[str], select: list[int]) -> Structure:
    cell, coord, types = get_averaged_structure(input, type_map, select)
    types = [type_map[t-1] for t in types]
    stru = Structure(Lattice(cell), types, coord, coords_are_cartesian=True)
    return stru

def calculate_displacement_cpp(input: str, select: list[int], nl: np.ndarray | str) -> np.ndarray:
    if isinstance(nl, str):
        nl = np.loadtxt(nl)
    elif isinstance(nl, np.ndarray):
        pass
    nl -= 1 # convert to 0-based index, for c++ backend
    disp = get_displacement(input, nl, select)
    disp = np.array(disp)
    return disp

def calculate_polarization_cpp(input: str, select: list[int], 
                               nl_ba: np.ndarray | str, nl_bx: np.ndarray | str, 
                               born_effective_charge: dict,
                               type_map: list[str]) -> np.ndarray:
    if isinstance(nl_ba, str):
        nl_ba = np.loadtxt(nl_ba)
    if isinstance(nl_bx, str):
        nl_bx = np.loadtxt(nl_bx)
    
    nl_ba -= 1
    nl_bx -= 1
    atomic_bec = []
    lmp = LAMMPSdump(input, type_map)
    stru = lmp.get_first_frame()
    for site in stru:
        atomic_bec.append(born_effective_charge[site.species_string])
    polar = get_polarization(input, nl_ba, nl_bx, atomic_bec, select)
    polar = np.array(polar)
    return polar

# ------------------- Other useful func ------------------- #
def convert_slice_to_list(input: str, select: slice) -> list[int]:
    '''
    convert slice object to list.
    If select is None, return a list of consecutive integers from nframes/2 to nframes.

    Parameters:
    ----------
    select: slice
        The slice object.

    Returns:
    -------
    list[int]:
        The list of integers, representing the index of selected frames.
    '''
    if select is not None:
        select = list(range(select.start, select.stop, select.step if select.step is not None else 1))
    else:
        if isinstance(input, list):
            nframes = len(input)
        else:
            lmp = LAMMPSdump(input)
            nframes = lmp.get_nframes()
        select = list(range(nframes//2, nframes))
    
    return select