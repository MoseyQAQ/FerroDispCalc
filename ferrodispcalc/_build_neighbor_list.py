import numpy as np
import os
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import read
from ase import Atoms
from ferrodispcalc.io.lammps import LAMMPSdump

class NeighborList:
    '''NeighborList class is used to build the neighbor list for a given atomic structure.

    Attributes:
    ----------
    input: str | pymatgen.core.Structure | ase.Atoms
        The input file name or the structure object.
    format: str
        The format of the input file. Default is None. For LAMMPS dump file, the format should be 'lmp-dump'.
    type_map: list[str]
        A list mapping the types of atoms as specified in the lammps dump file to their element symbols. Default is None.
    stru: pymatgen.core.Structure
        The structure object, used to build the neighbor list.
    nl: numpy.ndarray
        The neighbor list array where the first column is the index of the center element and the remaining columns are the indices of its neighbors. Indices are 1-based.

    Example:
    --------
    >>> from ferrodispcalc.build_neighbor_list import NeighborList
    >>> nl = NeighborList('POSCAR', format='vasp')
    >>> nl.build(['Pb', 'Sr'], ['O'], 4.0, 12)
    >>> nl.write('neighbor_list.dat')

    Methods:
    --------
    build(center_elements, neighbor_elements, cutoff, neighbor_num, defect)
        Constructs the neighbor list based on specified criteria.
    filter(nl, axis, rcut, neighbor_num)
        Filter the neighbor list based on the distance along the specified axis.
    write(output)
        Writes the constructed neighbor list to a file in a specified format.
    '''

    def __init__(self, input: str | Structure | Atoms, 
                 format: str=None, 
                 type_map: list[str]=None):
        '''
        Initializes the NeighborList object.

        Parameters
        ----------
        input : str | Structure | Atoms
            The input for creating the structure object. This could be a path to a file or a structure object
            from Pymatgen (`Structure`) or ASE (`Atoms`).
        format : str, optional
            The file format of the input if it is a path. If the input is a LAMMPS dump file, the format should be 'lmp-dump'.
        type_map : list[str], optional
            If using 'lmp-dump' format, which does not inherently contain atomic symbols, this should map types to elements.

        Raises
        ------
        FileNotFoundError
            If a file path is provided and the file does not exist.
        ValueError
            If the input format is not supported or cannot be properly read.
        '''
        self.input = input
        self.format = format
        self.type_map = type_map
        self.stru = self.__initialize_stru()

    def build(self, center_elements: list[str], neighbor_elements: list[str], 
              cutoff: float, neighbor_num: int, 
              defect: bool=False) -> np.ndarray:
        """
        Constructs the neighbor list based on specified criteria.

        Parameters
        ----------
        center_elements : list[str]
            Elements to consider as centers in the neighbor search.
        neighbor_elements : list[str]
            Elements to consider as neighbors in the search.
        cutoff : float
            The cutoff distance for considering an atom as a neighbor.
        neighbor_num : int
            The exact number of neighbors expected for each center element.
        defect : bool, optional
            If True, allows the number of neighbors to be different from `neighbor_num`.

        Returns
        -------
        numpy.ndarray
            The constructed neighbor list with 1-based indexing. The first column is the index of the center element and the remaining columns are the indices of its neighbors.

        Raises
        ------
        ValueError
            If the number of neighbors does not match `neighbor_num` and `defect` is False.
        """

        # initialize the index list
        center_elements_index = []
        neighbor_elements_index = []

        # use set to speed up the search
        center_elements_set = set(center_elements)
        neighbor_elements_set = set(neighbor_elements)

        # find the index of the center elements
        for idx in range(len(self.stru)):
            if str(self.stru[idx].specie) in center_elements_set:
                center_elements_index.append(idx)
        
        # build the neighbor list
        center_idx, point_idx, offset_vectors, distances = self.stru.get_neighbor_list(r=cutoff)

        # select the elements that are in the center_elements and neighbor_elements
        center_elements_mask = np.isin(center_idx, center_elements_index)
        neighbor_elements_mask = np.array([str(self.stru[idx].specie) in neighbor_elements_set for idx in point_idx])
        combined_mask = center_elements_mask & neighbor_elements_mask
        selected_center_elements_index = center_idx[combined_mask]
        selected_neighbor_elements_index = point_idx[combined_mask]
        selected_offset_vectors = offset_vectors[combined_mask]
        selected_distances = distances[combined_mask]

        # build the neighbor list in the format of {center: [neighbor1, neighbor2, ...]}
        result = {element_index: [] for element_index in center_elements_index}
        result_distance = {element_index: [] for element_index in center_elements_index}
        result_offset = {element_index: [] for element_index in center_elements_index}
        for center, point, offset, distance in zip(selected_center_elements_index, selected_neighbor_elements_index, selected_offset_vectors, selected_distances):
            result[center].append(point)
            result_distance[center].append(distance)
            result_offset[center].append(offset)
        
        # sort the neighbors by distance
        for center in center_elements_index:
            if len(result[center]) > 0:
                result[center] = np.array(result[center])
                result_distance[center] = np.array(result_distance[center])
                result_offset[center] = np.array(result_offset[center])
                result[center] = result[center][np.argsort(result_distance[center])]
        
        # check if the number of neighbors is correct
        # if defect is True, fill the missing neighbors with the center itself
        # if defect is False, raise an error
        # if the number of neighbors is more than neighbor_num, only keep the first neighbor_num neighbors
        for center in center_elements_index:
            if len(result[center]) < neighbor_num and not defect:
                raise ValueError(f"{center} {self.stru[center].specie} has {len(result[center])} neighbors, expected at least {neighbor_num}")
            elif len(result[center]) < neighbor_num and defect:
                print(f"Warning: {center} has {len(result[center])} neighbors, expected at least {neighbor_num}")
                neighbor_elements_index.append([center]*neighbor_num)
            elif len(result[center]) >= neighbor_num:
                neighbor_elements_index.append(result[center][:neighbor_num])
        
        center_elements_index = np.array(center_elements_index)
        neighbor_elements_index = np.array(neighbor_elements_index)
        self.nl = np.concatenate([center_elements_index[:,np.newaxis], neighbor_elements_index], axis=1)
        self.nl +=1 # convert the index to 1-based
        return self.nl
    
    def filter(self, nl: np.ndarray=None, axis: int=2, rcut: float=1, neighbor_num: int=3) -> np.ndarray:
        '''
        Filter the neighbor list based on the distance along the specified axis.
        It's useful for GaN system where the Ga atom is the center atom and the N atom is the neighbor atom.

        Parameters
        ----------
        nl : numpy.ndarray, optional
            The neighbor list to be filtered. If None, use the self.nl.
        axis : int, optional
            The axis along which the distance is calculated. Default is 2, which is the z-axis.
        rcut : float, optional
            The cutoff distance along the axis. Default is 1. Unit is Angstrom.
        neighbor_num : int, optional
            The number of neighbors to keep. Default is 3.
        
        Returns
        -------
        numpy.ndarray
            The filtered neighbor list with 1-based indexing. The first column is the index of the center element and the remaining columns are the indices of its neighbors.
        '''

        nl = self.nl.copy() if nl is None else nl
        center = nl[:,0]
        neighbors = nl[:,1:]
        center_pos = self.stru.cart_coords[center-1] # 0-based index
        new_nl = np.zeros((nl.shape[0], neighbor_num+1), dtype=int) # store the new neighbor list, 1-based index
        new_nl[:,0] = center
        for idx,n in enumerate(neighbors):
            neighbor_pos = self.stru.cart_coords[n-1] # 0-based index
            diff = neighbor_pos - center_pos[idx]

            # apply MIC
            box_inv = np.linalg.inv(self.stru.lattice.matrix)
            diff_frac = np.dot(diff, box_inv)
            diff_frac[diff_frac < -0.5] += 1
            diff_frac[diff_frac > 0.5] -= 1
            diff = np.dot(diff_frac, self.stru.lattice.matrix)
            diff = np.abs(diff[:,axis])
            mask = diff < rcut
            new_nl[idx,1:neighbor_num+1] = n[mask]

        self.nl = new_nl
        return self.nl
    
    def write(self, output: str):
        """
        Writes the neighbor list to a file.

        Parameters
        ----------
        output : str
            The file path where the neighbor list will be saved.

        Raises
        ------
        FileExistsError
            If the output file already exists to prevent accidental overwriting.
        """
        if os.path.exists(output):
            raise FileExistsError(f'{output} already exists.')
        else:
            np.savetxt(output, self.nl, fmt='%10d')

    def __initialize_stru(self):

        # initialize the structure object
        # if the input is a structure object, try to convert it to pymatgen structure
        if isinstance(self.input, str):
            if not os.path.exists(self.input):
                raise FileNotFoundError(f'{self.input} does not exist.')
            if self.format == 'lmp-dump':
                stru = LAMMPSdump(self.input, self.type_map).get_first_frame()
            else:
                try:
                    atom = read(self.input, format=self.format)
                    stru = AseAtomsAdaptor.get_structure(atom)
                except Exception as e:
                    print(f'ASE Error: {e}; Trying to read the structure with pymatgen...')
                    try:
                        stru = Structure.from_file(self.input)
                    except Exception as e:
                        print(f'Pymatgen Error: {e}')
                        raise ValueError('The input file format is not supported.')
        elif isinstance(self.input, Structure):
            stru = self.input
        elif isinstance(self.input, Atoms):
            stru = AseAtomsAdaptor.get_structure(self.input)
        return stru