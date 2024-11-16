import numpy as np
import matplotlib.pyplot as plt
from ase.geometry import get_layers
from ase.io import read
from ase import Atoms
import os

class ScalarPlotter:
    def __init__(self) -> None:
        pass 

     @staticmethod
    def __get_layers(stru: Atoms, 
                     element: list[str], 
                     tolerance: float,
                     axis: tuple[tuple]=((1, 0, 0), (0, 1, 0), (0, 0, 1))) -> tuple[np.ndarray, list[int]]:
        """
        Assign elements to layers based on the given axis

        Parameters:
        -----------
        stru: Atoms
            The input structure, which should be an `ASE` `Atoms` object
        element: list[str]
            Only the elements in the list will be considered
        tolerance: float
            The tolerance for the layer assignment
        axis: tuple[tuple]
            The axis for the layer assignment

        Returns:
        --------
        tag: np.ndarray
            The layer tag for each atom
        size: list[int]
            The size of the layers along the each axis
        """

        # 1. select the atoms with the given elements
        element_index = [idx for idx, i in enumerate(stru) if i.symbol in element]
        new_stru = stru[element_index]

        # 2. get the layer tag
        tag_x, _ = get_layers(new_stru, axis[0], tolerance=tolerance)
        tag_y, _ = get_layers(new_stru, axis[1], tolerance=tolerance)
        tag_z, _ = get_layers(new_stru, axis[2], tolerance=tolerance)

        # 3. combine the layer tag
        tag = np.concatenate((tag_x[:, np.newaxis], tag_y[:, np.newaxis], tag_z[:, np.newaxis]), axis=1)
        size = [len(set(tag_x)), len(set(tag_y)), len(set(tag_z))]

        return tag, size
    
    @staticmethod
    def __load_data(vector: str | np.ndarray,
                    tag: np.ndarray,
                    size: list[int]) -> np.ndarray:
        """
        Load vector data, and reshape it to the corresponding supercell shape.

        Parameters:
        -----------
        vector: str | np.ndarray
            The vector data, which should be a 2D array
        tag: np.ndarray
            The layer tag for each atom
        size: list[int]
            The size of the layers along the each axis
        
        Returns:
        --------
        data: np.ndarray
            The 4D array: in (nx, ny, nz, 3) format.
        """
        # 1. try to load the data
        if isinstance(vector, str) and vector.endswith(".npy"):
            raw_data = np.load(vector)
        elif isinstance(vector, np.ndarray):
            raw_data = vector.copy()
        else:
            raw_data = np.loadtxt(vector)
         
        # 2. check the shape of the data, should be 2D
        if raw_data.ndim != 2:
            raise ValueError("The shape of the vector should be 2D, but got {}".format(raw_data.shape))
        
        # 3. assign the data to the 3D array
        data = np.full((size[0], size[1], size[2], 3), np.nan)
        for i in range(len(raw_data)):
            data[tag[i, 0], tag[i, 1], tag[i, 2]] = raw_data[i]

        return data

    def plot(self) -> None:
        pass