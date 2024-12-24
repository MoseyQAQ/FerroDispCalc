import numpy as np
import matplotlib.pyplot as plt
from ase.geometry import get_layers
from ase.io import read
from ase import Atoms
import seaborn as sns
import os
from typing import Union

class ScalarPlotter:
    def __init__(self,
                 stru: Union[str, Atoms],
                 scalar: Union[str, np.ndarray],
                 element: list[str]=['Ti'],
                 tolerance: float=1.0,
                 axis: tuple[tuple]=((1, 0, 0), (0, 1, 0), (0, 0, 1))) -> None:
        
        self.stru = stru.copy() if isinstance(stru, Atoms) else read(stru) # use deepcopy to avoid changing the original structure
        self.tag, self.size = self.__get_layers(self.stru, element, tolerance, axis)
        print(f"Layer size: {self.size}")
        self.data = self.__load_data(scalar, self.tag, self.size)
        print("Loaded data successfully")

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
    def __load_data(vector: Union[str, np.ndarray],
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
        if raw_data.ndim != 1:
            raise ValueError("The shape of the vector should be 1D, but got {}".format(raw_data.shape))
        
        # 3. assign the data to the 3D array
        data = np.full((size[0], size[1], size[2]), np.nan)
        for i in range(len(raw_data)):
            data[tag[i, 0], tag[i, 1], tag[i, 2]] = raw_data[i]

        return data

    def plot(self, dir: str, cmap: str="crest") -> None:
        if not os.path.exists(dir):
            os.mkdir(dir)
        
        # plot xy plane
        for z_index in range(self.size[2]):
            print(f"Plotting XY plane, {z_index}th layer")
            fig, ax = plt.subplots(figsize=(1*self.size[0], 1*self.size[1]))
            data = self.data[:, :, z_index]
            sns.heatmap(data, cmap=cmap, annot=True, fmt=".2f", ax=ax, cbar=True)
            plt.title(f"XY plane, {z_index}th layer")
            plt.xlabel("[001]")
            plt.ylabel("[010]")
            plt.tight_layout()
            plt.savefig(f"{dir}/XY_{z_index}.png", bbox_inches='tight', dpi=300)
            plt.close()
        
        # plot yz plane
        for x_index in range(self.size[0]):
            print(f"Plotting YZ plane, {x_index}th layer")
            fig, ax = plt.subplots(figsize=(1*self.size[1], 1*self.size[2]))
            data = self.data[x_index, :, :]
            sns.heatmap(data, cmap=cmap, annot=True, fmt=".2f", ax=ax, cbar=True)
            plt.title(f"YZ plane, {x_index}th layer")
            plt.xlabel("[010]")
            plt.ylabel("[100]")
            plt.tight_layout()
            plt.savefig(f"{dir}/YZ_{x_index}.png", bbox_inches='tight', dpi=300)
            plt.close()

        # plot xz plane
        for y_index in range(self.size[1]):
            print(f"Plotting XZ plane, {y_index}th layer")
            fig, ax = plt.subplots(figsize=(1*self.size[0], 1*self.size[2]))
            data = self.data[:, y_index, :]
            sns.heatmap(data, cmap=cmap, annot=True, fmt=".2f", ax=ax, cbar=True)
            plt.title(f"XZ plane, {y_index}th layer")
            plt.xlabel("[001]")
            plt.ylabel("[100]")
            plt.tight_layout()
            plt.savefig(f"{dir}/XZ_{y_index}.png", bbox_inches='tight', dpi=300)
            plt.close()

