import numpy as np
import matplotlib.pyplot as plt
from ase.geometry import get_layers
from ase.io import read 
from ase import Atoms
import os 

class VectorPoltter:
    def __init__(self, 
                 stru: str | Atoms, 
                 vector: str | np.ndarray,
                 element: list[str]=['Ti'],
                 tolerance: float=1.0,
                 axis: tuple[tuple]=((1, 0, 0), (0, 1, 0), (0, 0, 1))) -> None:
        
        self.stru = stru.copy() if isinstance(stru, Atoms) else read(stru) # use deepcopy to avoid changing the original structure
        self.tag, self.size = self.__get_layers(self.stru, element, tolerance, axis)
        print(f"Layer size: {self.size}")
        self.data = self.__load_data(vector, self.tag, self.size)
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
    
    @staticmethod
    def __cal_angle(dx, dy) -> np.ndarray:
        """
        calculate the angle of the vector

        Parameters:
        -----------
        dx: np.ndarray
            The x component of the vector
        dy: np.ndarray
            The y component of the vector
        
        Returns:
        --------
        angle: np.ndarray
            The angle of the vector
        """
        pp = np.sqrt(dx * dx + dy * dy)
        angle = np.arccos(dx / pp) / np.pi * 180.0
        index = np.where(dy < 0.0)
        angle[index] = 360.0 - angle[index]
        return angle
    
    def plot(self, dir: str, relative: bool=True, select: dict=None) -> None:
        '''
        Plot the vector field in the 2D plane

        Parameters:
        -----------
        dir: str
            The directory to save the plots
        relative: bool
            Whether to use the relative scale for the vector, default is True
        select: dict
            The selected index for each plane, default is None.
            If None, all the planes will be plotted.
            If not None, it should be a dictionary with the following format:
            {'xy': [0, 1, 2], 'xz': [0, 1, 2], 'yz': [0, 1, 2]}
            which means we only plot the 0, 1, 2 layers in the xy, xz, yz plane, respectively.
        '''
        if not os.path.exists(dir):
            os.makedirs(dir)

        quiver_kwargs = {}
        if relative == False:
            quiver_kwargs['scale'] = 1
            quiver_kwargs['scale_units'] = 'xy'

        # plot xy 
        for z_index in range(self.size[2]):
            if select is not None and 'xy' not in select:
                break
            if select is not None and z_index not in select['xy']:
                continue
            # save data to txt at z_index
            print(f"Plotting XY plane, {z_index}th layer")
            fig, ax = plt.subplots(figsize=(1*self.size[0], 1*self.size[1]))
            dx = self.data[:, :, z_index, 0].T
            dy = self.data[:, :, z_index, 1].T
            angle = self.__cal_angle(dx, dy)
            ax.quiver(dx, dy, **quiver_kwargs)
            sc = ax.imshow(angle, cmap='hsv', vmax=360, vmin=0, aspect=1.0, origin='lower')
            plt.title(f"XY plane, {z_index}th layer")
            plt.xlabel("[100]")
            plt.ylabel("[010]")
            plt.savefig(f"{dir}/XY_{z_index}.png", bbox_inches='tight', dpi=300)
            plt.close()
        # plot xz
        for y_index in range(self.size[1]):
            if select is not None and 'xz' not in select:
                break 
            if select is not None and y_index not in select['xz']:
                continue
            print(f"Plotting XZ plane, {y_index}th layer")
            fig, ax = plt.subplots(figsize=(1*self.size[0], 1*self.size[2]))
            dx = self.data[:, y_index, :, 0].T
            dz = self.data[:, y_index, :, 2].T
            angle = self.__cal_angle(dx, dz)
            ax.quiver(dx, dz, **quiver_kwargs)
            sc = ax.imshow(angle, cmap='hsv', vmax=360, vmin=0, aspect=1.0, origin='lower')
            plt.title(f"XZ plane, {y_index}th layer")
            plt.xlabel("[100]")
            plt.ylabel("[001]")
            plt.savefig(f"{dir}/XZ_{y_index}.png", bbox_inches='tight', dpi=300)
            plt.close()

        # plot yz 
        for x_index in range(self.size[0]):
            if select is not None and 'yz' not in select:
                break
            if select is not None and x_index not in select['yz']:
                continue
            print(f"Plotting YZ plane, {x_index}th layer")
            fig, ax = plt.subplots(figsize=(1*self.size[1], 1*self.size[2]))
            dy = self.data[x_index, :, :, 1].T
            dz = self.data[x_index, :, :, 2].T
            angle = self.__cal_angle(dy, dz)
            ax.quiver(dy, dz, **quiver_kwargs)
            sc = ax.imshow(angle, cmap='hsv', vmax=360, vmin=0, aspect=1.0, origin='lower')
            plt.title(f"YZ plane, {x_index}th layer")
            plt.xlabel("[010]")
            plt.ylabel("[001]")
            plt.savefig(f"{dir}/YZ_{x_index}.png", bbox_inches='tight', dpi=300)
            plt.close()