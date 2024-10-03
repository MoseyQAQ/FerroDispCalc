import numpy as np
from pymatgen.core import Structure, Lattice

class LAMMPSdump:
    def __init__(self, file_name: str, type_map: list[str] = None):
        '''
        initialize the LAMMPSdump object.

        Args:
            file_name (str): The file name of the lammps dump file.
            type_map (list[str]): The list of element names in the lammps dump file. Default is None.
        '''
        self.file_name = file_name
        self.type_map = type_map

    def get_first_frame(self) -> Structure:
        '''
        get the first frame in the lammps dump file.

        Returns:
        -------
        Structure: 
            The structure object of the first frame.
        '''
        f = open(self.file_name, 'r')
        self.natoms = self.__get_natoms()
        cell, type_index, coord = self._read_lmp_traj(f)
        f.close()
        stru = Structure(Lattice(cell), type_index, coord, coords_are_cartesian=True)
        return stru
    
    def get_nframes(self) -> int:
        '''
        get the number of frames in the lammps dump file.

        Parameters:
        ----------
        None

        Returns:
        -------
        int: 
            The number of frames.
        '''

        f = open(self.file_name, 'r')
        nframes = 0
        for line in f:
            if 'ITEM: TIMESTEP' in line:
                nframes += 1
        f.close()
        return nframes
    
    def __get_natoms(self) -> int:
        '''
        Get the number of atoms in the lammps dump file.

        Returns:
            int: The number of atoms.
        '''
        f = open(self.file_name, 'r')
        for line in f:
            if 'ITEM: NUMBER OF ATOMS' in line:
                natoms = int(f.readline().strip())
                break
        f.close()
        return natoms
    
    def __read_cell(self,f) -> np.ndarray:
        line = f.readline().split()
        line = [ float(x) for x in line ]
        xlo_bound = line[0]
        xhi_bound = line[1]
        xy = line[2]
        line = f.readline().split()
        line = [ float(x) for x in line ]
        ylo_bound = line[0]
        yhi_bound = line[1]
        xz = line[2]
        line = f.readline().split()
        line = [ float(x) for x in line ]
        zlo_bound = line[0]
        zhi_bound = line[1]
        yz = line[2]
        xlo = xlo_bound - min(0.0, xy, xz, xy+xz)
        xhi = xhi_bound - max(0.0, xy, xz, xy+xz)
        ylo = ylo_bound - min(0.0, yz)
        yhi = yhi_bound - max(0.0, yz)
        zlo = zlo_bound
        zhi = zhi_bound
        xx = xhi - xlo
        yy = yhi - ylo
        zz = zhi - zlo
        cell = np.zeros((3,3))
        cell[0,:] = [xx, 0, 0]
        cell[1,:] = [xy, yy, 0]
        cell[2,:] = [xz, yz, zz]
        return cell
    
    def __read_atoms(self,f) -> tuple:
        type_index = [None]*self.natoms
        coord = np.zeros((self.natoms,3))
        for i in range(self.natoms):
            line = f.readline().split()
            type_index[i] = self.type_map[int(line[1])-1]
            tmp = [ float(x) for x in line[2:] ]
            coord[i,0] = tmp[0]
            coord[i,1] = tmp[1]
            coord[i,2] = tmp[2]
        return type_index, coord
    
    def __skip_blank_line(self,f,n:int) -> None:
        for i in range(n):
            f.readline()

    def _read_lmp_traj(self,f) -> tuple:
        self.__skip_blank_line(f,5)
        cell = self.__read_cell(f)
        self.__skip_blank_line(f,1)
        type_index, coord = self.__read_atoms(f)
        return cell, type_index, coord