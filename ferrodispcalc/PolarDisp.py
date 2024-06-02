import numpy as np
from pymatgen.core import Structure,Lattice
from tqdm import tqdm

class PolarLMP:
    def __init__(self, file_name:str,
                 type_map: list[str]=None,
                 natoms: int=None,
                 nframes: int=None) -> None:
        
        self.file_name = file_name
        self.type_map = type_map

        if natoms is None:
            self.natoms = self._get_natoms()
        else:
            self.natoms = natoms
        if nframes is None:
            self.nframes = self._get_nframes()
        else:
            self.nframes = nframes
        
    
    def summary(self) -> None:
        print("=====================================")
        print(f"File name: {self.file_name}")
        print(f"Type map: {self.type_map}")
        print(f"Number of atoms: {self.natoms}")
        print(f"Number of frames: {self.nframes}")
        print("=====================================")

    def _get_natoms(self) -> int:
        flag = False
        with open(self.file_name, 'r') as f:
            for line in f:
                if 'ITEM: NUMBER OF ATOMS' in line:
                    flag = True
                    continue
                if flag:
                    natoms = int(line)
                    break
        return natoms

    def _get_nframes(self) -> int:
        with open(self.file_name, 'r') as f:
            nframes = 0
            while True:
                line = f.readline()
                if not line:
                    break
                if 'ITEM: TIMESTEP' in line:
                    nframes += 1
        return nframes
    def parse_first_frame(self, st: Structure, ele:list[str]=['Ti'], nn_ele:list[str]=['O'],
                          r:float=4.0, num:int=6,
                          vacancy: bool=False) -> tuple:
        ele_idx = []
        nn_idx = []

        ele_set = set(ele)
        nn_ele_set = set(nn_ele)

        # get element index
        for idx in range(len(st)):
            if str(st[idx].specie) in ele_set:
                ele_idx.append(idx)

        # build neighbor list
        center_idx, point_idx, offset_vectors, distances = st.get_neighbor_list(r)

        # select the neighbors of the element
        center_mask = np.isin(center_idx, ele_idx)
        species_mask = np.array([str(st[idx].specie) in nn_ele_set for idx in point_idx])
        combined_mask = center_mask & species_mask
        selected_center_idx = center_idx[combined_mask]
        selected_point_idx = point_idx[combined_mask]

        # build the neighbor list
        temp_result = {ele : [] for ele in ele_idx}
        for center, point in tqdm(zip(selected_center_idx, selected_point_idx)):
            temp_result[center].append(point)

        # check if the number of neighbors is correct
        for idx in ele_idx:
            if len(temp_result[idx]) != num and not vacancy:
                raise ValueError(f"Number of neighbors for {st[idx].specie, idx} is not correct")                
            else:
                nn_idx.append(temp_result[idx])
        
        return ele_idx, nn_idx
    
    def _read_cell(self,f) -> np.ndarray:
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
    def _read_atoms(self,f) -> tuple:
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
    def _skip_blank_line(self,f,n:int) -> None:
        for i in range(n):
            f.readline()
    def _read_lmp_traj(self,f) -> tuple:
        self._skip_blank_line(f,5)
        cell = self._read_cell(f)
        self._skip_blank_line(f,1)
        type_index, coord = self._read_atoms(f)
        return cell, type_index, coord
    def _get_disp(self,coords: np.ndarray,cells: np.ndarray,ele_idx: list[int], nn_idx: list[list[int]],
                  vacancy: bool=False, O_num: int=6) -> tuple:
        disp = np.zeros((len(ele_idx),3))
        origin = np.zeros((len(ele_idx),3))

        for idx,(ele,nn) in enumerate(zip(ele_idx,nn_idx)):
            xyz = np.array([0.0,0.0,0.0])
            # if vacancy, add a zero vector
            if len(nn) != O_num and vacancy:
                disp[idx] = np.array([0.0,0.0,0.0])
                origin[idx] = coords[ele]
                continue

            # calculate the center of the nn
            for n in nn:
                for i in range(3):
                    # check pbc
                    if coords[n][i] - coords[ele][i] > 5:
                        coords[n][i] -= cells[i][i]
                    elif coords[n][i] - coords[ele][i] < -5:
                        coords[n][i] += cells[i][i]
                xyz += coords[n]
            xyz /= len(nn)
            d = coords[ele]-xyz
            disp[idx] = d
            origin[idx] = coords[ele]

        return origin,disp
    def _get_polar(self,
                    ele:list[str]=['Ti'],
                    nn_ele:list[str]=['O'],
                    r:float=4.0,
                    num:int=6,
                    vacancy:bool=False) -> None:
        '''
        calculate the polar displacement

        Args:
            ele: center element to calculate the polar displacement
            nn_ele: neighbor element
            r: cutoff radius
            num: number of neighbors
            vacancy: whether to consider vacancy
        
        Returns:
            origin: origin coordinates
            disp: polar displacement
            cells: cell vectors
            coords: coordinates

        '''

        # read the first frame
        f = open(self.file_name, 'r')
        cell, type_index, coord = self._read_lmp_traj(f)
        st = Structure(Lattice(cell), type_index, coord,coords_are_cartesian=True)
        ele_idx, nn_idx = self.parse_first_frame(st, ele, nn_ele,r, num, vacancy)
        
        # initialize the arrays
        disp = np.zeros((self.nframes,len(ele_idx),3))
        origin = np.zeros((self.nframes,len(ele_idx),3))
        cells = np.zeros((self.nframes,3,3))
        coords = np.zeros((self.nframes,self.natoms,3))

        # read the rest of the frames
        for i in tqdm(range(self.nframes-1)):

            # deal with the last frame's data
            origin_, disp_ = self._get_disp(coord,cell,ele_idx,nn_idx,vacancy,num)
            cells[i] = cell
            coords[i] = coord
            origin[i] = origin_
            disp[i] = disp_

            # read in new frame
            cell, type_index, coord = self._read_lmp_traj(f)
  
        
        # deal with the last frame's data
        origin_, disp_ = self._get_disp(coord,cell,ele_idx,nn_idx,vacancy,num)
        origin[-1] = origin_
        disp[-1] = disp_
        cells[-1] = cell
        coords[-1] = coord

        # close the file
        f.close()

        return origin,disp,cells,coords

    def get_polar_displacement(self,
                  prefix: str,
                  ele:list[str]=['Ti'],
                  nn_ele:list[str]=['O'],
                  r:float=4.0,
                  num:int=6,
                  vacancy: bool=False,
                  save: bool=True) -> None:
        '''
        calculate the polar displacement, it will call the _get_polar function

        Args:
            prefix: prefix of the output file
            ele: center element to calculate the polar displacement
            nn_ele: neighbor element
            r: cutoff radius
            num: number of neighbors
            vacancy: whether to consider vacancy
            save: whether to save the output
        '''
        origin, disp, cells, coords = self._get_polar(ele=ele,nn_ele=nn_ele,r=r,num=num,vacancy=vacancy)
        
        if save:
            np.save(f'{prefix}_origin.npy', origin)
            np.save(f'{prefix}_disp.npy', disp)
    
    def get_avgeraged_structure(self,timeStep: int=500, output: str='out.xsf') -> None:
        '''
        calculate the averaged structure

        Args:
            timeStep: time step to average
            output: output file name
        '''
        f = open(self.file_name, 'r')
        cell, type_index, coord = self._read_lmp_traj(f)
        cells = np.zeros((self.nframes,3,3))
        coords = np.zeros((self.nframes,self.natoms,3))

        for i in tqdm(range(self.nframes-1)):
            cells[i] = cell
            coords[i] = coord
            cell, type_index, coord = self._read_lmp_traj(f)
        
        cells[-1] = cell
        coords[-1] = coord

        cells = np.mean(cells[timeStep:],axis=0)
        coords = np.mean(coords[timeStep:],axis=0)
        f.close()

        # write to the xsf file 
        self.write_xsf(coords,cells,type_index,output)

    def write_xsf(self, coord: np.ndarray, cell: np.ndarray, element: list[str], filename: str) -> None:
        '''
        write the xsf file

        Args:
            coord: coordinates
            cell: cell vectors
            element: element list
            filename: output file name
        '''
        with open(filename, 'w') as f:
            f.write(f'CRYSTAL\n')
            f.write(f'PRIMVEC\n')
            for i in range(3):
                f.write(f'{cell[i,0]:.6f} {cell[i,1]:.6f} {cell[i,2]:.6f}\n')
            f.write(f'PRIMCOORD\n')
            f.write(f'{len(coord)} 1\n')
            for i in range(len(coord)):
                f.write(f'{element[i]} {coord[i,0]:.6f} {coord[i,1]:.6f} {coord[i,2]:.6f}\n')