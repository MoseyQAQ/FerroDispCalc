import numpy as np
from pymatgen.core import Structure,Lattice
from tqdm import tqdm

class Polar:
    def __init__(self, file_name:str,
                 type_map: list[str]=None,
                 natoms: int=5000,
                 nframes: int=2501,
                 fmt: str='traj') -> None:
        
        if fmt not in ['SingleFrame','traj']:
            raise NotImplementedError(f'Format {fmt} is not supported.')
        if fmt == 'SingleFrame' and nframes != 1:
            raise ValueError(f'Number of frames should be 1 for SingleFrame format.')
        
        print("========Info========")
        print("File name: ", file_name)
        print("Number of atoms: ", natoms)
        print("Number of frames: ", nframes)
        print("Type map: ", *type_map)
        print("Format: ", fmt)
        print("====================")
        self.file_name = file_name
        self.natoms = natoms
        self.nframes = nframes
        self.type_map = type_map
        self.fmt = fmt

    def parse_first_frame(self, st: Structure, ele:list[str]=['Ti'],
                          r:float=4.0, O_num:int=6,
                          vacancy: bool=False) -> tuple:
        ele_idx = []
        nn_idx = []
        for idx, site in enumerate(st):
            # only consider the specific element
            if str(site.specie) not in  ele:
                continue
            nn = st.get_neighbors(site,r)
            nn_idx.append([n.index for n in nn if str(n.specie) == 'O'])
            ele_idx.append(idx)
            # check if the number of oxygen atoms is correct
            if not vacancy and len(nn_idx[-1]) != O_num:
                raise ValueError(f'Number of oxygen atoms is not {O_num}.')
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
    def _get_one_frame_polar(self,prefix:str,
                             ele:str='Ti',
                             r:float=4.0,
                             O_num:int=6,
                             vacancy:bool=False) -> None:
        st = Structure.from_file(self.file_name)
        ele_idx, nn_idx = self.parse_first_frame(st,ele=ele,r=r,O_num=O_num,vacancy=vacancy)
        disp = np.zeros((self.nframes,len(ele_idx),3))
        origin = np.zeros((self.nframes,len(ele_idx),3))
        origin_,disp_ = self._get_disp(coords=st.frac_coords.copy(),
                                       cells=st.lattice.matrix.copy(),
                                       ele_idx=ele_idx,
                                       nn_idx=nn_idx,
                                       vacancy=vacancy,
                                       O_num=O_num)
        origin[0] = origin_
        disp[0] = disp_

        np.save(f'{prefix}origin.npy',origin)
        np.save(f'{prefix}disp.npy',disp)
    def _get_traj_polar(self,prefix:str,
                        ele:list[str]=['Ti'],
                        r:float=4.0,
                        O_num:int=6,
                        vacancy:bool=False) -> None:
        f = open(self.file_name, 'r')
        cell, type_index, coord = self._read_lmp_traj(f)
        st = Structure(Lattice(cell), type_index, coord,coords_are_cartesian=True)
        ele_idx, nn_idx = self.parse_first_frame(st, ele, r, O_num, vacancy)
        
        disp = np.zeros((self.nframes,len(ele_idx),3))
        origin = np.zeros((self.nframes,len(ele_idx),3))

        # read the rest of the frames
        for i in tqdm(range(self.nframes-1)):
            # deal with the last frame's data
            origin_, disp_ = self._get_disp(coord,cell,ele_idx,nn_idx,vacancy,O_num)
            origin[i] = origin_
            disp[i] = disp_

            # read in new frame
            cell, type_index, coord = self._read_lmp_traj(f)
        
        # deal with the last frame's data
        origin_, disp_ = self._get_disp(coord,cell,ele_idx,nn_idx,vacancy,O_num)
        origin[-1] = origin_
        disp[-1] = disp_

        # close the file
        f.close()

        # save file
        np.save(f'{prefix}origin.npy',origin)
        np.save(f'{prefix}disp.npy',disp)

    def get_polar(self,
                  prefix: str,
                  ele:list[str]=['Ti'],
                  r:float=4.0,
                  O_num:int=6,
                  vacancy: bool=False) -> None:
        
        if self.fmt == 'SingleFrame':
            self._get_one_frame_polar(prefix,ele=ele,r=r,O_num=O_num,vacancy=vacancy)
            return None 
        elif self.fmt == 'traj':
            self._get_traj_polar(prefix,ele=ele,r=r,O_num=O_num,vacancy=vacancy)
            return None
        