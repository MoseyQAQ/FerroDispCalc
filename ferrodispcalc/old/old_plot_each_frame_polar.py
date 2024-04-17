from pymatgen.core import Structure
from dpdata import System
import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed
from glob import glob 
import matplotlib.pyplot as plt

type_map = ['Ba','Pb','Ca','Sr','Bi',
            'K','Na','Hf','Ti','Zr',
            'Nb','Mg','In','Zn','O']

def parse_first_frame(st: Structure,ele:str='Ti',
                      r:float=4.0, O_num:int=6,
                      vacancy: bool=False) -> tuple:
    ele_idx = []
    nn_idx = []


    for idx, site in enumerate(st):

        # only consider the specific element
        if str(site.specie) != ele:
            continue
        
        nn = st.get_neighbors(site,r)
        nn_idx.append([n.index for n in nn if str(n.specie) == 'O'])
        ele_idx.append(idx)

        # check if the number of oxygen atoms is correct
        if not vacancy and len(nn_idx[-1]) != O_num:
            raise ValueError(f'Number of oxygen atoms is not {O_num}.')
        
    return ele_idx, nn_idx


def get_polar(idx:int, st: Structure,ele_idx:list[int],nn_dix:list[list[int]],
              vacancy: bool=False, O_num: int=6) -> tuple:
    lattice = st.lattice.matrix
    disp = []
    origin = []

    for ele, nn in zip(ele_idx,nn_dix):
        xyz = np.array([0.0,0.0,0.0])

        # if vacancy, add a zero vector
        if len(nn) != O_num and vacancy:
            disp.append(np.array([0.0,0.0,0.0]))
            continue
        
        # calculate the center of the nn
        for n in nn:
            # check pbc
            for i in range(3):
                if st[n].coords[i] - st[ele].coords[i] > 5:
                    st[n].coords[i] -= lattice[i][i]
                elif st[n].coords[i] - st[ele].coords[i] < -5:
                    st[n].coords[i] += lattice[i][i]
            xyz += st[n].coords
        xyz /= len(nn)
        d = st[ele].coords-xyz
        disp.append(d)
        origin.append(st[ele].coords)

    return idx, np.array(origin), np.array(disp)

def main():

    pass

if __name__ == "__main__":
    #plot_unitcell_evolution(50)
    #plot_dxyz()
    #plot_distribution()
    #plot_phase_percentage()
    pass