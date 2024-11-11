from ase.io import read
from ase import Atoms
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from ferrodispcalc.io.lammps import LAMMPSdump
import os

class Traj:
    def __init__(self,
                 input_stru: str | list[Structure] | list[Atoms] | Structure | Atoms,
                 format: str=None,
                 type_map: list[str]=None,) -> None:
        self.input_stru = input_stru
        self.format = format
        self.type_map = type_map
        self.type = None
        self.stru = None


    def __initialize_stru(self):

        # if input_stru is a string, try to read the file.
        if isinstance(self.input_stru, str): 

            # check if the file exists.
            if not os.path.exists(self.input_stru):
                raise FileNotFoundError(f'{self.input_stru} does not exist.')
            
            # read lammps
            if self.format == "lmp-dump":
                stru = LAMMPSdump(self.input_stru, self.type_map).get_first_frame()
                traj = [stru]
            elif self.format == "lmp-data":
                raise NotImplementedError('This feature is not implemented yet')
            
            # use ase or pymatgen to read the structure.
            else:
                try:
                    atom = read(self.input_stru, format=self.format, index=":")
                    traj = [AseAtomsAdaptor.get_structure(a) for a in atom]
                except Exception as e:
                    print(f'ASE Error: {e}; Trying to read the structure with pymatgen...')
                    try:
                        stru = Structure.from_file(self.input_stru)
                    except Exception as e:
                        print(f'Pymatgen Error: {e}')
                        raise ValueError('The input file format is not supported.')
                    
        elif isinstance(self.input_stru, list):
            pass

        elif isinstance(self.input_stru, Structure):
            traj = [self.input_stru]
        elif isinstance(self.input_stru, Atoms):
            traj = [AseAtomsAdaptor.get_structure(self.input_stru)]
        else:
            raise ValueError(f'Invalid input type: {type(self.input_stru)}')

        return traj
