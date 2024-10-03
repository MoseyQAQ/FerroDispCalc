from ferrodispcalc.compute import Compute
from ferrodispcalc.type_map import UniPero
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io import write
c = Compute(input='../test/BaTiO3/traj.lammpstrj', format='lmp-dump', type_map=UniPero)

stru = c.get_averaged_structure(select=slice(0, 10))
print('finished')
#stru.to('POSCAR', 'stru.vasp')
atoms = AseAtomsAdaptor.get_atoms(stru)
write("1.vasp",atoms)