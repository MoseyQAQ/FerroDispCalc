from build_neighbor_list import NeighborList
from type_map import UniPero
from ferrodispcalc.compute import Compute
nl = NeighborList(
    input='../test/BaTiO3/traj.lammpstrj',
    format='lmp-dump',
    type_map=UniPero
)

nl_ba=nl.build(
    center_elements=['Ti'],
    neighbor_elements=['Ba'],
    cutoff=4,
    neighbor_num=8,
    defect=False
)

nl_bx=nl.build(
    center_elements=['Ti'],
    neighbor_elements=['O'],
    cutoff=4,
    neighbor_num=6,
    defect=False
)
import numpy as np
c = Compute('../test/BaTiO3/traj.lammpstrj',format='lmp-dump',type_map=UniPero)
from ase.io import write
polar=c.get_polarization(slice(0,10),nl_ba,nl_bx,{'Ba': 2.77, 'Ti': 6.26, 'O': -3.01})
p = np.mean(polar,axis=0)
p = np.mean(p,axis=0)
print(p)
print(polar.shape)
