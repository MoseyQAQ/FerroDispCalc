from build_neighbor_list import NeighborList
from type_map import UniPero

nl = NeighborList(
    input='Sc8_0.vasp',
    format='vasp',
    type_map=None
)

nl.build(
    center_elements=['Ga','Sc'],
    neighbor_elements=['N'],
    cutoff=4,
    neighbor_num=4,
    defect=False
)

nl.filter(nl=None,
          axis=2,
          rcut=1,
          neighbor_num=3)

nl.write("nl.dat")