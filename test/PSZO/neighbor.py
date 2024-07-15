from ferrodispcalc.build_neighbor_list import NeighborListABO3
from ferrodispcalc.typemap import UniPero
nl = NeighborListABO3(file_name='1.lammpstrj',
                      format='lmp-dump',
                      type_map=UniPero,)

nl.build(center_elements=['Pb','Sr'],
         neighbor_elements=['O'],
         cutoff=4.0,
         neighbor_num=12,
         defect=False)

nl.write('neighbor.dat')
