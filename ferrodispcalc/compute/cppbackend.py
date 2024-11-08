from ferrodispcalc.compute.backend import ComputeBackend
from ase import Atoms
import numpy as np
from ferrodispcalc.core import get_averaged_structure

class CppCompute(ComputeBackend):
    def get_averaged_structure(self, select: list[int]) -> Atoms:
        data = get_averaged_structure(self.traj, self.type_map, select)

        cell = np.array(data[0])
        coords = np.array(data[1])
        symbols = data[2]

        return Atoms(
            symbols=symbols,
            positions=coords,
            cell=cell,
            pbc=True
        )

    def get_polarization(self, input, type_map, select, nl_ba, nl_bx, born_effective_charge):
        pass

    def get_displacement(self, input, type_map, select, nl):
        pass