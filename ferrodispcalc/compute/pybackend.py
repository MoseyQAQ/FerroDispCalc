from ferrodispcalc.compute.backend import ComputeBackend
from ase import Atoms
import numpy as np

class PyCompute(ComputeBackend):
    def get_averaged_structure(self, traj: list[Atoms], select: list[int]) -> Atoms:

        # shape: (nframes, natoms, 3)
        coords = np.array([atoms.get_positions() for atoms in traj])
        cells = np.array([atoms.get_cell().array for atoms in traj])
        coords = coords[select]
        cells = cells[select]
        coords_frac = np.array([np.dot(coords[i], np.linalg.inv(cells[i])) for i in range(len(coords))])
        coords_diff = coords_frac - coords_frac[0]
        coords_frac[coords_diff > 0.5] -= 1
        coords_frac[coords_diff < -0.5] += 1
        coords = np.array([np.dot(coords_frac[i], cells[i]) for i in range(len(coords))])
        avg_cell = np.mean(cells, axis=0)
        avg_coords = np.mean(coords, axis=0)
        symbols = [atom.symbol for atom in traj[0]]
        return Atoms(
            symbols=symbols,
            positions=avg_coords,
            cell=avg_cell,
            pbc=True
        )
        
    def get_polarization(self, input, type_map, select, nl_ba, nl_bx, born_effective_charge):
        pass

    def get_displacement(self, traj: list[Atoms], select: list[int], nl: np.ndarray) -> np.ndarray:
        coords = np.array([atoms.get_positions() for atoms in traj])
        cells = np.array([atoms.get_cell().array for atoms in traj])
        coords = coords[select]
        cells = cells[select]
        nframes = coords.shape[0]
        natoms = nl.shape[0]
        nl -=1 # convert to 0-based index
        displacement = np.full((nframes, natoms, 3), np.nan) # shape: (nframes, nneighbors, 3)
        
        # walk through frames
        for i in range(nframes):
            # select center atoms and their coordinates
            center_id = nl[:, 0]
            center_coords = coords[i, center_id]
            cell = cells[i]

            # walk through neighbors
            for j, neighbors in enumerate(nl):
                neighbors_id = neighbors[1:]
                neighbors_coords = coords[i,neighbors_id]
                neighbors_coords_diff = center_coords[j] - neighbors_coords
                neighbors_coords_diff_frac = np.dot(neighbors_coords_diff, np.linalg.inv(cell))
                neighbors_coords_diff_frac[neighbors_coords_diff_frac > 0.5] -= 1
                neighbors_coords_diff_frac[neighbors_coords_diff_frac < -0.5] += 1
                neighbors_coords_diff = np.dot(neighbors_coords_diff_frac, cell)
                displacement[i, j] = np.mean(neighbors_coords_diff, axis=0)
        
        return displacement