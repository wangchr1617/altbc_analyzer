"""Neighbor geometry using ASE row-vector cells and explicit image identities."""

from dataclasses import dataclass

import numpy as np
from ase.geometry import find_mic
from ase.neighborlist import neighbor_list
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class NeighborData:
    """Directed bonds: vector = r[j] - r[i] + shift @ cell, sorted by i.

    A neighbor is identified by (j, shift), not j alone. Nonzero images of
    the center itself are physical neighbors in small periodic cells.
    """

    i: np.ndarray
    j: np.ndarray
    distance: np.ndarray
    vector: np.ndarray
    shift: np.ndarray


class NeighborList:
    """Neighbors with strict distance < cutoff.

    apply_mic=True enables the input PBC flags and returns ALL images in
    range, including in small cells. False treats atoms as a finite cluster.
    apply_half is for pair calculations only; ALTBC requires full lists.
    """

    def __init__(self, cutoff, apply_mic=False, apply_half=False):
        if not np.isfinite(cutoff) or cutoff <= 0:
            raise ValueError("cutoff must be finite and positive")
        self.cutoff_neighbor = float(cutoff)
        self.cutoff_square = cutoff**2
        self.apply_mic = bool(apply_mic)
        self.apply_half = bool(apply_half)

    def get_determinant(self, box):
        return float(np.linalg.det(np.asarray(box).reshape(3, 3)))

    def get_area(self, a, b):
        return float(np.linalg.norm(np.cross(a, b)))

    def get_thickness(self, box):
        cell = np.asarray(box).reshape(3, 3)
        areas = np.linalg.norm(np.cross(cell[[1, 2, 0]], cell[[2, 0, 1]]), axis=1)
        if np.any(areas == 0) or abs(np.linalg.det(cell)) == 0:
            raise ValueError("cell must be nonsingular")
        return (abs(np.linalg.det(cell)) / areas).tolist()

    @staticmethod
    def apply_minimum_image_convention_one(value):
        return value - np.floor(value + 0.5)

    def apply_minimum_image_convention(self, box, x12, y12, z12):
        """True triclinic MIC, rather than fractional component rounding."""
        vector, _ = find_mic([x12, y12, z12], np.asarray(box).reshape(3, 3))
        return vector

    def geometry(self, atoms):
        positions = atoms.get_positions()
        if not np.all(np.isfinite(positions)) or not np.all(np.isfinite(atoms.cell.array)):
            raise ValueError("positions and cell must be finite")
        if self.apply_mic and np.any(atoms.pbc):
            periodic_cell = atoms.cell.array[atoms.pbc]
            if np.linalg.matrix_rank(periodic_cell) != len(periodic_cell):
                raise ValueError("periodic cell vectors must be nonzero and independent")
            i, j, d, vectors, shifts = neighbor_list(
                "ijdDS", atoms, self.cutoff_neighbor, self_interaction=False
            )
        else:
            pairs = cKDTree(positions).query_pairs(self.cutoff_neighbor, output_type="ndarray")
            i = np.concatenate((pairs[:, 0], pairs[:, 1]))
            j = np.concatenate((pairs[:, 1], pairs[:, 0]))
            vectors = positions[j] - positions[i]
            d = np.linalg.norm(vectors, axis=1)
            shifts = np.zeros((len(i), 3), dtype=int)
        keep = d < self.cutoff_neighbor
        if self.apply_half:
            positive = (shifts[:, 0] > 0) | (
                (shifts[:, 0] == 0)
                & ((shifts[:, 1] > 0) | ((shifts[:, 1] == 0) & (shifts[:, 2] > 0)))
            )
            keep &= (i < j) | ((i == j) & positive)
        if np.any(d[keep] <= 1e-12):
            raise ValueError("coincident atoms/images: bond angles are undefined")
        order = np.flatnonzero(keep)
        order = order[np.argsort(i[order], kind="stable")]
        return NeighborData(i[order], j[order], d[order], vectors[order], shifts[order])

    def find_neighbors(self, atoms):
        """Compatibility interface; geometry() also retains image identity."""
        data = self.geometry(atoms)
        return data.i, data.j, data.distance
