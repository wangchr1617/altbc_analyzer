"""Independent exhaustive image enumeration, with no ASE neighbor/MIC routines."""

from itertools import permutations, product

import numpy as np
import pytest
from ase import Atoms

from altbc_analyzer import ALTBCAnalyzer, NeighborList


def reference_bonds(atoms, cutoff, periodic):
    ranges = [range(-4, 5) if periodic and p else [0] for p in atoms.pbc]
    result = {}
    for i, ri in enumerate(atoms.positions):
        for j, rj in enumerate(atoms.positions):
            for shift in product(*ranges):
                if i == j and shift == (0, 0, 0):
                    continue
                vector = rj - ri + np.array(shift) @ atoms.cell.array
                distance = np.linalg.norm(vector)
                if 1e-12 < distance < cutoff:
                    result[(i, j, *shift)] = (distance, vector)
    return result


def reference_triplets(atoms, analyzer):
    bonds = reference_bonds(atoms, analyzer.cutoff, analyzer.apply_mic)
    result = {}
    for center in range(len(atoms)):
        if analyzer.center_element and atoms[center].symbol != analyzer.center_element:
            continue
        neighbors = [
            (key, value)
            for key, value in bonds.items()
            if key[0] == center
            and (not analyzer.neighbor_element or atoms[key[1]].symbol == analyzer.neighbor_element)
        ]
        for (a, (r1, v1)), (c, (r2, v2)) in permutations(neighbors, 2):
            angle = np.rad2deg(np.arccos(np.clip(np.dot(v1, v2) / (r1 * r2), -1, 1)))
            if analyzer.theta_min - 1e-10 <= angle <= analyzer.theta_max + 1e-10:
                result[(a[1], center, c[1], *a[2:], *c[2:])] = (r1, r2, angle)
    return result


@pytest.mark.parametrize("pbc", [[True] * 3, [True, False, True], [False] * 3])
@pytest.mark.parametrize("periodic", [False, True])
@pytest.mark.parametrize("seed", [10, 42])
def test_exhaustive_skew_cell_neighbors_and_triplets(pbc, periodic, seed):
    rng = np.random.default_rng(seed)
    cell = np.array([[5.0, 0, 0], [2.4, 4.2, 0], [1.8, 0.7, 4.6]])
    atoms = Atoms("GeTeGeTeGe", scaled_positions=rng.random((5, 3)), cell=cell, pbc=pbc)
    analyzer = ALTBCAnalyzer("unused", cutoff=3.1, apply_mic=periodic, theta_min=75)
    expected = reference_bonds(atoms, analyzer.cutoff, periodic)
    bonds = NeighborList(analyzer.cutoff, periodic).geometry(atoms)
    actual = {
        (int(i), int(j), *map(int, shift)): (d, v)
        for i, j, d, v, shift in zip(bonds.i, bonds.j, bonds.distance, bonds.vector, bonds.shift)
    }
    assert len(actual) == len(bonds.i)  # no duplicate image identities
    assert actual.keys() == expected.keys()
    for key in actual:
        np.testing.assert_allclose(actual[key][0], expected[key][0], atol=1e-12)
        np.testing.assert_allclose(actual[key][1], expected[key][1], atol=1e-12)
    expected_tri = reference_triplets(atoms, analyzer)
    data = analyzer._compute_altbc(atoms)
    actual_tri = {
        (
            int(r.A),
            int(r.B),
            int(r.C),
            int(r.A_sx),
            int(r.A_sy),
            int(r.A_sz),
            int(r.C_sx),
            int(r.C_sy),
            int(r.C_sz),
        ): (r.AB, r.BC, r.angle_ABC)
        for r in data.itertuples()
    }
    assert len(actual_tri) == len(data)
    assert actual_tri.keys() == expected_tri.keys()
    for key in actual_tri:
        np.testing.assert_allclose(actual_tri[key], expected_tri[key], atol=2e-6)


def test_analytic_short_long_bonds():
    atoms = Atoms(
        "GeTe", positions=[[0, 0, 0], [2.7, 0, 0]], cell=[6.2, 9, 9], pbc=[True, False, False]
    )
    analyzer = ALTBCAnalyzer("unused", apply_mic=True, center_element="Ge", neighbor_element="Te")
    data = analyzer._compute_altbc(atoms).sort_values("AB")
    np.testing.assert_allclose(data[["AB", "BC", "angle_ABC"]], [[2.7, 3.5, 180], [3.5, 2.7, 180]])
    # The same unit-cell Te index represents two distinct physical neighbors.
    assert (data.A == data.C).all()
    assert (data.A_sx != data.C_sx).all()


def test_primitive_self_images_and_half_list():
    atoms = Atoms("Ge", cell=[3, 3, 3], pbc=True)
    assert len(NeighborList(3.1, True).geometry(atoms).i) == 6
    assert len(NeighborList(3.1, True, True).geometry(atoms).i) == 3
    assert len(ALTBCAnalyzer("unused", cutoff=3.1, apply_mic=True)._compute_altbc(atoms)) == 6


def test_no_cell_nonperiodic_and_strict_cutoff():
    atoms = Atoms("GeTeTe", positions=[[0, 0, 0], [1, 0, 0], [-1, 0, 0]])
    assert len(NeighborList(1, False).geometry(atoms).i) == 0
    assert len(NeighborList(1.01, False).geometry(atoms).i) == 4


@pytest.mark.parametrize("angle,expected", [(154.999, 0), (155, 2), (180, 2)])
def test_angle_boundaries(angle, expected):
    rad = np.deg2rad(angle)
    atoms = Atoms(
        "GeTeTe", positions=[[0, 0, 0], [2.7, 0, 0], [3.5 * np.cos(rad), 3.5 * np.sin(rad), 0]]
    )
    analyzer = ALTBCAnalyzer("unused", center_element="Ge")
    assert len(analyzer._compute_altbc(atoms)) == expected


def test_overlap_is_an_error():
    with pytest.raises(ValueError, match="coincident"):
        NeighborList(4).geometry(Atoms("GeTe", positions=np.zeros((2, 3))))


def test_unwrapped_rotation_and_translation():
    atoms = Atoms("GeTe", positions=[[0, 0, 0], [2.7, 0, 0]], cell=[6.2, 9, 9], pbc=True)
    analyzer = ALTBCAnalyzer("unused", apply_mic=True, center_element="Ge")
    expected = analyzer._compute_altbc(atoms)[["AB", "BC", "angle_ABC"]].sort_values("AB")
    atoms.positions[1] += np.array([3, -2, 1]) @ atoms.cell.array
    atoms.rotate(37, "z", rotate_cell=True)
    atoms.positions += [11, -7, 23]
    actual = analyzer._compute_altbc(atoms)[["AB", "BC", "angle_ABC"]].sort_values("AB")
    np.testing.assert_allclose(actual, expected, atol=2e-6)


@pytest.mark.parametrize(
    "kw",
    [
        dict(cutoff=0),
        dict(cutoff=np.nan),
        dict(grid_size=0),
        dict(theta_min=181),
        dict(xmin=4),
        dict(frame_interval=0),
        dict(frame_interval=1.5),
        dict(apply_half=True),
        dict(histogram_mode="wrong"),
        dict(filetype="wrong"),
        dict(grid_size=1e-9),
        dict(pair_block_size=0),
    ],
)
def test_invalid_parameters(kw):
    with pytest.raises(ValueError):
        ALTBCAnalyzer("unused", **kw)


def test_invalid_periodic_cell():
    with pytest.raises(ValueError, match="periodic cell"):
        NeighborList(4, True).geometry(Atoms("Ge", pbc=True))


def test_dense_real_poscar_boundary_centers():
    """All atoms in 27 independently enumerated images, sampled edge centers."""
    from pathlib import Path

    from ase.io import read

    atoms = read(Path(__file__).resolve().parents[1] / "example/POSCAR")
    x = ALTBCAnalyzer("unused", apply_mic=True, center_element="Ge", grid_size=0.02)
    bonds = NeighborList(x.cutoff, True).geometry(atoms)
    offsets = np.r_[0, np.cumsum(np.bincount(bonds.i, minlength=len(atoms)))]
    # Include atoms at different boundaries, not a convenient interior-only sample.
    ge = np.flatnonzero(atoms.numbers == 32)
    frac = atoms.get_scaled_positions()
    selected = set()
    for axis in range(3):
        selected.add(int(ge[np.argmin(frac[ge, axis])]))
        selected.add(int(ge[np.argmax(frac[ge, axis])]))
    selected.add(int(ge[len(ge) // 2]))
    shifts = np.array(list(product(range(-1, 2), repeat=3)))
    for center in selected:
        expected = {}
        for shift in shifts:
            vectors = atoms.positions - atoms.positions[center] + shift @ atoms.cell.array
            distances = np.linalg.norm(vectors, axis=1)
            for j in np.flatnonzero((distances > 1e-12) & (distances < x.cutoff)):
                expected[(j, *shift)] = distances[j]
        a, b = offsets[center : center + 2]
        actual = {
            (int(j), *map(int, s)): d
            for j, s, d in zip(bonds.j[a:b], bonds.shift[a:b], bonds.distance[a:b])
        }
        assert actual.keys() == expected.keys()
        for key in expected:
            assert actual[key] == pytest.approx(expected[key], abs=1e-12)
        assert len(actual) == 6
