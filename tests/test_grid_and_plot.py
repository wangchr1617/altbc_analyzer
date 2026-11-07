import json
from pathlib import Path

import matplotlib
from ase import Atoms
from ase.io import write

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from PIL import Image

from altbc_analyzer import ALTBCAnalyzer
from altbc_analyzer.__main__ import main


@pytest.fixture
def chain(tmp_path):
    atoms = Atoms(
        "GeTe", positions=[[0, 0, 0], [2.7, 0, 0]], cell=[6.2, 9, 9], pbc=[True, False, False]
    )
    path = tmp_path / "chain.extxyz"
    write(path, [atoms, atoms.repeat((2, 1, 1)), atoms])
    return path


def analyzer(path, **kw):
    return ALTBCAnalyzer(path, center_element="Ge", apply_mic=True, grid_size=0.02, **kw)


def test_counts_symmetry_repeat_and_density(chain, monkeypatch):
    x = analyzer(chain, filetype="MultiSystem", frame_interval=1)
    # Grid path must not construct any triplet DataFrame.
    monkeypatch.setattr(x, "_compute_altbc", lambda _: pytest.fail("table created"))
    grid = x.compute_grid().copy()
    assert grid.sum() == 8
    np.testing.assert_array_equal(grid, grid.T)
    np.testing.assert_array_equal(x.compute_grid(), grid)
    area = np.outer(np.diff(x.bin_edges), np.diff(x.bin_edges))
    assert np.sum(x.normalized_grid("density") * area) == pytest.approx(1)
    assert np.sum(x.normalized_grid("per_center") * area) == pytest.approx(2)
    assert x.stats["center_observations"] == 4
    assert x.stats["frames"] == 3
    assert x.stats["out_of_range_triplets"] == 0
    np.testing.assert_array_equal(x.grid, grid)


def test_occupancy_is_per_frame_not_per_center(chain):
    x = analyzer(chain, filetype="MultiSystem", frame_interval=1, histogram_mode="occupancy")
    assert x.compute_grid().sum() == 6
    with pytest.raises(ValueError, match="requires"):
        x.normalized_grid("per_center")


def test_stride_and_table_histogram_agree(chain):
    x = analyzer(chain, filetype="MultiSystem", frame_interval=2, neighbor_element="Te")
    data = x.compute_altbc()
    assert sorted(data.frame.unique()) == [0, 2]
    assert len(data) == 4
    expected, _, _ = np.histogram2d(data.BC, data.AB, bins=(x.bin_edges, x.bin_edges))
    np.testing.assert_array_equal(x.compute_grid(), expected)
    assert data.A.dtype.kind == "i"


def test_bin_boundaries_and_nonintegral_width():
    x = ALTBCAnalyzer("unused", xmin=1, xmax=2, grid_size=0.3)
    data = pd.DataFrame({"AB": [0.999, 1, 1.3, 2, 2.001, 1.9], "BC": [1, 1, 1.3, 2, 1, 1.9]})
    x._compute_grid(data)
    expected, _, _ = np.histogram2d(data.BC, data.AB, bins=(x.bin_edges, x.bin_edges))
    np.testing.assert_array_equal(x.grid, expected)
    assert x.grid.sum() == 4
    assert x.bin_edges[-1] == 2
    assert np.diff(x.bin_edges)[-1] == pytest.approx(0.1)


def test_primitive_supercell_per_center_invariance(chain, tmp_path):
    from ase.io import read

    primitive = read(chain, index=0)
    grids = []
    for k in (1, 4):
        path = tmp_path / f"{k}.xyz"
        write(path, primitive.repeat((k, 1, 1)))
        x = analyzer(path)
        x.compute_grid()
        grids.append(x.normalized_grid("per_center"))
    np.testing.assert_array_equal(*grids)


@pytest.mark.parametrize("log_scale", [False, True])
def test_png_raw_counts_axes_and_empty(chain, tmp_path, log_scale):
    x = analyzer(chain)
    raw = x.compute_grid().copy()
    original_rc = matplotlib.rcParams.copy()
    path = tmp_path / "sub" / "plot.png"
    fig, ax = x.plot_grid(path, dpi=400, log_scale=log_scale)
    assert ax.figure is fig
    assert ax.get_aspect() == 1
    with Image.open(path) as im:
        assert im.size == (2560, 2200)
        assert im.info["dpi"][0] == pytest.approx(400, abs=0.1)
    fig2, _ = x.plot_grid()
    np.testing.assert_array_equal(x.grid, raw)
    assert matplotlib.rcParams["font.size"] == original_rc["font.size"]
    plt.close(fig)
    plt.close(fig2)
    x.center_element = "Xe"
    fig, _ = x.plot_grid(tmp_path / "empty.png", log_scale=log_scale, recompute=True)
    assert np.isfinite(x.normalized_grid()).all()
    assert x.compute_altbc().empty
    assert "AB" in x.compute_altbc()
    plt.close(fig)


def test_npz_and_cli(chain, tmp_path, capsys):
    out = tmp_path / "cli.png"
    main(
        [
            str(chain),
            "--trajectory",
            "--center",
            "Ge",
            "--neighbor",
            "Te",
            "--stride",
            "2",
            "-o",
            str(out),
            "--dpi",
            "150",
        ]
    )
    assert out.is_file()
    with np.load(out.with_suffix(".npz"), allow_pickle=False) as archive:
        assert archive["grid"].sum() == 4
        meta = json.loads(str(archive["metadata"]))
        assert meta["frames"] == 2
        assert meta["apply_mic"] is True
        assert meta["neighbor_element"] == "Te"
    summary = json.loads(out.with_suffix(".json").read_text())
    assert summary["in_range_triplets"] == 4
    assert summary["total_seconds"] > 0
    assert "ordered_triplets" in capsys.readouterr().out


def test_repo_example_against_independent_dense_distances():
    # Independent O(N^2) finite-cluster reference on the real first frame.
    # Avoid the broken historical cell list as the correctness oracle.
    from ase.io import read

    root = Path(__file__).resolve().parents[1]
    atoms = read(root / "example/test.xyz", index=0)
    x = ALTBCAnalyzer("unused", center_element="Ge", apply_mic=False, grid_size=0.02)
    positions = atoms.positions
    r1, r2 = [], []
    for i in np.flatnonzero(atoms.numbers == 32):
        v = positions - positions[i]
        d = np.linalg.norm(v, axis=1)
        neighbors = np.flatnonzero((d > 1e-12) & (d < x.cutoff))
        for a in neighbors:
            for c in neighbors:
                if a == c:
                    continue
                angle = np.rad2deg(np.arccos(np.clip(np.dot(v[a], v[c]) / (d[a] * d[c]), -1, 1)))
                if x.theta_min <= angle <= x.theta_max:
                    r1.append(d[a])
                    r2.append(d[c])
    data = x._compute_altbc(atoms)
    assert len(data) == len(r1)
    ref = np.histogram2d(r2, r1, bins=(x.bin_edges, x.bin_edges))[0]
    x._compute_grid(data)
    np.testing.assert_array_equal(x.grid, ref)
