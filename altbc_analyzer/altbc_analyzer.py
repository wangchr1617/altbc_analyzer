"""Angular-limited ordered bond-pair distributions, without implicit smoothing."""

import hashlib
import json
import warnings
from itertools import islice
from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import iread, read

from .neighbor_list import NeighborData, NeighborList

TRIPLET_COLUMNS = ["A", "B", "C", "AB", "BC", "angle_ABC"]
SHIFT_COLUMNS = ["A_sx", "A_sy", "A_sz", "C_sx", "C_sy", "C_sz"]


class ALTBC_Analyzer:
    """Near-collinear A-B-C triplets, B as center, both ordered contributions.

    Legacy defaults (including apply_mic=False) are retained. For bulk data
    set apply_mic=True, which respects the input PBC flags. Each neighbor is
    an atom-image identity, so primitive and repeated cells are consistent.

    histogram_mode='counts' counts all ordered triplets; 'occupancy' records
    at most one occurrence per bin per frame (old statistic, fixed geometry).
    Neither is a volume-normalized g^(3).
    """

    def __init__(
        self,
        filename,
        filetype="SingleSystem",
        frame_interval=10,
        center_element=None,
        xmin=2.5,
        xmax=3.8,
        grid_size=0.001,
        cutoff=4.0,
        theta_min=155,
        theta_max=180,
        apply_mic=False,
        apply_half=False,
        savegrid=False,
        *,
        neighbor_element=None,
        histogram_mode="counts",
        pair_block_size=256,
    ):
        if filetype not in {"SingleSystem", "MultiSystem"}:
            raise ValueError("filetype must be 'SingleSystem' or 'MultiSystem'")
        for name, value in [
            ("frame_interval", frame_interval),
            ("pair_block_size", pair_block_size),
        ]:
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        if not np.all(np.isfinite([xmin, xmax, grid_size, cutoff, theta_min, theta_max])):
            raise ValueError("distance and angle parameters must be finite")
        if xmin < 0 or xmax <= xmin or grid_size <= 0 or cutoff <= 0:
            raise ValueError("require 0 <= xmin < xmax, grid_size > 0, cutoff > 0")
        if not 0 <= theta_min <= theta_max <= 180:
            raise ValueError("require 0 <= theta_min <= theta_max <= 180 degrees")
        if apply_half:
            raise ValueError("ALTBC requires all neighbors: apply_half must be False")
        if histogram_mode not in {"counts", "occupancy"}:
            raise ValueError("histogram_mode must be 'counts' or 'occupancy'")
        self.filename, self.filetype = filename, filetype
        self.frame_interval = int(frame_interval)
        self.center_element, self.neighbor_element = center_element, neighbor_element
        self.nframe = 0
        self.xmin, self.xmax, self.grid_size = float(xmin), float(xmax), float(grid_size)
        ratio = (self.xmax - self.xmin) / self.grid_size
        if not np.isfinite(ratio) or ratio > 4096 + 1e-9:
            raise ValueError("grid exceeds 4096 bins per axis; increase grid_size")
        nearest = round(ratio)
        self.xbins = max(
            1,
            int(nearest if np.isclose(ratio, nearest, rtol=1e-12, atol=1e-12) else np.ceil(ratio)),
        )
        if self.xbins > 4096:
            raise ValueError("grid exceeds 4096 bins per axis; increase grid_size")
        self.bin_edges = np.minimum(
            self.xmin + np.arange(self.xbins + 1) * self.grid_size, self.xmax
        )
        self.bin_edges[-1] = self.xmax
        self.grid = np.zeros((self.xbins, self.xbins), dtype=np.int64)
        self.cutoff = float(cutoff)
        self.theta_min, self.theta_max = float(theta_min), float(theta_max)
        self.apply_mic, self.apply_half = bool(apply_mic), False
        self.savegrid, self.histogram_mode = savegrid, histogram_mode
        self.pair_block_size = int(pair_block_size)
        self.stats = {}
        self._grid_ready = False

    def _iter_structures(self):
        """Yield (zero-based source-frame index, atoms), without keeping frames."""
        self.nframe = 0
        if self.filetype == "SingleSystem":
            self.nframe = 1
            yield 0, read(self.filename)
        else:
            for selected, atoms in enumerate(
                islice(iread(self.filename), 0, None, self.frame_interval)
            ):
                self.nframe += 1
                yield selected * self.frame_interval, atoms

    def read_structures(self):
        """Compatibility method; compute methods use the streaming iterator."""
        return [atoms for _, atoms in self._iter_structures()]

    def generate_neighbor_list(self, atoms):
        return NeighborList(self.cutoff, self.apply_mic).find_neighbors(atoms)

    def _centers(self, atoms):
        if self.center_element is None:
            return np.arange(len(atoms))
        return np.flatnonzero(np.asarray(atoms.get_chemical_symbols()) == self.center_element)

    def _iter_pairs(self, atoms):
        """Selected unordered pairs; temporary matrices use bounded row blocks."""
        bonds = NeighborList(self.cutoff, self.apply_mic).geometry(atoms)
        if self.neighbor_element is not None:
            symbols = np.asarray(atoms.get_chemical_symbols())
            keep = symbols[bonds.j] == self.neighbor_element
            bonds = NeighborData(
                *(a[keep] for a in (bonds.i, bonds.j, bonds.distance, bonds.vector, bonds.shift))
            )
        offsets = np.concatenate(([0], np.cumsum(np.bincount(bonds.i, minlength=len(atoms)))))
        lower, upper = np.cos(np.deg2rad([self.theta_max, self.theta_min]))
        for center in self._centers(atoms):
            start, stop = offsets[center : center + 2]
            degree = stop - start
            if degree < 2:
                continue
            unit = bonds.vector[start:stop] / bonds.distance[start:stop, None]
            columns = np.arange(degree)[None, :]
            for row in range(0, degree, self.pair_block_size):
                end = min(row + self.pair_block_size, degree)
                cosine = np.clip(unit[row:end] @ unit.T, -1.0, 1.0)
                selected = (
                    (columns > np.arange(row, end)[:, None])
                    & (cosine >= lower - 1e-14)
                    & (cosine <= upper + 1e-14)
                )
                a, c = np.nonzero(selected)
                if len(a):
                    yield center, bonds, a + row + start, c + start, cosine[a, c]

    @staticmethod
    def _empty_frame():
        return pd.DataFrame(
            {
                name: pd.Series(dtype="float64" if name in {"AB", "BC", "angle_ABC"} else "int64")
                for name in TRIPLET_COLUMNS + SHIFT_COLUMNS
            }
        )

    def _compute_altbc(self, atoms):
        chunks = {name: [] for name in TRIPLET_COLUMNS + SHIFT_COLUMNS}
        for center, bonds, a, c, cosine in self._iter_pairs(atoms):
            left, right = np.concatenate((a, c)), np.concatenate((c, a))
            data = {
                "A": bonds.j[left],
                "B": np.full(len(left), center, dtype=int),
                "C": bonds.j[right],
                "AB": bonds.distance[left],
                "BC": bonds.distance[right],
                "angle_ABC": np.tile(np.rad2deg(np.arccos(cosine)), 2),
            }
            shifts = np.column_stack((bonds.shift[left], bonds.shift[right]))
            data.update({name: shifts[:, k] for k, name in enumerate(SHIFT_COLUMNS)})
            for name, values in data.items():
                chunks[name].append(values)
        if not chunks["A"]:
            return self._empty_frame()
        return pd.DataFrame({name: np.concatenate(values) for name, values in chunks.items()})

    def compute_altbc(self):
        """Materialize ordered triplets; use compute_grid for large trajectories.

        A/C are unit-cell indices; *_s* columns identify periodic images
        relative to B. frame is the zero-based source-file frame index.
        """
        chunks = []
        for frame, atoms in self._iter_structures():
            data = self._compute_altbc(atoms)
            data["frame"] = frame
            chunks.append(data)
        if chunks:
            return pd.concat(chunks, ignore_index=True)
        data = self._empty_frame()
        data["frame"] = pd.Series(dtype="int64")
        return data

    def _bin_pairs(self, x, y):
        # Snap numerical ties within 1e-12 Å to edges. Otherwise replication
        # and wrapping can move physically identical bonds between bins.
        tolerance = 1e-12
        valid = (
            (x >= self.xmin - tolerance)
            & (x <= self.xmax + tolerance)
            & (y >= self.xmin - tolerance)
            & (y <= self.xmax + tolerance)
        )

        def indices(values):
            scaled = (np.clip(values[valid], self.xmin, self.xmax) - self.xmin) / self.grid_size
            nearest = np.rint(scaled)
            scaled = np.where(
                np.abs(scaled - nearest) * self.grid_size <= tolerance, nearest, scaled
            )
            return np.minimum(np.floor(scaled).astype(np.int64), self.xbins - 1)

        xbin, ybin = indices(x), indices(y)
        return ybin * self.xbins + xbin

    def _compute_grid(self, df):
        """Compatibility helper: add one frame's table to the raw histogram."""
        bins = self._bin_pairs(df["AB"].to_numpy(), df["BC"].to_numpy())
        if self.histogram_mode == "occupancy":
            bins = np.unique(bins)
        ids, counts = np.unique(bins, return_counts=True)
        self.grid.ravel()[ids] += counts
        self._grid_ready = True

    def compute_grid(self):
        """Recompute raw counts from scratch, streaming one frame at a time.

        Histogram accumulation bypasses DataFrames and arccos entirely.
        Repeated calls and plots never double count or normalize raw data.
        """
        self._grid_ready = False
        self.grid.fill(0)
        self.stats = {}
        total = in_range = centers = 0
        for _, atoms in self._iter_structures():
            centers += len(self._centers(atoms))
            occupied = set() if self.histogram_mode == "occupancy" else None
            for _, bonds, a, c, _ in self._iter_pairs(atoms):
                total += 2 * len(a)
                ids = self._bin_pairs(bonds.distance[a], bonds.distance[c])
                reverse = (ids % self.xbins) * self.xbins + ids // self.xbins
                ids = np.concatenate((ids, reverse))
                in_range += len(ids)
                if occupied is not None:
                    occupied.update(ids.tolist())
                else:
                    unique, counts = np.unique(ids, return_counts=True)
                    self.grid.ravel()[unique] += counts
            if occupied:
                self.grid.ravel()[list(occupied)] += 1
        self.stats = {
            "frames": self.nframe,
            "center_observations": centers,
            "ordered_triplets": total,
            "in_range_triplets": in_range,
            "out_of_range_triplets": total - in_range,
            "histogram_mode": self.histogram_mode,
        }
        self._grid_ready = True
        if self.savegrid:
            path = Path(self.savegrid if isinstance(self.savegrid, (str, Path)) else "grid.txt")
            path.parent.mkdir(parents=True, exist_ok=True)
            np.savetxt(path, self.grid, fmt="%d", delimiter=",")
        return self.grid

    def normalized_grid(self, normalization="max"):
        """Transform a copy: density integrates to 1 over the displayed window.

        per_center is counts / summed selected centers / bin area. Frames
        are weighted by selected-center counts, rather than equally.
        """
        if not self._grid_ready:
            self.compute_grid()
        values = self.grid.astype(float)
        area = np.outer(np.diff(self.bin_edges), np.diff(self.bin_edges))
        if normalization == "max":
            denominator = values.max(initial=0)
        elif normalization == "density":
            denominator = values.sum() * area
        elif normalization == "per_center":
            if self.histogram_mode != "counts":
                raise ValueError("per_center normalization requires histogram_mode='counts'")
            denominator = self.stats.get("center_observations", 0) * area
        elif normalization == "none":
            return values
        else:
            raise ValueError("normalization must be max, density, per_center or none")
        return np.divide(values, denominator, out=np.zeros_like(values), where=denominator != 0)

    def save_grid(self, filename):
        """Compressed, pickle-free raw histogram with numerical provenance."""
        if not self._grid_ready:
            self.compute_grid()
        path = Path(filename)
        path.parent.mkdir(parents=True, exist_ok=True)
        metadata = {
            **self.stats,
            "source": str(self.filename),
            "filetype": self.filetype,
            "frame_interval": self.frame_interval,
            "center_element": self.center_element,
            "neighbor_element": self.neighbor_element,
            "cutoff": self.cutoff,
            "theta_min": self.theta_min,
            "theta_max": self.theta_max,
            "apply_mic": self.apply_mic,
            "grid_size": self.grid_size,
            "pair_convention": "ordered neighbor-image pairs",
            "version": "0.2.0",
        }
        source = Path(self.filename)
        if source.is_file():
            digest = hashlib.sha256()
            with source.open("rb") as stream:
                for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                    digest.update(chunk)
            metadata["source_sha256"] = digest.hexdigest()
        metadata["bin_edge_tolerance_angstrom"] = 1e-12
        metadata["angle_cosine_tolerance"] = 1e-14
        np.savez_compressed(
            path,
            grid=self.grid,
            bin_edges=self.bin_edges,
            metadata=json.dumps(metadata, ensure_ascii=False),
        )
        return path

    def plot_grid(
        self,
        output=None,
        *,
        dpi=400,
        normalization="max",
        cmap="Blues",
        ax=None,
        title=None,
        log_scale=False,
        recompute=False,
    ):
        """Render unsmoothed data, returning (figure, axes). Cached grid reused.

        Use output='figure.png' to save; recompute=True after changing input
        or analysis settings. Raw grid counts are never changed by plotting.
        """
        from .plotting import plot_grid

        if recompute or not self._grid_ready:
            self.compute_grid()
        return plot_grid(
            self,
            output=output,
            dpi=dpi,
            normalization=normalization,
            cmap=cmap,
            ax=ax,
            title=title,
            log_scale=log_scale,
        )

    def plot_scatter(self, output=None, *, dpi=400, ax=None, title=None):
        """Render every selected triplet (no hidden downsampling)."""
        from .plotting import plot_scatter

        data = self.compute_altbc()
        if len(data) > 1_000_000:
            warnings.warn(
                "Large scatter table; plot_grid avoids materializing triplets.", stacklevel=2
            )
        return plot_scatter(self, data, output=output, dpi=dpi, ax=ax, title=title)


ALTBCAnalyzer = ALTBC_Analyzer
