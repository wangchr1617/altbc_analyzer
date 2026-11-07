"""Explicit axes, perceptually ordered colors and high-resolution exports."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import MaxNLocator

STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
}


def _axes(ax):
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 5.5), layout="constrained")
    else:
        fig = ax.figure
    return fig, ax


def _decorate(analyzer, ax, title):
    ax.set(
        xlabel=r"Bond length $r_1$ (Å)",
        ylabel=r"Bond length $r_2$ (Å)",
        xlim=(analyzer.xmin, analyzer.xmax),
        ylim=(analyzer.xmin, analyzer.xmax),
        aspect="equal",
    )
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.tick_params(direction="out", top=False, right=False)
    if title is None:
        title = f"{analyzer.center_element}-centered ALTBC" if analyzer.center_element else "ALTBC"
    ax.set_title(title, loc="left", weight="semibold", pad=12)


def _save(fig, output, dpi):
    if not np.isfinite(dpi) or dpi < 72:
        raise ValueError("dpi must be finite and >= 72")
    if output is not None:
        path = Path(output)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=dpi, facecolor="white", transparent=False)


def plot_grid(analyzer, *, output, dpi, normalization, cmap, ax, title, log_scale):
    values = analyzer.normalized_grid(normalization)
    positive = values[values > 0]
    if log_scale and positive.size:
        vmin, vmax = float(positive.min()), float(positive.max())
        # A single-valued positive grid still has a meaningful log colorbar.
        norm = LogNorm(vmin=vmin, vmax=vmax if vmax > vmin else vmin * 10)
    else:
        norm = Normalize(vmin=0, vmax=max(float(values.max(initial=0)), 1e-12))
    labels = {
        "max": "Relative intensity (maximum = 1)",
        "density": "Window probability density (Å⁻²)",
        "per_center": "Ordered pairs / center / Å²",
        "none": "Ordered triplet count",
    }
    if analyzer.histogram_mode == "occupancy":
        labels["none"] = "Occupied frame count"
        labels["density"] = "Window occupancy density (Å⁻²)"
    with mpl.rc_context(STYLE):
        fig, ax = _axes(ax)
        colors = mpl.colormaps[cmap].copy()
        colors.set_bad("white")
        data = np.ma.masked_where(values == 0, values)
        if np.allclose(np.diff(analyzer.bin_edges), np.diff(analyzer.bin_edges)[0]):
            artist = ax.imshow(
                data,
                origin="lower",
                extent=(analyzer.xmin, analyzer.xmax, analyzer.xmin, analyzer.xmax),
                cmap=colors,
                norm=norm,
                interpolation="nearest",
                rasterized=True,
            )
        else:
            artist = ax.pcolormesh(
                analyzer.bin_edges,
                analyzer.bin_edges,
                data,
                cmap=colors,
                norm=norm,
                shading="flat",
                rasterized=True,
            )
        colorbar = fig.colorbar(artist, ax=ax, fraction=0.047, pad=0.04)
        colorbar.set_label(labels[normalization], fontsize=9)
        colorbar.ax.tick_params(labelsize=8)
        colorbar.outline.set_linewidth(0.6)
        _decorate(analyzer, ax, title)
        if not positive.size:
            ax.text(
                0.5,
                0.5,
                "No triplets in the displayed window",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=9,
                color="0.35",
                wrap=True,
            )
        _save(fig, output, dpi)
    return fig, ax


def plot_scatter(analyzer, data, *, output, dpi, ax, title):
    with mpl.rc_context(STYLE):
        fig, ax = _axes(ax)
        ax.scatter(
            data["AB"], data["BC"], s=3, color="#315c8c", alpha=0.4, linewidths=0, rasterized=True
        )
        _decorate(analyzer, ax, title)
        if data.empty:
            ax.text(
                0.5,
                0.5,
                "No selected triplets",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=9,
            )
        _save(fig, output, dpi)
    return fig, ax
