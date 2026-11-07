"""Generate figures from repository inputs; all 200 trajectory frames are used.

Run after installing: python example/make_figures.py
No smoothing/downsampling. Bin width 0.02 Å, 600 dpi, maximum-normalized counts.
These are structural distributions, not an electronic-mechanism diagnosis.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from altbc_analyzer import ALTBCAnalyzer

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)

for name, source, mode, title in [
    ("poscar_ge", "POSCAR", "SingleSystem", "Ge-centered bond-pair distribution"),
    ("trajectory_ge", "test.xyz", "MultiSystem", "Ge-centered bond-pair distribution"),
]:
    analyzer = ALTBCAnalyzer(
        ROOT / "example" / source,
        filetype=mode,
        frame_interval=1,
        center_element="Ge",
        apply_mic=True,
        grid_size=0.02,
    )
    analyzer.compute_grid()
    analyzer.save_grid(OUT / f"{name}.npz")
    fig, ax = analyzer.plot_grid(OUT / f"{name}.png", dpi=600, title=title)
    # Preserve editable labels in companion vector formats.
    with matplotlib.rc_context(
        {"pdf.fonttype": 42, "svg.fonttype": "none", "font.family": "sans-serif", "font.size": 10}
    ):
        fig.savefig(OUT / f"{name}.pdf", facecolor="white")
        fig.savefig(OUT / f"{name}.svg", facecolor="white")
    (OUT / f"{name}.json").write_text(json.dumps(analyzer.stats, indent=2) + "\n")
    print(name, analyzer.stats, flush=True)
    plt.close(fig)
