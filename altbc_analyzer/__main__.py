"""Command-line analysis and high-resolution PNG export."""

import argparse
import json
import time
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compute ordered ALTBC histograms and save PNGs.")
    parser.add_argument("input", type=Path)
    parser.add_argument("-o", "--output", type=Path, default=Path("altbc.png"))
    parser.add_argument(
        "--trajectory", action="store_true", help="Read all frames (default stride 1)"
    )
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--center", help="Center element, e.g. Ge")
    parser.add_argument("--neighbor", help="Neighbor element, e.g. Te for Te-Ge-Te")
    parser.add_argument("--no-pbc", action="store_true", help="Treat input as a finite cluster")
    parser.add_argument("--cutoff", type=float, default=4.0)
    parser.add_argument("--theta-min", type=float, default=155)
    parser.add_argument("--theta-max", type=float, default=180)
    parser.add_argument("--xmin", type=float, default=2.5)
    parser.add_argument("--xmax", type=float, default=3.8)
    parser.add_argument("--bin-width", type=float, default=0.02, help="Histogram resolution in Å")
    parser.add_argument(
        "--dpi", type=int, default=600, help="Image resolution, independent of bin width"
    )
    parser.add_argument("--mode", choices=["counts", "occupancy"], default="counts")
    parser.add_argument(
        "--normalization", choices=["max", "density", "per_center", "none"], default="max"
    )
    parser.add_argument("--log", action="store_true", help="Log color scale, zero bins masked")
    args = parser.parse_args(argv)
    if args.output.suffix.lower() != ".png":
        parser.error("--output must have a .png suffix")
    if args.dpi < 72:
        parser.error("--dpi must be >= 72")
    if not args.trajectory and args.stride != 1:
        parser.error("--stride requires --trajectory")
    # Headless, and before pyplot is imported. Library calls never force a backend.
    import matplotlib

    matplotlib.use("Agg")
    from . import ALTBCAnalyzer

    start = time.perf_counter()
    try:
        analyzer = ALTBCAnalyzer(
            args.input,
            filetype="MultiSystem" if args.trajectory else "SingleSystem",
            frame_interval=args.stride,
            center_element=args.center,
            neighbor_element=args.neighbor,
            xmin=args.xmin,
            xmax=args.xmax,
            grid_size=args.bin_width,
            cutoff=args.cutoff,
            theta_min=args.theta_min,
            theta_max=args.theta_max,
            apply_mic=not args.no_pbc,
            histogram_mode=args.mode,
        )
        analyzer.compute_grid()
        compute_seconds = time.perf_counter() - start
        fig, _ = analyzer.plot_grid(
            args.output, dpi=args.dpi, normalization=args.normalization, log_scale=args.log
        )
        analyzer.save_grid(args.output.with_suffix(".npz"))
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    import matplotlib.pyplot as plt

    plt.close(fig)
    report = {
        **analyzer.stats,
        "compute_seconds": compute_seconds,
        "total_seconds": time.perf_counter() - start,
        "output": str(args.output),
        "dpi": args.dpi,
        "normalization": args.normalization,
        "log_scale": args.log,
    }
    args.output.with_suffix(".json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    if not analyzer.stats["in_range_triplets"]:
        print(
            "No triplets in the displayed window; check PBC, elements, angles and distance limits."
        )


if __name__ == "__main__":
    main()
