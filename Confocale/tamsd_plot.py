"""
Plot TAMSD (time-averaged MSD for a single trajectory) on linear axes and save to file.

- Reads trajectories from CSV (TrackMate-like headers supported by data_reader).
- Computes TAMSD using the function in msd_analyzer (no external averaging).
- Plots MSD vs time lag (seconds) on linear scale.

Labels:
  X: "Time Lag (τ) [seconds]"
  Y: "Time-Averaged MSD (⟨r²⟩) [Length Units²]"
Title: "Time-Averaged MSD (Single Trajectory)"

Usage example:
  python tamsd_plot.py sferette240nm_spots.csv --track-id <ID>

Dependencies: matplotlib
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt

from data_reader import read_trajectories_from_csv, Trajectory
from msd_analyzer import calculate_time_averaged_msd_per_track


def plot_linear_and_save(tau, msd, output_path: Path) -> None:
    plt.figure(figsize=(7.5, 5.0))
    plt.plot(tau, msd, marker="o", linestyle="-", color="C1", label="TAMSD")
    plt.ylabel(r"Time-Averaged MSD")
    plt.xlabel(r"Time Lag ($\tau$) [seconds]")
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


essential_description = (
    "Compute time-averaged MSD (single trajectory) and save a linear-scale plot."
)


def main() -> None:
    parser = argparse.ArgumentParser(description=essential_description)
    parser.add_argument("csv", type=str, help="Path to the trajectories CSV (e.g., sferette240nm_spots.csv)")
    parser.add_argument(
        "--track-id",
        type=str,
        default=None,
        help="Track ID to plot; if omitted, the first track by sorted ID is used.",
    )
    parser.add_argument(
        "--max-lag-fraction",
        type=float,
        default=None,
        help=(
            "Fraction (0<f<=1] of the selected track to cap maximum lag; if omitted, use full range (N-1)."
        ),
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output image filename (default: tamsd_plot_<track>.svg)",
    )
    args = parser.parse_args()

    trajectories = read_trajectories_from_csv(args.csv)
    if not trajectories:
        print("No trajectories found in the CSV.")
        return

    if args.track_id is None:
        # Choose the first track in sorted order for reproducibility
        selected_id = sorted(trajectories.keys(), key=lambda k: str(k))[0]
        print(f"No --track-id provided. Using first track: {selected_id}")
    else:
        selected_id = args.track_id
        if selected_id not in trajectories:
            # Attempt string coercion since keys are strings in reader
            if str(selected_id) in trajectories:
                selected_id = str(selected_id)
            else:
                print(f"Track ID '{args.track_id}' not found. Available IDs (first 10): "
                      f"{list(sorted(trajectories.keys(), key=lambda k: str(k)))[:10]}")
                return

    track: Trajectory = trajectories[selected_id]

    result = calculate_time_averaged_msd_per_track(track, max_lag_fraction=args.max_lag_fraction)

    if result.tau.size == 0:
        print("No data to plot for the selected trajectory (insufficient points).")
        return

    # Print metadata
    print(f"Track ID: {selected_id}")
    print(f"dt (track/global used): {result.dt}")
    tau_max = result.tau[-1] if result.tau.size else float('nan')
    print(f"n_max (steps): {result.n_max}  |  tau_max (seconds): {tau_max}")
    print(f"Trajectory length (points): {result.longest_trajectory_points}")

    # Resolve output path
    if args.output:
        output_path = Path(args.output)
        # Ensure parent directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        # Create default output directory if it doesn't exist
        output_dir = Path("tamsd_plots")
        output_dir.mkdir(exist_ok=True)
        output_filename = f"tamsd_plot_{str(selected_id).replace(' ', '_')}.svg"
        output_path = output_dir / output_filename

    plot_linear_and_save(result.tau, result.msd, output_path)
    print(f"Saved plot to: {output_path.resolve()}")


if __name__ == "__main__":
    main()
