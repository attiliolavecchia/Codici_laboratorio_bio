"""
Plot eaMSD or taMSD with error bars and save to file.

Usage:
    python plot_msd.py <csv> --mode eamsd [--max-lag-fraction F] [--output PATH]
    python plot_msd.py <csv> --mode tamsd [--track-id ID] [--max-lag-fraction F] [--output PATH]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from data_reader import read_trajectories_from_csv, Trajectory
from msd_analyzer import calculate_ensemble_msd, calculate_time_averaged_msd_per_track


def plot_and_save(tau, msd, output_path: Path, msd_sem=None, label="MSD", color="C0") -> None:
    """Plot MSD vs τ on linear axes with optional SEM error bars."""
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    if msd_sem is not None:
        ax.errorbar(tau, msd, yerr=msd_sem, fmt="o-", color=color,
                    capsize=4, capthick=1.2, elinewidth=1.2, label=label)
    else:
        ax.plot(tau, msd, marker="o", linestyle="-", color=color, label=label)
    ax.set_xlabel(r"Time Lag ($\tau$) [seconds]")
    ax.set_ylabel(label)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend()
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot eaMSD or taMSD from a trajectories CSV.")
    parser.add_argument("csv", help="Path to the trajectories CSV file")
    parser.add_argument("--mode", choices=["eamsd", "tamsd"], required=True,
                        help="eamsd = ensemble-averaged MSD, tamsd = time-averaged MSD")
    parser.add_argument("--track-id", default=None,
                        help="Track ID for taMSD (default: first track by sorted ID)")
    parser.add_argument("--max-lag-fraction", type=float, default=0.25,
                        help="Fraction (0<f<=1] of track to cap max lag; default 0.25 (first 25%%)")
    parser.add_argument("--output", default=None,
                        help="Output image path (default: auto-generated in eamsd_plots/ or tamsd_plots/)")
    args = parser.parse_args()

    trajectories = read_trajectories_from_csv(args.csv)
    if not trajectories:
        print("No trajectories found in the CSV.")
        return

    csv_stem = Path(args.csv).stem

    if args.mode == "eamsd":
        result = calculate_ensemble_msd(trajectories, max_lag_fraction=args.max_lag_fraction)
        if result.tau.size == 0:
            print("No data to plot (empty or insufficient points).")
            return
        print(f"dt (global): {result.dt}")
        print(f"n_max (steps): {result.n_max}  |  tau_max: {result.tau[-1]:.4g} s")
        print(f"Trajectories: {result.total_trajectories}  |  Longest: {result.longest_trajectory_points} pts")

        out = Path(args.output) if args.output else Path("eamsd_plots") / f"eamsd_{csv_stem}.svg"
        plot_and_save(result.tau, result.msd, out, msd_sem=result.msd_sem, label="EA-MSD", color="C0")
        print(f"Saved: {out.resolve()}")

    else:  # tamsd
        if args.track_id is None:
            selected_id = sorted(trajectories.keys(), key=lambda k: str(k))[0]
            print(f"No --track-id given. Using first track: {selected_id}")
        else:
            selected_id = args.track_id
            if selected_id not in trajectories:
                if str(selected_id) in trajectories:
                    selected_id = str(selected_id)
                else:
                    print(f"Track '{args.track_id}' not found. Available (first 10): "
                          f"{sorted(trajectories.keys(), key=str)[:10]}")
                    return

        track: Trajectory = trajectories[selected_id]
        result = calculate_time_averaged_msd_per_track(track, max_lag_fraction=args.max_lag_fraction,
                                                       drift_corrected=False)
        if result.tau.size == 0:
            print("No data to plot.")
            return
        print(f"Track ID: {selected_id}  |  dt: {result.dt}")
        print(f"n_max (steps): {result.n_max}  |  tau_max: {result.tau[-1]:.4g} s")
        print(f"Trajectory length: {result.longest_trajectory_points} pts")

        tid = str(selected_id).replace(" ", "_")
        out = Path(args.output) if args.output else Path("tamsd_plots") / f"tamsd_{tid}.svg"
        plot_and_save(result.tau, result.msd, out, msd_sem=result.msd_sem, label="TA-MSD", color="C1")
        print(f"Saved: {out.resolve()}")


if __name__ == "__main__":
    main()
