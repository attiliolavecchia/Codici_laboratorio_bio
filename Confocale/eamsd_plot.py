"""
Plot eaMSD (initial-displacement only) on linear axes and save to file.

Requirements implemented:
  - Read trajectories, compute eaMSD with no time-averaging (msd_analyzer).
  - Plot MSD vs time (seconds) on linear scale for both axes.
  - Labels:
      X: "Time Lag (τ) [seconds]"
      Y: "Ensemble-Averaged MSD (⟨r²⟩) [Length Units²]"
  - Title: "Ensemble-Averaged MSD (Initial Displacement)"
  - Save figure to eamsd_initial_plot.svg

"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from data_reader import read_trajectories_from_csv
from msd_analyzer import calculate_ensemble_msd


def plot_linear_and_save(tau, msd, output_path: Path, msd_sem=None) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    if msd_sem is not None:
        ax.errorbar(tau, msd, yerr=msd_sem, fmt='o-', color='C0',
                    capsize=4, capthick=1.2, elinewidth=1.2, label='eaMSD')
    else:
        ax.plot(tau, msd, marker='o', linestyle='-', color='C0', label='eaMSD')
    ax.set_ylabel('EA-MSD')
    ax.set_xlabel(r'Time Lag ($\tau$) [seconds]')
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


essential_description = (
    "Compute ensemble-averaged MSD (initial displacement only) and save a linear-scale plot."
)


def main() -> None:
    parser = argparse.ArgumentParser(description=essential_description)
    parser.add_argument("csv", type=str, help="Path to the trajectories CSV (e.g., sferette240nm_spots.csv)")
    parser.add_argument(
        "--max-lag-fraction",
        type=float,
        default=None,
        help=(
            "Fraction (0<f<=1] of the longest track to cap the maximum lag; if omitted, use full range (N_max-1)."
        ),
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output image filename (default: eamsd_plot_<csv_name>.svg, saved in eamsd_plots/ folder)",
    )
    args = parser.parse_args()

    trajectories = read_trajectories_from_csv(args.csv)
    result = calculate_ensemble_msd(trajectories, max_lag_fraction=args.max_lag_fraction)

    if result.tau.size == 0:
        print("No data to plot (empty or insufficient points).")
        return

    # Print requested metadata
    print(f"dt (global): {result.dt}")
    tau_max = result.tau[-1] if result.tau.size else float('nan')
    print(f"n_max (steps): {result.n_max}  |  tau_max (seconds): {tau_max}")
    print(f"Total trajectories (M): {result.total_trajectories}")
    print(f"Longest trajectory length (points): {result.longest_trajectory_points}")

    # Determine output path
    if args.output:
        output_path = Path(args.output)
        # Ensure parent directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        # Create output directory if it doesn't exist
        output_dir = Path("eamsd_plots")
        output_dir.mkdir(exist_ok=True)
        csv_name = Path(args.csv).stem
        output_filename = f"eamsd_plot_{csv_name}.svg"
        output_path = output_dir / output_filename
    
    plot_linear_and_save(result.tau, result.msd, output_path, msd_sem=result.msd_sem)
    print(f"Saved plot to: {output_path.resolve()}")


if __name__ == "__main__":
    main()
