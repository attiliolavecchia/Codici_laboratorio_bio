"""
Fit MSD data to one of four diffusion models and save a plot with results.

Models:
    linear           — MSD = 4D·τ
    nonlinear        — MSD = 4D·τ + v²·τ²
    anomalous        — MSD = 4D_α·τ^α
    anomalous_drift  — MSD = 4D_α·τ^α + v²·τ²

Usage:
    python fit_msd.py <csv> --model linear   [--max-lag-fraction F] [--fit-fraction F] [--output-dir DIR]
    python fit_msd.py <csv> --model nonlinear [--max-lag-fraction F] [--output-dir DIR]
    python fit_msd.py <csv> --model anomalous [--max-lag-fraction F] [--output-dir DIR]
    python fit_msd.py <csv> --model anomalous_drift [--max-lag-fraction F] [--output-dir DIR]
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from data_reader import read_trajectories_from_csv
from msd_analyzer import calculate_ensemble_msd
from msd_fitting import (
    fit_msd_linear, linear_msd_model,
    fit_msd_nonlinear, nonlinear_msd_model, analyze_velocities,
    fit_msd_anomalous, anomalous_msd_model,
    fit_msd_anomalous_drift, anomalous_drift_msd_model,
    calculate_r_squared,
)


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _plot_fit(
    tau_fit, msd_fit, msd_predicted, textstr: str, fit_label: str,
    output_path: Path, msd_sigma=None,
) -> None:
    """Generic fit plot: data + fit line + text box."""
    fig, ax = plt.subplots(figsize=(8, 6))
    if msd_sigma is not None:
        ax.errorbar(tau_fit, msd_fit, yerr=msd_sigma, fmt="o", color="C0",
                    markersize=8, alpha=0.7, capsize=4, capthick=1.2,
                    elinewidth=1.2, label="MSD Data", zorder=2)
    else:
        ax.plot(tau_fit, msd_fit, "o", color="C0", markersize=8, alpha=0.7,
                label="MSD Data", zorder=2)
    ax.plot(tau_fit, msd_predicted, "-", color="C3", linewidth=2.5,
            label=fit_label, zorder=3)
    ax.set_xlabel(r"Time Lag $\tau$ [s]", fontsize=12)
    ax.set_ylabel(r"MSD [$\mu$m$^2$]", fontsize=12)
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(loc="upper left", fontsize=10, framealpha=0.9)
    props = dict(boxstyle="round", facecolor="white", alpha=0.95,
                 edgecolor="black", linewidth=1.2)
    ax.text(0.98, 0.97, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment="top", horizontalalignment="right", bbox=props)
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Plot saved: {output_path.resolve()}")


# ---------------------------------------------------------------------------
# Per-model runners
# ---------------------------------------------------------------------------

def _run_linear(msd_result, args, output_dir: Path) -> None:
    fit = fit_msd_linear(
        msd_result.tau, msd_result.msd, msd_result.n_max, msd_result.dt,
        fit_fraction=args.fit_fraction,
        msd_sigma=msd_result.msd_sem,
    )
    print(f"\n  D  = ({fit.D:.4e} ± {fit.D_error:.2e}) μm²/s")
    print(f"  χ²_ν = {fit.chi_squared_red:.4f}")

    txt = "\n".join([
        r"$D = (%.2e \pm %.1e)\ \mu m^2/s$" % (fit.D, fit.D_error),
        r"$\chi^2_\nu = %.4f$" % fit.chi_squared_red,
    ])
    out = output_dir / f"{Path(args.csv).stem}_linear_fit.svg"
    _plot_fit(fit.tau_fit, fit.msd_fit, fit.msd_predicted, txt,
              r"Linear: MSD = 4D$\tau$", out, fit.msd_sigma_fit)


def _run_nonlinear(msd_result, trajectories, args, output_dir: Path) -> None:
    vstats = analyze_velocities(trajectories)
    fit = fit_msd_nonlinear(
        msd_result.tau, msd_result.msd, msd_result.n_max, msd_result.dt,
        velocity_stats=vstats, msd_sigma=msd_result.msd_sem,
    )
    print(f"\n  D  = ({fit.D:.4e} ± {fit.D_error:.2e}) μm²/s")
    print(f"  v  = ({fit.v:.4e} ± {fit.v_error:.2e}) μm/s")
    print(f"  χ²_ν = {fit.chi_squared_red:.4f}  (interval {fit.optimal_fraction:.0%})")

    txt = "\n".join([
        r"$D = (%.2e \pm %.1e)\ \mu m^2/s$" % (fit.D, fit.D_error),
        r"$v = (%.2e \pm %.1e)\ \mu m/s$" % (fit.v, fit.v_error),
        r"$\chi^2_\nu = %.4f$" % fit.chi_squared_red,
    ])
    out = output_dir / f"{Path(args.csv).stem}_nonlinear_fit.svg"
    _plot_fit(fit.tau_fit, fit.msd_fit, fit.msd_predicted, txt,
              r"Nonlinear: MSD = 4D$\tau$ + $v^2\tau^2$", out, fit.msd_sigma_fit)


def _run_anomalous(msd_result, args, output_dir: Path) -> None:
    fit = fit_msd_anomalous(
        msd_result.tau, msd_result.msd, msd_result.n_max, msd_result.dt,
        msd_sigma=msd_result.msd_sem,
    )
    print(f"\n  D_α = ({fit.D_alpha:.4e} ± {fit.D_alpha_error:.2e}) μm²/s^α")
    print(f"  α   = {fit.alpha:.4f} ± {fit.alpha_error:.4f}")
    print(f"  χ²_ν = {fit.chi_squared_red:.4f}  (interval {fit.optimal_fraction:.0%})")

    txt = "\n".join([
        r"$D_\alpha = (%.2e \pm %.1e)\ \mu m^2/s^\alpha$" % (fit.D_alpha, fit.D_alpha_error),
        r"$\alpha = %.4f \pm %.4f$" % (fit.alpha, fit.alpha_error),
        r"$\chi^2_\nu = %.4f$" % fit.chi_squared_red,
    ])
    out = output_dir / f"{Path(args.csv).stem}_anomalous_fit.svg"
    _plot_fit(fit.tau_fit, fit.msd_fit, fit.msd_predicted, txt,
              r"Anomalous: MSD = 4$D_\alpha \tau^\alpha$", out)


def _run_anomalous_drift(msd_result, trajectories, args, output_dir: Path) -> None:
    vstats = analyze_velocities(trajectories)
    fit = fit_msd_anomalous_drift(
        msd_result.tau, msd_result.msd, msd_result.n_max, msd_result.dt,
        velocity_stats=vstats, msd_sigma=msd_result.msd_sem,
    )
    print(f"\n  D_α = ({fit.D_alpha:.4e} ± {fit.D_alpha_error:.2e}) μm²/s^α")
    print(f"  α   = {fit.alpha:.4f} ± {fit.alpha_error:.4f}")
    print(f"  v   = ({fit.v:.4e} ± {fit.v_error:.2e}) μm/s")
    print(f"  χ²_ν = {fit.chi_squared_red:.4f}  (interval {fit.optimal_fraction:.0%})")

    txt = "\n".join([
        r"$D_\alpha = (%.2e \pm %.1e)\ \mu m^2/s^\alpha$" % (fit.D_alpha, fit.D_alpha_error),
        r"$\alpha = %.4f \pm %.4f$" % (fit.alpha, fit.alpha_error),
        r"$v = (%.2e \pm %.1e)\ \mu m/s$" % (fit.v, fit.v_error),
        r"$\chi^2_\nu = %.4f$" % fit.chi_squared_red,
    ])
    out = output_dir / f"{Path(args.csv).stem}_anomalous_drift_fit.svg"
    _plot_fit(fit.tau_fit, fit.msd_fit, fit.msd_predicted, txt,
              r"Anomalous+drift: MSD = 4$D_\alpha \tau^\alpha$ + $v^2\tau^2$", out)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Fit MSD data to a diffusion model.")
    parser.add_argument("csv", help="Path to the trajectories CSV file")
    parser.add_argument("--model", required=True,
                        choices=["linear", "nonlinear", "anomalous", "anomalous_drift"],
                        help="Diffusion model to fit")
    parser.add_argument("--max-lag-fraction", type=float, default=None,
                        help="Fraction of longest track to cap max lag for MSD calculation")
    parser.add_argument("--fit-fraction", type=float, default=0.10,
                        help="Fraction of lag steps for linear fit (default: 0.10)")
    parser.add_argument("--output-dir", type=str, default="fits",
                        help="Output directory for plots (default: fits/)")
    args = parser.parse_args()

    print(f"Reading: {args.csv}")
    trajectories = read_trajectories_from_csv(args.csv)
    print(f"Trajectories: {len(trajectories)}")

    msd_result = calculate_ensemble_msd(trajectories, max_lag_fraction=args.max_lag_fraction)
    if msd_result.tau.size == 0:
        print("Error: No MSD data.")
        return

    print(f"dt={msd_result.dt:.6f} s  n_max={msd_result.n_max}  "
          f"τ_max={msd_result.tau[-1]:.2f} s  M={msd_result.total_trajectories}")

    script_dir = Path(__file__).parent
    output_dir = script_dir / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.model == "linear":
        _run_linear(msd_result, args, output_dir)
    elif args.model == "nonlinear":
        _run_nonlinear(msd_result, trajectories, args, output_dir)
    elif args.model == "anomalous":
        _run_anomalous(msd_result, args, output_dir)
    elif args.model == "anomalous_drift":
        _run_anomalous_drift(msd_result, trajectories, args, output_dir)


if __name__ == "__main__":
    main()
