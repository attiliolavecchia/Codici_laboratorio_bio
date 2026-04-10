"""
Ergodicity check — compare EA-MSD and ensemble-averaged TA-MSD.

For an ergodic process the two MSD estimators must converge to the same
curve.  This script overlays them on a single plot for each CSV file in
the non-anomalous dataset, at lag fractions of 10 % and 25 %.

Output:
    Results/no_anomalous/ergodicity/
        <stem>_ergodicity_f010.svg   — overlay plot (10 % lag fraction)
        <stem>_ergodicity_f025.svg   — overlay plot (25 % lag fraction)
        <stem>_ergodicity_f050.svg   — overlay plot (50 % lag fraction)
        <stem>_eb_f010.svg           — EB parameter plot (10 %)
        <stem>_eb_f025.svg           — EB parameter plot (25 %)
        <stem>_eb_f050.svg           — EB parameter plot (50 %)

Usage:
    python check_ergodicity.py [--compare-fraction 0.10]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Suppress the harmless ddof warning when few tracks contribute at large lags
warnings.filterwarnings("ignore", message="Degrees of freedom <= 0 for slice")

from data_reader import read_trajectories_from_csv
from msd_analyzer import (
    calculate_ensemble_msd,
    calculate_time_averaged_msd_per_track,
    average_across_trajectories,
    build_tau_array,
    determine_maximum_lag_steps,
    estimate_global_time_step,
)

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "Data" / "31_10_no_anomalous"
OUT_DIR = SCRIPT_DIR / "Results" / "no_anomalous" / "ergodicity"

LAG_FRACTIONS = [0.10, 0.25, 0.50]  # Lag fractions to test (10%, 25%, 50%)


def compute_ensemble_tamsd(trajectories, *, max_lag_fraction=None, global_dt=None):
    """Compute the ensemble-averaged TA-MSD: average TA-MSD across all tracks.

    For each trajectory we compute the full time-averaged MSD, then we
    average those curves across trajectories (same procedure as the EA-MSD
    averaging, but starting from TA-MSD per track instead of single initial
    displacements).

    Returns (tau, msd_mean, msd_sem, per_track_matrix) where
    *per_track_matrix* is a (M, K) float array of individual TA-MSD curves
    (NaN-padded for tracks shorter than K).
    """
    if not trajectories:
        raise ValueError("No trajectories")

    max_points = max(t.n_points for t in trajectories.values())
    eff_frac = 1.0 if max_lag_fraction is None else float(max_lag_fraction)
    K = determine_maximum_lag_steps(max_points, eff_frac)

    if global_dt is None:
        global_dt = estimate_global_time_step(trajectories)

    tau = build_tau_array(K, global_dt)

    per_track_tamsd = []
    for track in trajectories.values():
        res = calculate_time_averaged_msd_per_track(
            track, max_lag_fraction=max_lag_fraction, dt_override=global_dt,
        )
        # Pad / trim to K entries
        arr = np.full(K, np.nan, dtype=float)
        n = min(K, res.msd.size)
        arr[:n] = res.msd[:n]
        per_track_tamsd.append(arr)

    per_track_matrix = np.array(per_track_tamsd)  # (M, K)

    msd_mean, msd_std, counts = average_across_trajectories(per_track_tamsd)
    with np.errstate(invalid="ignore"):
        msd_sem = msd_std / np.sqrt(counts.astype(float))

    return tau, msd_mean, msd_sem, per_track_matrix


def compute_eb_parameter(per_track_matrix):
    """Ergodicity breaking parameter EB(τ).

    EB(τ) = [⟨δ²(τ)²⟩ − ⟨δ²(τ)⟩²] / ⟨δ²(τ)⟩²

    i.e. the relative variance of the per-track TA-MSD values at each lag.
    For an ergodic process EB → 0 when observation time T ≫ τ.

    Parameters
    ----------
    per_track_matrix : (M, K) ndarray
        Individual TA-MSD curves (NaN for missing lags).

    Returns
    -------
    eb : (K,) ndarray
        EB parameter at each lag (NaN where fewer than 2 tracks contribute).
    n_tracks : (K,) ndarray of int
        Number of tracks contributing at each lag.
    """
    with np.errstate(invalid="ignore"):
        mean = np.nanmean(per_track_matrix, axis=0)        # ⟨δ²⟩
        mean_sq = np.nanmean(per_track_matrix ** 2, axis=0) # ⟨δ²²⟩
        eb = (mean_sq - mean ** 2) / (mean ** 2)

    n_tracks = np.sum(~np.isnan(per_track_matrix), axis=0)
    eb[n_tracks < 2] = np.nan
    return eb, n_tracks


def plot_eb_parameter(tau, eb, n_tracks, output_path, title=""):
    """Plot EB(τ) vs τ with a secondary axis showing contributing tracks."""
    fig, ax1 = plt.subplots(figsize=(8, 5.5))

    valid = ~np.isnan(eb)
    ax1.plot(tau[valid], eb[valid], "o-", color="C2", markersize=4,
             label="EB(τ)", zorder=3)
    ax1.axhline(0, color="grey", linewidth=0.8, linestyle="--", alpha=0.6)
    ax1.set_xlabel(r"Time lag $\tau$ [s]", fontsize=12)
    ax1.set_ylabel(r"EB($\tau$)", fontsize=12, color="C2")
    ax1.tick_params(axis="y", labelcolor="C2")
    ax1.grid(True, linestyle=":", alpha=0.5)

    # Secondary axis: number of tracks
    ax2 = ax1.twinx()
    ax2.fill_between(tau, 0, n_tracks, alpha=0.15, color="C7", step="mid")
    ax2.set_ylabel("# tracks", fontsize=11, color="C7")
    ax2.tick_params(axis="y", labelcolor="C7")

    if title:
        ax1.set_title(title, fontsize=11)

    ax1.legend(loc="upper right", fontsize=10)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_ergodicity(tau_ea, msd_ea, sem_ea,
                    tau_ta, msd_ta, sem_ta,
                    compare_n, output_path, title=""):
    """Overlay EA-MSD and ensemble-averaged TA-MSD with comparison region shaded."""
    fig, ax = plt.subplots(figsize=(8, 5.5))

    ax.errorbar(tau_ea, msd_ea, yerr=sem_ea, fmt="o-", color="C0",
                capsize=3, capthick=1, elinewidth=1, markersize=5,
                label="EA-MSD", alpha=0.85, zorder=3)
    ax.errorbar(tau_ta, msd_ta, yerr=sem_ta, fmt="s-", color="C3",
                capsize=3, capthick=1, elinewidth=1, markersize=5,
                label="TA-MSD (ensemble avg)", alpha=0.85, zorder=3)

    # Shade comparison region
    #tau_limit = tau_ea[min(compare_n, len(tau_ea)) - 1] if compare_n <= len(tau_ea) else tau_ea[-1]
    #ax.axvspan(0, tau_limit, color="grey", alpha=0.10, label=f"Comparison region")

    ax.set_xlabel(r"Time lag $\tau$ [s]", fontsize=12)
    ax.set_ylabel(r"MSD [$\mu$m$^2$]", fontsize=12)
    if title:
        ax.set_title(title, fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(True, linestyle=":", alpha=0.5)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


#def compute_deviation(msd_ea, msd_ta, compare_n):
#    """Relative deviation metrics over the first *compare_n* lag points.

#    Returns (mean_rel_dev_pct, max_rel_dev_pct).
#    """
#    n = min(compare_n, len(msd_ea), len(msd_ta))
#    ea = msd_ea[:n]
#    ta = msd_ta[:n]
#    # Avoid division by zero — only consider points where EA-MSD > 0
#    valid = ea > 0
#    if not np.any(valid):
#        return float("nan"), float("nan")
#    rel = np.abs(ea[valid] - ta[valid]) / ea[valid]
#    return float(np.mean(rel)) * 100.0, float(np.max(rel)) * 100.0


def process_file(csv_file, compare_fraction):
    """Process one CSV: for each lag fraction, compute MSD curves and make plot."""
    stem = csv_file.stem
    print(f"  {csv_file.name}")

    trajectories = read_trajectories_from_csv(str(csv_file))
    if not trajectories:
        print("    no trajectories")
        return

    for lag_fraction in LAG_FRACTIONS:
        pct = int(lag_fraction * 100)
        tag = f"f{pct:03d}"

        # EA-MSD
        ea = calculate_ensemble_msd(trajectories, max_lag_fraction=lag_fraction)
        if ea.tau.size == 0:
            print(f"    {pct:2d}% — empty EA-MSD")
            continue

        # Ensemble-averaged TA-MSD (using the same lag fraction and global dt)
        tau_ta, msd_ta, sem_ta, per_track_matrix = compute_ensemble_tamsd(
            trajectories, max_lag_fraction=lag_fraction, global_dt=ea.dt,
        )

        # Trim both to common length
        n_common = min(ea.tau.size, tau_ta.size)
        tau_ea = ea.tau[:n_common]
        msd_ea = ea.msd[:n_common]
        sem_ea = ea.msd_sem[:n_common]
        tau_ta = tau_ta[:n_common]
        msd_ta = msd_ta[:n_common]
        sem_ta = sem_ta[:n_common]

        compare_n = max(2, int(compare_fraction * n_common))
        #mean_dev, max_dev = compute_deviation(msd_ea, msd_ta, compare_n)

        out_path = OUT_DIR / f"{stem}_ergodicity_{tag}.svg"
        plot_ergodicity(
            tau_ea, msd_ea, sem_ea,
            tau_ta, msd_ta, sem_ta,
            compare_n, out_path,
        )

        # Ergodicity breaking parameter
        eb, eb_ntracks = compute_eb_parameter(per_track_matrix[:, :n_common])
        eb_path = OUT_DIR / f"{stem}_eb_{tag}.svg"
        plot_eb_parameter(tau_ea, eb, eb_ntracks, eb_path)


def main():
    parser = argparse.ArgumentParser(description="Ergodicity check: EA-MSD vs TA-MSD (non-anomalous).")
    parser.add_argument("--compare-fraction", type=float, default=0.10,
                        help="Fraction of lag range used for deviation metric (default 0.10)")
    args = parser.parse_args()

    files = sorted(DATA_DIR.glob("*.csv"))
    if not files:
        print(f"No CSV files found in {DATA_DIR}")
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fracs_str = ", ".join(f"{int(f*100)}%" for f in LAG_FRACTIONS)
    print(f"Ergodicity check: {len(files)} files × lag fractions [{fracs_str}]")
    print(f"  compare_fraction = {args.compare_fraction}\n")

    for csv_file in files:
        process_file(csv_file, args.compare_fraction)

    print("\nDone.")


if __name__ == "__main__":
    main()
