"""
Ergodicity check — compare D from EA-MSD and ⟨TA-MSD⟩.

For an ergodic process the two MSD estimators must converge to the same
curve, yielding the same diffusion coefficient D.  This script extracts D
from both estimators (linear+offset fit on the first 10 % of the lag range)
for each 40 min-step experiment and produces a paired dot-plot comparison.

Note on the ergodicity-breaking (EB) parameter
-----------------------------------------------
For a simple Brownian system like beads in glycerol/water, the direct
comparison of D is the most informative ergodicity test.  A full EB(τ)
analysis is omitted here because at long time-lags the MSD is dominated
by drift (v²τ²), which inflates the EB variance without reflecting
genuine non-ergodic behaviour.

Output:
    Results/no_anomalous/ergodicity/
        D_comparison_paired.svg  — paired dot plot D_eaMSD vs D_⟨taMSD⟩

Usage:
    python check_ergodicity.py
"""

from __future__ import annotations

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
    compute_ensemble_drift,
    average_across_trajectories,
    build_tau_array,
    determine_maximum_lag_steps,
    estimate_global_time_step,
)
from msd_fitting import fit_msd_linear_offset

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "Data" / "31_10_no_anomalous"
OUT_DIR = SCRIPT_DIR / "Results" / "no_anomalous" / "ergodicity"

D_FIT_FRACTION = 0.10  # Fit D on first 10% of lag points


def compute_ensemble_tamsd(trajectories, max_lag_fraction=None, global_dt=None,
                           drift_corrected=False):
    """Compute the ensemble-averaged time-averaged MSD  ⟨TA-MSD⟩.

    For each trajectory we compute the full time-averaged MSD, then we
    average those curves across trajectories (same procedure as the EA-MSD
    averaging, but starting from TA-MSD per track instead of single initial
    displacements).

    Parameters
    ----------
    drift_corrected : bool
        If True, pass ``drift_corrected=True`` to each per-track TA-MSD
        computation (linear drift subtraction via regression).

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

    # Compute the common-mode (ensemble) drift once for this file, so each
    # per-track TAMSD subtracts the SAME drift instead of doing per-track OLS
    # (which biases the TAMSD by ~τ/T; Qian-Sheetz-Elson 1991, Vestergaard 2014).
    ens_drift = None
    if drift_corrected:
        ens_drift = compute_ensemble_drift(
            trajectories, global_dt=global_dt, min_tracks_per_step=5,
        )

    per_track_tamsd = []
    for track in trajectories.values():
        res = calculate_time_averaged_msd_per_track(
            track, max_lag_fraction=max_lag_fraction, dt_override=global_dt,
            drift_corrected=drift_corrected, ensemble_drift=ens_drift,
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


# ===================================================================
# D comparison: EA-MSD vs ⟨TA-MSD⟩ (paired dot plot)
# ===================================================================

def fit_D_from_msd(tau, msd, msd_sem, dt):
    """Fit linear_offset model on first 10 % of lag points.

    Returns (D, D_error) or (nan, nan) on failure.
    """
    n_max = len(tau)
    try:
        res = fit_msd_linear_offset(
            tau, msd, n_max=n_max, dt=dt,
            fit_fraction=D_FIT_FRACTION, msd_sigma=msd_sem,
        )
        return res.D, res.D_error
    except Exception:
        return float("nan"), float("nan")


def run_d_comparison(files):
    """For each file, fit D from eaMSD and ⟨taMSD⟩, then make paired dot plot."""
    labels = []
    d_ea_vals, d_ea_errs = [], []
    d_ta_vals, d_ta_errs = [], []

    # Compute MSD up to 50% so fit_msd_linear_offset can use 10% of that range
    msd_lag_fraction = 0.50

    for csv_file in files:
        stem = csv_file.stem
        short = stem.replace("_spots_40minstep", "")
        trajectories = read_trajectories_from_csv(str(csv_file))
        if not trajectories:
            continue

        # EA-MSD
        ea = calculate_ensemble_msd(trajectories, max_lag_fraction=msd_lag_fraction)
        if ea.tau.size == 0:
            continue

        # Ensemble-averaged TA-MSD
        tau_ta, msd_ta, sem_ta, _ = compute_ensemble_tamsd(
            trajectories, max_lag_fraction=msd_lag_fraction, global_dt=ea.dt,
        )

        # Trim to common length
        n = min(ea.tau.size, tau_ta.size)

        D_ea, D_ea_err = fit_D_from_msd(ea.tau[:n], ea.msd[:n], ea.msd_sem[:n], ea.dt)
        D_ta, D_ta_err = fit_D_from_msd(tau_ta[:n], msd_ta[:n], sem_ta[:n], ea.dt)

        labels.append(short)
        d_ea_vals.append(D_ea)
        d_ea_errs.append(D_ea_err)
        d_ta_vals.append(D_ta)
        d_ta_errs.append(D_ta_err)

    if not labels:
        print("  D comparison: no valid files")
        return

    d_ea = np.array(d_ea_vals)
    d_ta = np.array(d_ta_vals)
    d_ea_e = np.array(d_ea_errs)
    d_ta_e = np.array(d_ta_errs)

    plot_d_comparison(labels, d_ea, d_ea_e, d_ta, d_ta_e)
    print_d_statistics(labels, d_ea, d_ea_e, d_ta, d_ta_e)


def plot_d_comparison(labels, d_ea, d_ea_err, d_ta, d_ta_err):
    """Paired dot plot: two points per experiment connected by a line."""
    n = len(labels)
    x = np.arange(n)

    fig, ax = plt.subplots(figsize=(8, 5.5))

    # EA-MSD points
    ax.errorbar(x, d_ea, yerr=d_ea_err, fmt="o", color="C0",
                capsize=4, capthick=1.2, markersize=8, label="EA-MSD",
                zorder=3)
    # ⟨TA-MSD⟩ points
    ax.errorbar(x, d_ta, yerr=d_ta_err, fmt="o", color="C3",
                capsize=4, capthick=1.2, markersize=8,
                label=r"$\langle$TA-MSD$\rangle$", zorder=3)

    ax.set_xticks(x)
    ax.set_xticklabels([f"Experiment {i+1}" for i in range(n)], fontsize=10)
    ax.set_ylabel(r"$D$ [$\mu$m$^2$/s]", fontsize=12)
    #ax.set_title("D comparison: EA-MSD vs "
                 #r"$\langle$TA-MSD$\rangle$"
                 #f" (linear+offset, {int(D_FIT_FRACTION*100)}% lag)",
                 #fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(True, axis="y", linestyle=":", alpha=0.5)
    fig.tight_layout()

    out = OUT_DIR / "D_comparison_paired.svg"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"  -> {out.name}")


def print_d_statistics(labels, d_ea, d_ea_err, d_ta, d_ta_err):
    """Print detailed D comparison statistics to the terminal."""
    n = len(labels)

    print(f"\n  {'File':<15s}  {'D_eaMSD':>12s} {'± δD':>10s}  "
          f"{'D_⟨taMSD⟩':>12s} {'± δD':>10s}  {'ΔD/D̄ [%]':>10s}")
    print("  " + "─" * 80)

    rel_diffs = []
    for i in range(n):
        d_bar = 0.5 * (d_ea[i] + d_ta[i])
        rel = (d_ea[i] - d_ta[i]) / d_bar * 100 if d_bar > 0 else float("nan")
        rel_diffs.append(rel)
        print(f"  {labels[i]:<15s}  {d_ea[i]:12.4e} {d_ea_err[i]:10.2e}  "
              f"{d_ta[i]:12.4e} {d_ta_err[i]:10.2e}  {rel:+10.1f}")

    rel_arr = np.array(rel_diffs)
    abs_rel = np.abs(rel_arr[~np.isnan(rel_arr)])

    print("  " + "─" * 80)
    print(f"\n  Summary ({n} experiments, 40 min-step):")
    print(f"    ⟨D_eaMSD⟩  = {np.mean(d_ea):.4e} ± {np.std(d_ea, ddof=1):.2e} µm²/s")
    print(f"    ⟨D_taMSD⟩  = {np.mean(d_ta):.4e} ± {np.std(d_ta, ddof=1):.2e} µm²/s")
    print(f"    Mean |ΔD/D̄| = {np.mean(abs_rel):.1f} %")
    print(f"    Max  |ΔD/D̄| = {np.max(abs_rel):.1f} %  ({labels[np.argmax(abs_rel)]})")

    # Check consistency within combined fit errors
    consistent = 0
    for i in range(n):
        sigma_combined = np.sqrt(d_ea_err[i]**2 + d_ta_err[i]**2)
        if sigma_combined > 0 and abs(d_ea[i] - d_ta[i]) < 2 * sigma_combined:
            consistent += 1
    print(f"    Consistent within 2σ_fit: {consistent}/{n} experiments")
    print()


def main():
    files = sorted(DATA_DIR.glob("*_40minstep.csv"))
    if not files:
        print(f"No 40 min-step CSV files found in {DATA_DIR}")
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Ergodicity check: D comparison on {len(files)} × 40 min-step files")
    print(f"  Fit model: linear+offset on first {int(D_FIT_FRACTION*100)}% of lag range\n")

    run_d_comparison(files)

    print("Done.")


if __name__ == "__main__":
    main()
