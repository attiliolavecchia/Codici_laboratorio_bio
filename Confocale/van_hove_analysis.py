"""
Complementary van Hove analysis for anomalous trajectories.

This script computes the self-part van Hove distribution from 2D trajectories
and plots the 1D displacement PDF Gs(dx, tau) at representative lag times.
It also reports the 2D non-Gaussian parameter:

    beta(tau) = <r^4> / (2 <r^2>^2) - 1

for each lag. These diagnostics are useful to characterize dynamic
heterogeneity and non-Gaussian transport, as a complement to EA-MSD vs <TA-MSD>.

Output:
    Results/anomalous/van_hove/
        van_hove_dx_<group>_lag<steps>_tau<seconds>s.svg
        beta_vs_tau_<group>.svg
    Docu/van_hove_summary.csv

Usage:
    python van_hove_analysis.py
    python van_hove_analysis.py --group glic50
    python van_hove_analysis.py --group glic200
    python van_hove_analysis.py --group all --no-drift-correction
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from data_reader import read_trajectories_from_csv, estimate_global_time_step
from msd_analyzer import compute_ensemble_drift, determine_maximum_lag_steps


SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "Data" / "14_11_anomalous"
OUT_DIR = SCRIPT_DIR / "Results" / "anomalous" / "van_hove"
DOC_DIR = SCRIPT_DIR / "Docu"

GLIC50_SUBDIR = "PEG_15_Glicerolo_50_H2O_35"
GLIC200_SUBDIR = "PEG_18_Glicerolo_40_H2O_42"

MAX_LAG_FRACTION = 0.50
LAG_STEP_FRACTIONS = [0.05, 0.15, 0.30]
BETA_N_TAU_POINTS = 12
N_BINS = 90


def _apply_clean_axes_style(ax):
    ax.set_facecolor("#f7f9fc")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#4b5563")
    ax.spines["bottom"].set_color("#4b5563")
    ax.tick_params(axis="both", labelsize=10, colors="#374151")
    ax.grid(True, linestyle=":", linewidth=0.8, color="#cbd5e1", alpha=0.9)


def group_csv_files(data_dir: Path):
    glic50_dir = data_dir / GLIC50_SUBDIR
    glic200_dir = data_dir / GLIC200_SUBDIR
    return {
        "glic50": sorted(glic50_dir.glob("*.csv")) if glic50_dir.is_dir() else [],
        "glic200": sorted(glic200_dir.glob("*.csv")) if glic200_dir.is_dir() else [],
    }


def select_lag_steps(n_max: int, fractions):
    steps = []
    for frac in fractions:
        k = int(round(float(frac) * n_max))
        k = max(1, min(n_max, k))
        steps.append(k)
    return sorted(set(steps))


def select_evenly_spaced_lag_steps(n_max: int, n_points: int):
    if n_max < 1:
        return []
    if n_points <= 1:
        return [1]
    grid = np.linspace(1, int(n_max), int(n_points))
    steps = np.rint(grid).astype(int)
    steps = np.clip(steps, 1, int(n_max))
    return sorted(set(int(s) for s in steps))


def collect_displacements_for_file(trajectories, lag_steps, global_dt, drift_corrected):
    per_lag = {
        lag: {"dx": [], "dy": [], "tau_s": lag * float(global_dt)}
        for lag in lag_steps
    }

    ens_drift = None
    t_min = None
    if drift_corrected:
        ens_drift = compute_ensemble_drift(
            trajectories, global_dt=global_dt, min_tracks_per_step=5,
        )
        t_min = float(ens_drift.time[0]) if ens_drift.time.size else 0.0

    for tr in trajectories.values():
        x = np.asarray(tr.x, dtype=float)
        y = np.asarray(tr.y, dtype=float)

        if drift_corrected and ens_drift is not None and ens_drift.time.size:
            t = np.asarray(tr.time, dtype=float)
            frames = np.rint((t - t_min) / float(global_dt)).astype(int)
            frames = np.clip(frames, 0, ens_drift.time.size - 1)
            x = x - ens_drift.Rx[frames]
            y = y - ens_drift.Ry[frames]

        n = x.size
        for lag in lag_steps:
            if n <= lag:
                continue
            dx = x[lag:] - x[:-lag]
            dy = y[lag:] - y[:-lag]
            if dx.size == 0:
                continue
            per_lag[lag]["dx"].append(dx)
            per_lag[lag]["dy"].append(dy)

    out = {}
    for lag in lag_steps:
        if per_lag[lag]["dx"]:
            dx = np.concatenate(per_lag[lag]["dx"])
            dy = np.concatenate(per_lag[lag]["dy"])
        else:
            dx = np.asarray([], dtype=float)
            dy = np.asarray([], dtype=float)
        out[lag] = {
            "dx": dx,
            "dy": dy,
            "tau_s": per_lag[lag]["tau_s"],
        }
    return out


def _format_tau_for_filename(tau_s: float) -> str:
    # Compact, filesystem-safe tau string (e.g., 1p23e-02)
    s = f"{float(tau_s):.3e}"
    return s.replace("+", "").replace("-", "m").replace(".", "p")


def plot_van_hove_dx_separate(group_name: str, plot_items, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)

    for item in sorted(plot_items, key=lambda x: x["tau_s"]):
        lag_steps = int(item["lag_steps"])
        dx = item["dx"]
        tau_s = float(item["tau_s"])
        beta = float(item["beta_2D"])

        fig, ax = plt.subplots(figsize=(6.2, 4.8))
        fig.patch.set_facecolor("white")
        _apply_clean_axes_style(ax)

        q = np.nanpercentile(np.abs(dx), 99.5)
        if (not np.isfinite(q)) or q <= 0:
            q = np.nanmax(np.abs(dx)) if dx.size else 1.0
        q = float(max(q, 1e-8))

        bins = np.linspace(-q, q, N_BINS + 1)
        hist, edges = np.histogram(dx, bins=bins, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])

        ax.plot(centers, hist, color="C0", lw=2.0, label=r"$G_s(\Delta x,\tau)$")

        mu = float(np.mean(dx))
        sigma = float(np.std(dx, ddof=0))
        if np.isfinite(sigma) and sigma > 0:
            gauss = (1.0 / (np.sqrt(2.0 * np.pi) * sigma)) * np.exp(
                -0.5 * ((centers - mu) / sigma) ** 2
            )
            ax.plot(centers, gauss, color="C3", lw=1.8, linestyle="--",
                    label="Gaussian match")

        # No figure title by user request.
        ax.set_xlabel(r"$\Delta x$ [$\mu$m]", fontsize=11)
        ax.set_ylabel("PDF", fontsize=11)

        text = (
            rf"$N_{{\mathrm{{jumps}}}} = {dx.size}$"
            + "\n"
            + rf"$\beta = {beta:.3g}$"
            + "\n"
            + rf"$\tau = {tau_s:.3g}\,\mathrm{{s}}$"
        )
        ax.text(0.03, 0.97, text, transform=ax.transAxes,
                va="top", ha="left", fontsize=9.5,
                bbox=dict(boxstyle="round,pad=0.28", facecolor="white",
                          edgecolor="#cbd5e1", alpha=0.96))
        ax.legend(fontsize=9)

        tau_tag = _format_tau_for_filename(tau_s)
        out_path = out_dir / f"van_hove_dx_{group_name}_lag{lag_steps:03d}_tau{tau_tag}s.svg"
        plt.tight_layout()
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)


def plot_beta_vs_tau(group_name: str, summary_rows, out_path: Path):
    rows = sorted(summary_rows, key=lambda r: r["tau_s"])
    tau = np.asarray([r["tau_s"] for r in rows], dtype=float)
    beta = np.asarray([r["beta_2D"] for r in rows], dtype=float)

    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    fig.patch.set_facecolor("white")
    _apply_clean_axes_style(ax)

    ax.plot(tau, beta, "o-", color="C2", lw=1.8, ms=7, label=r"$\beta(\tau)$")
    ax.axhline(0.0, color="#111827", lw=1.2, linestyle="--", label="Gaussian baseline")

    ax.set_xlabel(r"$\tau$ [s]", fontsize=11)
    ax.set_ylabel(r"$\beta(\tau)$", fontsize=11)
    #ax.set_title(f"2D non-Gaussian parameter - {group_name}", fontsize=12)
    ax.legend(fontsize=9)

    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def run_group(group_name: str, csv_files, drift_corrected: bool):
    print(f"\n{'='*62}\n[{group_name}] {len(csv_files)} file(s)\n{'='*62}")

    file_entries = []
    for csv_path in csv_files:
        trajectories = read_trajectories_from_csv(str(csv_path))
        if not trajectories:
            print(f"  -> {csv_path.name}: skipped (no trajectories)")
            continue

        max_points = max(t.n_points for t in trajectories.values())
        if max_points < 3:
            print(f"  -> {csv_path.name}: skipped (too short)")
            continue

        dt = estimate_global_time_step(trajectories)
        n_max = determine_maximum_lag_steps(max_points, MAX_LAG_FRACTION)
        if n_max < 1:
            print(f"  -> {csv_path.name}: skipped (n_max < 1)")
            continue

        file_entries.append({
            "name": csv_path.name,
            "trajectories": trajectories,
            "dt": float(dt),
            "n_max": int(n_max),
        })

    if not file_entries:
        print(f"[{group_name}] no usable files.")
        return []

    common_n_max = min(entry["n_max"] for entry in file_entries)
    lag_steps = select_lag_steps(common_n_max, LAG_STEP_FRACTIONS)
    beta_lag_steps = select_evenly_spaced_lag_steps(common_n_max, BETA_N_TAU_POINTS)
    all_lag_steps = sorted(set(lag_steps) | set(beta_lag_steps))
    print(
        f"  common n_max = {common_n_max}; "
        f"van Hove lag steps = {lag_steps}; beta points = {len(beta_lag_steps)}"
    )

    accum = {
        lag: {"dx": [], "dy": [], "tau_s": []}
        for lag in all_lag_steps
    }

    for entry in file_entries:
        print(f"  -> {entry['name']}")
        displ = collect_displacements_for_file(
            entry["trajectories"], all_lag_steps, entry["dt"], drift_corrected,
        )
        for lag in all_lag_steps:
            dx = displ[lag]["dx"]
            dy = displ[lag]["dy"]
            if dx.size == 0:
                continue
            accum[lag]["dx"].append(dx)
            accum[lag]["dy"].append(dy)
            accum[lag]["tau_s"].append(displ[lag]["tau_s"])

    rows_by_lag = {}
    for lag in all_lag_steps:
        if not accum[lag]["dx"]:
            continue

        dx = np.concatenate(accum[lag]["dx"])
        dy = np.concatenate(accum[lag]["dy"])
        tau_s = float(np.median(np.asarray(accum[lag]["tau_s"], dtype=float)))

        r2 = dx * dx + dy * dy
        r4 = r2 * r2
        mean_r2 = float(np.mean(r2))
        mean_r4 = float(np.mean(r4))
        beta = (mean_r4 / (2.0 * mean_r2 * mean_r2) - 1.0) if mean_r2 > 0 else np.nan

        var_dx = float(np.var(dx, ddof=0))
        if var_dx > 0:
            m4x = float(np.mean((dx - np.mean(dx)) ** 4))
            excess_kurt_x = m4x / (var_dx * var_dx) - 3.0
        else:
            excess_kurt_x = np.nan

        rows_by_lag[int(lag)] = {
            "group": group_name,
            "lag_steps": int(lag),
            "tau_s": tau_s,
            "n_jumps": int(dx.size),
            "mean_r2": mean_r2,
            "mean_r4": mean_r4,
            "beta_2D": float(beta),
            "excess_kurtosis_dx": float(excess_kurt_x),
            "drift_corrected": bool(drift_corrected),
            "dx": dx,
        }

    summary_rows = [rows_by_lag[lag] for lag in lag_steps if lag in rows_by_lag]
    beta_rows = [rows_by_lag[lag] for lag in beta_lag_steps if lag in rows_by_lag]
    plot_items = [
        {
            "lag_steps": int(row["lag_steps"]),
            "tau_s": float(row["tau_s"]),
            "dx": row["dx"],
            "beta_2D": float(row["beta_2D"]),
        }
        for row in summary_rows
    ]

    if not summary_rows:
        print(f"[{group_name}] no displacement samples at selected lags.")
        return []

    plot_van_hove_dx_separate(
        group_name,
        plot_items,
        OUT_DIR,
    )
    plot_beta_vs_tau(
        group_name,
        beta_rows,
        OUT_DIR / f"beta_vs_tau_{group_name}.svg",
    )

    for row in sorted(summary_rows, key=lambda x: x["tau_s"]):
        print(
            f"  [tau={row['tau_s']:.3g}s] "
            f"n_jumps={row['n_jumps']}, beta={row['beta_2D']:+.3f}, "
            f"kurt_dx={row['excess_kurtosis_dx']:+.3f}"
        )

    return summary_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", choices=["glic50", "glic200", "all"], default="all")
    parser.add_argument("--no-drift-correction", action="store_true",
                        help="Use raw displacements without ensemble drift subtraction.")
    args = parser.parse_args()

    #drift_corrected = not args.no_drift_correction
    drift_corrected = True # Force drift correction for better results, as per user request.

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DOC_DIR.mkdir(parents=True, exist_ok=True)

    files_by_group = group_csv_files(DATA_DIR)
    selected = ["glic50", "glic200"] if args.group == "all" else [args.group]

    all_rows = []
    for group_name in selected:
        rows = run_group(group_name, files_by_group[group_name], drift_corrected)
        all_rows.extend(rows)

    if all_rows:
        out_csv = DOC_DIR / "van_hove_summary.csv"
        pd.DataFrame(all_rows).to_csv(out_csv, index=False)
        print(f"\nWrote {out_csv}")


if __name__ == "__main__":
    main()
