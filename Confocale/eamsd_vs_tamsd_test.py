"""
Ergodicity / anomalous-diffusion diagnostic — EA-MSD vs ⟨TA-MSD⟩.

For each CSV file (per group: glic50, glic200) we fit the same
anomalous+offset model

    MSD(τ) = 4 D_α τ^α + c

on the SAME τ-range to two independent estimators of the MSD:

  * EA-MSD (drift-corrected, variance method)            → (α_EA, D_EA)
  * ⟨TA-MSD⟩ (ensemble average of per-track TA-MSDs,
    drift-corrected with the common-mode ensemble drift) → (α_TA, D_TA)

The test compares the two pairs:

  * Δα = α_TA − α_EA           (compatible with 0 if ergodic / Brownian)
  * ρ_D = D_TA / D_EA          (compatible with 1 if ergodic / Brownian)

Group-level statistics are reported as median + IQR over the files
of each group.  Sensitivity to the fit range is checked at 10 % and
25 % of the common lag range.

Output:
    Results/anomalous/eamsd_vs_tamsd/
        eamsd_vs_tamsd_<group>_f<pct>.svg   — log-log overlay (per file)
        alpha_paired_<group>_f<pct>.svg     — α paired dot plot
        D_paired_<group>_f<pct>.svg         — D paired dot plot
    Docu/eamsd_vs_tamsd_per_file.csv        — one row per file × fraction
    Docu/eamsd_vs_tamsd_group_summary.csv   — median + IQR per group × fraction

Usage:
    python eamsd_vs_tamsd_test.py                  # both groups
    python eamsd_vs_tamsd_test.py --group glic50
    python eamsd_vs_tamsd_test.py --group glic200
"""

from __future__ import annotations

import argparse
import re
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", message="Degrees of freedom <= 0 for slice")

from check_ergodicity import compute_ensemble_tamsd
from data_reader import read_trajectories_from_csv, estimate_global_time_step
from msd_analyzer import calculate_ensemble_msd
from msd_fitting import fit_msd_anomalous_offset


# ── Configuration ──────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "Data" / "14_11_anomalous"
OUT_DIR = SCRIPT_DIR / "Results" / "anomalous" / "eamsd_vs_tamsd"
DOC_DIR = SCRIPT_DIR / "Docu"

GLIC50_SUBDIR = "PEG_15_Glicerolo_50_H2O_35"
GLIC200_SUBDIR = "PEG_18_Glicerolo_40_H2O_42"

# Compute MSD curves up to this fraction of the longest trajectory; the
# fit is then restricted to FIT_FRACTIONS of this curve.
MSD_LAG_FRACTION = 0.50
FIT_FRACTIONS = [0.10, 0.25]
DRIFT_CORRECTED = True
# Show median curves only where at least this many files contribute.
# This avoids edge artifacts at the first lag when only a few files are valid.
MIN_VALID_FILES_FOR_MEDIAN = 4
# Show only an early-lag window around the fit region to improve readability.
PLOT_WINDOW_FACTOR = 1.6


# ── Plot style helpers ────────────────────────────────────────────────

def _short_file_label(stem: str, idx: int) -> str:
    """Compact label for x-axis to avoid long filenames in paired plots."""
    m = re.search(r"series(\d+)", stem, flags=re.IGNORECASE)
    if m:
        return f"S{int(m.group(1))}"
    m = re.search(r"_serie(\d+)", stem, flags=re.IGNORECASE)
    if m:
        return f"S{int(m.group(1))}"
    return f"F{idx + 1}"


def _apply_clean_axes_style(ax, log_grid=False):
    """Consistent, report-friendly visual style used across all plots."""
    ax.set_facecolor("#f7f9fc")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#4b5563")
    ax.spines["bottom"].set_color("#4b5563")
    ax.tick_params(axis="both", labelsize=10, colors="#374151")
    if log_grid:
        ax.grid(True, which="major", linestyle=":", linewidth=0.9,
                color="#cbd5e1", alpha=0.9)
        ax.grid(True, which="minor", linestyle=":", linewidth=0.6,
                color="#e2e8f0", alpha=0.9)
    else:
        ax.grid(True, axis="y", linestyle=":", linewidth=0.9,
                color="#cbd5e1", alpha=0.9)


# ── File grouping ──────────────────────────────────────────────────────

def group_csv_files(data_dir: Path):
    glic50_dir = data_dir / GLIC50_SUBDIR
    glic200_dir = data_dir / GLIC200_SUBDIR
    return {
        "glic50": sorted(glic50_dir.glob("*.csv")) if glic50_dir.is_dir() else [],
        "glic200": sorted(glic200_dir.glob("*.csv")) if glic200_dir.is_dir() else [],
    }


# ── Per-file analysis ──────────────────────────────────────────────────

def analyse_file(csv_path: Path):
    """Compute EA-MSD and ⟨TA-MSD⟩ for one file on the same τ-range,
    return a dict with both curves and the per-fraction fits.
    """
    trajectories = read_trajectories_from_csv(str(csv_path))
    if not trajectories:
        return None

    global_dt = estimate_global_time_step(trajectories)

    # EA-MSD (drift-corrected)
    ea = calculate_ensemble_msd(
        trajectories, max_lag_fraction=MSD_LAG_FRACTION,
        global_dt=global_dt, drift_corrected=DRIFT_CORRECTED,
    )
    if ea.tau.size < 4:
        return None

    # ⟨TA-MSD⟩ (drift-corrected with common-mode ensemble drift)
    tau_ta, msd_ta, sem_ta, _ = compute_ensemble_tamsd(
        trajectories, max_lag_fraction=MSD_LAG_FRACTION,
        global_dt=global_dt, drift_corrected=DRIFT_CORRECTED,
    )

    # Common τ-range
    n = min(ea.tau.size, tau_ta.size)
    tau_common = ea.tau[:n]
    ea_msd, ea_sem = ea.msd[:n], ea.msd_sem[:n]
    ta_msd, ta_sem = msd_ta[:n], sem_ta[:n]

    fits = {}
    for frac in FIT_FRACTIONS:
        try:
            fit_ea = fit_msd_anomalous_offset(
                tau_common, ea_msd, n, global_dt,
                fit_fraction=frac, msd_sigma=ea_sem,
            )
            alpha_ea, D_ea = fit_ea.alpha, fit_ea.D_alpha
        except (ValueError, RuntimeError):
            alpha_ea, D_ea = np.nan, np.nan

        try:
            fit_ta = fit_msd_anomalous_offset(
                tau_common, ta_msd, n, global_dt,
                fit_fraction=frac, msd_sigma=ta_sem,
            )
            alpha_ta, D_ta = fit_ta.alpha, fit_ta.D_alpha
        except (ValueError, RuntimeError):
            alpha_ta, D_ta = np.nan, np.nan

        fits[frac] = dict(
            alpha_ea=alpha_ea, D_ea=D_ea,
            alpha_ta=alpha_ta, D_ta=D_ta,
        )

    return dict(
        stem=csv_path.stem,
        n_tracks=len(trajectories),
        dt=global_dt,
        tau=tau_common,
        ea_msd=ea_msd, ea_sem=ea_sem,
        ta_msd=ta_msd, ta_sem=ta_sem,
        fits=fits,
    )


# ── Plotting ──────────────────────────────────────────────────────────

def plot_loglog_overlay(group_name, results, frac, out_path):
    """Per-file log-log curves of EA-MSD and ⟨TA-MSD⟩ on the same axes,
    with the fit range highlighted."""
    fig, ax = plt.subplots(figsize=(8.2, 6.0))
    fig.patch.set_facecolor("white")
    _apply_clean_axes_style(ax, log_grid=False)

    # Fit-range shading: the largest τ used at this frac, taken from the
    # longest curve (gives a representative band).
    long_idx = int(np.argmax([len(r["tau"]) for r in results]))
    longest_tau = results[long_idx]["tau"]
    n_fit = max(int(np.floor(frac * len(longest_tau))), 4)
    tau_fit_max = longest_tau[n_fit - 1]
    tau_fit_min = longest_tau[0]

    # Do not show the entire tail: keep a zoomed window near the fit region.
    tau_plot_max = tau_fit_max * PLOT_WINDOW_FACTOR

    for r in results:
        m = (r["tau"] >= tau_fit_min) & (r["tau"] <= tau_plot_max)
        ax.plot(r["tau"][m], r["ea_msd"][m], color="C0", alpha=0.30, lw=1.1)
        ax.plot(r["tau"][m], r["ta_msd"][m], color="C3", alpha=0.30, lw=1.1)

    # Median across files on a common log-spaced grid
    tau_all = np.unique(np.concatenate([
        r["tau"][(r["tau"] >= tau_fit_min) & (r["tau"] <= tau_plot_max)]
        for r in results
    ]))
    tau_all = tau_all[tau_all > 0]
    ea_stack, ta_stack = [], []
    for r in results:
        ea_stack.append(np.interp(tau_all, r["tau"], r["ea_msd"],
                                  left=np.nan, right=np.nan))
        ta_stack.append(np.interp(tau_all, r["tau"], r["ta_msd"],
                                  left=np.nan, right=np.nan))
    ea_stack = np.vstack(ea_stack)
    ta_stack = np.vstack(ta_stack)

    # Keep only lags where enough files contribute to BOTH estimators.
    min_valid = min(MIN_VALID_FILES_FOR_MEDIAN, len(results))
    n_valid_both = np.sum(np.isfinite(ea_stack) & np.isfinite(ta_stack), axis=0)
    valid_cols = n_valid_both >= min_valid

    tau_med = tau_all[valid_cols]
    ea_med = np.nanmedian(ea_stack[:, valid_cols], axis=0)
    ta_med = np.nanmedian(ta_stack[:, valid_cols], axis=0)

    ax.plot(tau_med, ea_med, color="C0", lw=2.7, label="EA-MSD (median)")
    ax.plot(tau_med, ta_med, color="C3", lw=2.7, label=r"$\langle$TA-MSD$\rangle$ (median)")

    ax.axvspan(tau_fit_min, tau_fit_max, color="#94a3b8", alpha=0.18,
               label=f"Fit range ({int(frac*100)}%)")

    # Show fitted alpha values directly on the log-log panel for reporting.
    alpha_ea_vals = np.array([r["fits"][frac]["alpha_ea"] for r in results],
                             dtype=float)
    alpha_ta_vals = np.array([r["fits"][frac]["alpha_ta"] for r in results],
                             dtype=float)
    med_alpha_ea = np.nanmedian(alpha_ea_vals)
    med_alpha_ta = np.nanmedian(alpha_ta_vals)
    if np.isfinite(med_alpha_ea) and np.isfinite(med_alpha_ta):
        alpha_text = (
            rf"$\alpha_{{EA}}$ = {med_alpha_ea:.3f}" + "\n" +
            rf"$\alpha_{{\langle TA \rangle}}$  = {med_alpha_ta:.3f}"
        )
        ax.text(
            0.02, 0.87, alpha_text,
            transform=ax.transAxes,
            va="top", ha="left", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="#cbd5e1", alpha=0.96),
        )

    ax.set_xlabel(r"Lag time $\tau$ [s]", fontsize=12)
    ax.set_ylabel(r"MSD [$\mu$m$^2$]", fontsize=12)
    # Match the linear scale used in linear_offset_fits-style plots.
    ax.set_xscale("linear")
    ax.set_yscale("linear")
    ax.set_xlim(tau_fit_min * 0.95, tau_plot_max)

    y_pool = np.concatenate([ea_med[np.isfinite(ea_med)], ta_med[np.isfinite(ta_med)]])
    if y_pool.size:
        y_min = np.min(y_pool)
        y_max = np.max(y_pool)
        if np.isfinite(y_min) and np.isfinite(y_max) and y_max > y_min:
            pad = 0.12 * (y_max - y_min)
            ax.set_ylim(max(0.0, y_min - pad), y_max + pad)
    ax.legend(loc="lower right", fontsize=10, framealpha=0.95,
              facecolor="white", edgecolor="#cbd5e1")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_paired(group_name, labels, vals_ea, vals_ta, ylabel, ref_value,
                ref_label, out_path):
    """Paired dot plot: EA vs TA, one couple per file connected by a line."""
    n = len(labels)
    x = np.arange(n)
    fig, ax = plt.subplots(figsize=(max(8, 0.58 * n + 4.3), 5.4))
    fig.patch.set_facecolor("white")
    _apply_clean_axes_style(ax, log_grid=False)

    # Alternating vertical bands help readability on dense paired plots.
    for i in range(n):
        if i % 2 == 0:
            ax.axvspan(i - 0.5, i + 0.5, color="white", alpha=0.38, zorder=0)

    for xi, vea, vta in zip(x, vals_ea, vals_ta):
        if np.isfinite(vea) and np.isfinite(vta):
            ax.plot([xi - 0.14, xi + 0.14], [vea, vta], color="#6b7280",
                    lw=1.1, alpha=0.75, zorder=2)
    ax.plot(x - 0.14, vals_ea, "o", color="C0", markersize=8.5,
            markeredgecolor="white", markeredgewidth=0.8,
            label="EA-MSD", zorder=3)
    ax.plot(x + 0.14, vals_ta, "s", color="C3", markersize=8.2,
            markeredgecolor="white", markeredgewidth=0.8,
            label=r"$\langle$TA-MSD$\rangle$", zorder=3)

    if ref_value is not None:
        ax.axhline(ref_value, color="#111827", lw=1.3, linestyle="--",
                   label=ref_label)

    # Median bands
    med_ea = np.nanmedian(vals_ea)
    med_ta = np.nanmedian(vals_ta)
    ax.axhline(med_ea, color="C0", lw=1.1, linestyle=":",
               label=f"Median EA = {med_ea:.3g}")
    ax.axhline(med_ta, color="C3", lw=1.1, linestyle=":",
               label=f"Median TA = {med_ta:.3g}")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=0, ha="center",
                       fontsize=(9 if n <= 10 else 8.5))
    ax.set_xlabel("File", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.margins(x=0.02)
    ax.legend(loc="best", fontsize=9, framealpha=0.95,
              facecolor="white", edgecolor="#cbd5e1")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ── Statistics helpers ────────────────────────────────────────────────

def median_iqr(values):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.nan, np.nan, np.nan, 0
    q1, med, q3 = np.percentile(arr, [25, 50, 75])
    return float(med), float(q1), float(q3), int(arr.size)


# ── Main per-group routine ────────────────────────────────────────────

def run_group(group_name: str, csv_files):
    print(f"\n{'='*60}\n[{group_name}] {len(csv_files)} file(s)\n{'='*60}")
    results = []
    for csv_path in csv_files:
        print(f"  -> {csv_path.name}")
        r = analyse_file(csv_path)
        if r is None:
            print("     skipped (no/short trajectories)")
            continue
        results.append(r)

    if not results:
        print(f"[{group_name}] no usable files.")
        return [], []

    # Per-file rows
    per_file_rows = []
    for r in results:
        for frac, f in r["fits"].items():
            per_file_rows.append(dict(
                group=group_name,
                file=r["stem"],
                n_tracks=r["n_tracks"],
                dt=r["dt"],
                fit_fraction=frac,
                alpha_EA=f["alpha_ea"], D_EA=f["D_ea"],
                alpha_TA=f["alpha_ta"], D_TA=f["D_ta"],
                delta_alpha=f["alpha_ta"] - f["alpha_ea"],
                rho_D=(f["D_ta"] / f["D_ea"]
                       if (f["D_ea"] not in (0, np.nan) and np.isfinite(f["D_ea"]))
                       else np.nan),
            ))

    # Group summary rows + plots per fraction
    group_rows = []
    for frac in FIT_FRACTIONS:
        sub = [row for row in per_file_rows if row["fit_fraction"] == frac]
        if not sub:
            continue

        labels = [_short_file_label(row["file"], i) for i, row in enumerate(sub)]
        a_ea = np.array([row["alpha_EA"] for row in sub], dtype=float)
        a_ta = np.array([row["alpha_TA"] for row in sub], dtype=float)
        d_ea = np.array([row["D_EA"] for row in sub], dtype=float)
        d_ta = np.array([row["D_TA"] for row in sub], dtype=float)
        dalpha = a_ta - a_ea
        rho = np.where(np.isfinite(d_ea) & (d_ea != 0), d_ta / d_ea, np.nan)

        for name, arr in [("alpha_EA", a_ea), ("alpha_TA", a_ta),
                          ("D_EA", d_ea), ("D_TA", d_ta),
                          ("delta_alpha", dalpha), ("rho_D", rho)]:
            med, q1, q3, n = median_iqr(arr)
            group_rows.append(dict(
                group=group_name, fit_fraction=frac, quantity=name,
                median=med, q1=q1, q3=q3, iqr=q3 - q1, n=n,
            ))

        # Plots
        pct = int(frac * 100)
        plot_loglog_overlay(group_name, results, frac,
                            OUT_DIR / f"eamsd_vs_tamsd_{group_name}_f{pct:03d}.svg")
        plot_paired(group_name, labels, a_ea, a_ta, r"$\alpha$",
                    1.0, r"$\alpha = 1$ (Brownian)",
                    OUT_DIR / f"alpha_paired_{group_name}_f{pct:03d}.svg")
        plot_paired(group_name, labels, d_ea, d_ta,
                    r"$D_\alpha$ [$\mu$m$^2$/s$^\alpha$]",
                    None, None,
                    OUT_DIR / f"D_paired_{group_name}_f{pct:03d}.svg")

        med_da, q1_da, q3_da, _ = median_iqr(dalpha)
        med_rho, q1_rho, q3_rho, _ = median_iqr(rho)
        print(f"  [{group_name} | fit {pct}%]"
              f"  median delta_alpha = {med_da:+.3f} (IQR [{q1_da:+.3f},{q3_da:+.3f}])"
              f"   median rho_D = {med_rho:.3f} (IQR [{q1_rho:.3f},{q3_rho:.3f}])")

    return per_file_rows, group_rows


# ── Entry point ───────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", choices=["glic50", "glic200", "all"],
                        default="all")
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DOC_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Drift correction enabled: {DRIFT_CORRECTED}")

    files_by_group = group_csv_files(DATA_DIR)
    selected = (["glic50", "glic200"] if args.group == "all" else [args.group])

    all_per_file = []
    all_group = []
    for g in selected:
        per_file, group = run_group(g, files_by_group[g])
        all_per_file.extend(per_file)
        all_group.extend(group)

    if all_per_file:
        per_file_path = DOC_DIR / "eamsd_vs_tamsd_per_file.csv"
        pd.DataFrame(all_per_file).to_csv(per_file_path, index=False)
        print(f"\nWrote {per_file_path}")
    if all_group:
        group_path = DOC_DIR / "eamsd_vs_tamsd_group_summary.csv"
        pd.DataFrame(all_group).to_csv(group_path, index=False)
        print(f"Wrote {group_path}")


if __name__ == "__main__":
    main()
