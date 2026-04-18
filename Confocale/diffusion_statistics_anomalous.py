"""
Rigorous statistical analysis of anomalous diffusion parameters D_α and α.

Mirrors ``diffusion_statistics.py`` for the anomalous dataset, using
anomalous MSD models instead of linear/nonlinear.

Workflow
--------
1. Per-track taMSD fitting  → individual (D_α, α) values (pooled across all CSV files)
2. Per-file  eaMSD fitting  → one (D_α, α) per CSV file
3. Per-file  ⟨taMSD⟩ fitting → one (D_α, α) per CSV file
4. Histograms of D_α and α for each estimator
5. Summary comparison and statistics

Models (anomalous equivalents of linear / linear_offset / nonlinear):
    anomalous:        MSD(τ) = 4D_α τ^α
    anomalous_offset: MSD(τ) = 4D_α τ^α + c   (c = 4σ², localization error)
    anomalous_drift:  MSD(τ) = 4D_α τ^α + v²τ²

Fit fractions: 10 % and 25 %
Dataset: anomalous only (Data/14_11_anomalous)
    ─ glicerolo50:   240nm beads in 50 % glycerol
    ─ glicerolo200:  240nm beads in concentrated glycerol (200×)

Output (SVG for paper, CSV for data):
    Results/anomalous/diffusion_statistics/     — histograms
    Results/anomalous/linear_fits/              — anomalous + anomalous_offset fit plots
    Results/anomalous/nonlinear_fits/           — anomalous_drift fit plots
    Docu/diffusion_statistics_anomalous_*.csv   — tabular results

Usage:
    python diffusion_statistics_anomalous.py                  # full analysis
    python diffusion_statistics_anomalous.py --group glic50   # only glicerolo50
    python diffusion_statistics_anomalous.py --group glic200  # only glicerolo200
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from check_ergodicity import compute_ensemble_tamsd
from data_reader import read_trajectories_from_csv, estimate_global_time_step
from msd_analyzer import calculate_ensemble_msd, calculate_time_averaged_msd_per_track
from msd_fitting import (
    fit_msd_anomalous_offset,
    anomalous_offset_msd_model,
    _fit_anomalous_at_fraction,
    anomalous_msd_model,
    _fit_anomalous_drift_at_fraction,
    anomalous_drift_msd_model,
    analyze_velocities,
)

# ── Configuration ──────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "Data" / "14_11_anomalous"
STATS_DIR = SCRIPT_DIR / "Results" / "anomalous" / "diffusion_statistics"
LINEAR_DIR = SCRIPT_DIR / "Results" / "anomalous" / "linear_fits"
LINEAR_OFFSET_DIR = SCRIPT_DIR / "Results" / "anomalous" / "linear_offset_fits"
NONLINEAR_DIR = SCRIPT_DIR / "Results" / "anomalous" / "nonlinear_fits"
DOC_DIR = SCRIPT_DIR / "Docu"

MIN_TRACK_POINTS = 30
FIT_FRACTIONS = [0.10, 0.25]

# File grouping
GLIC50_PATTERN = "240nm_glicerolo50"
GLIC200_PATTERN = "240nm_glicerolo200"


# ── Plotting helpers ───────────────────────────────────────────────────

def _save_fit_plot(tau_fit, msd_fit, msd_predicted, textstr, fit_label,
                   output_path, msd_sigma=None, title=None, data_color="C0"):
    """Save a single MSD fit plot (data + curve + annotation box)."""
    fig, ax = plt.subplots(figsize=(8, 6))
    if msd_sigma is not None:
        ax.errorbar(tau_fit, msd_fit, yerr=msd_sigma, fmt="o", color=data_color,
                    markersize=8, alpha=0.7, capsize=4, capthick=1.2,
                    elinewidth=1.2, label="MSD Data", zorder=2)
    else:
        ax.plot(tau_fit, msd_fit, "o", color=data_color, markersize=8, alpha=0.7,
                label="MSD Data", zorder=2)
    ax.plot(tau_fit, msd_predicted, "-", color="C3", linewidth=2.5,
            label=fit_label, zorder=3)
    ax.set_xlabel(r"Time Lag $\tau$ [s]", fontsize=12)
    ax.set_ylabel(r"MSD [$\mu$m$^2$]", fontsize=12)
    if title:
        ax.set_title(title, fontsize=13)
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(loc="upper left", fontsize=10, framealpha=0.9)
    props = dict(boxstyle="round", facecolor="white", alpha=0.95,
                 edgecolor="black", linewidth=1.2)
    ax.text(0.98, 0.97, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment="top", horizontalalignment="right", bbox=props)
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_histogram(values, xlabel, output_path,
                    mean_val=None, median_val=None, title=None,
                    ref_line=None, ref_label=None):
    """Plot histogram with optional reference, mean, and median lines."""
    fig, ax = plt.subplots(figsize=(8, 5.5))
    n_bins = max(10, int(np.sqrt(len(values))))
    ax.hist(values, bins=n_bins, color="steelblue", edgecolor="white",
            alpha=0.8, density=True, label="Data")
    if ref_line is not None and ref_label is not None:
        ax.axvline(ref_line, color="red", linewidth=2, linestyle="--",
                   label=ref_label)
    if mean_val is not None:
        ax.axvline(mean_val, color="darkorange", linewidth=2, linestyle="-",
                   label=rf"Mean = {mean_val:.3e}")
    if median_val is not None:
        ax.axvline(median_val, color="green", linewidth=2, linestyle="-.",
                   label=rf"Median = {median_val:.3e}")
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    if title:
        ax.set_title(title, fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, linestyle=":", alpha=0.4)
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ── File grouping ──────────────────────────────────────────────────────

def group_csv_files(data_dir: Path):
    """Split CSV files into glicerolo50 and glicerolo200 groups."""
    all_csv = sorted(data_dir.glob("*.csv"))
    glic50 = [f for f in all_csv if GLIC50_PATTERN in f.name]
    glic200 = [f for f in all_csv if GLIC200_PATTERN in f.name]
    return {"glic50": glic50, "glic200": glic200}


# ── Core analysis — per-track taMSD ────────────────────────────────────

def extract_per_track_D(csv_files):
    """Fit taMSD per track → collect (D_α, α) for every model×fraction.

    Returns a DataFrame with one row per (track, file, model, fraction).
    """
    rows = []

    for csv_path in sorted(csv_files):
        stem = csv_path.stem
        print(f"\n{'='*60}")
        print(f"[taMSD] Processing: {csv_path.name}")
        trajectories = read_trajectories_from_csv(str(csv_path))
        if not trajectories:
            print("  No trajectories — skipping.")
            continue

        global_dt = estimate_global_time_step(trajectories)

        # Global velocity stats for anomalous_drift fits (per-file)
        try:
            vstats = analyze_velocities(trajectories, min_points=MIN_TRACK_POINTS)
        except ValueError:
            vstats = None

        eligible = {tid: t for tid, t in trajectories.items()
                    if t.n_points >= MIN_TRACK_POINTS}
        print(f"  Tracks >= {MIN_TRACK_POINTS} pts: {len(eligible)} / {len(trajectories)}")

        for tid, track in eligible.items():
            for frac in FIT_FRACTIONS:
                pct = int(frac * 100)

                # Compute taMSD for this track
                tamsd = calculate_time_averaged_msd_per_track(
                    track, max_lag_fraction=frac, dt_override=global_dt,
                )
                if tamsd.tau.size < 4:
                    continue

                base_row = dict(
                    file=stem, track_id=str(tid), n_points=track.n_points,
                    fraction=frac, dt=global_dt,
                )

                # ── Anomalous fit (4D_α τ^α) ──────────────────
                try:
                    Da, Da_err, a, a_err, _pcov, tau_f, msd_f, msd_p, chi2, _rss, _sig = \
                        _fit_anomalous_at_fraction(
                            tamsd.tau, tamsd.msd, tamsd.n_max, tamsd.dt,
                            1.0,  # already capped by max_lag_fraction
                            1e-2, (1e-6, 1e2),
                            1.0, (0.01, 2.0),
                            tamsd.msd_sem,
                        )
                    rows.append({
                        **base_row,
                        "model": "anomalous",
                        "D_alpha": Da, "D_alpha_error": Da_err,
                        "alpha": a, "alpha_error": a_err,
                        "chi2_red": chi2,
                        "v": np.nan, "v_error": np.nan,
                    })
                except (ValueError, RuntimeError):
                    pass

                # ── Anomalous+offset fit (4D_α τ^α + c) ───────
                try:
                    fit_ao = fit_msd_anomalous_offset(
                        tamsd.tau, tamsd.msd, tamsd.n_max, tamsd.dt,
                        fit_fraction=1.0,
                        msd_sigma=tamsd.msd_sem,
                    )
                    rows.append({
                        **base_row,
                        "model": "anomalous_offset",
                        "D_alpha": fit_ao.D_alpha, "D_alpha_error": fit_ao.D_alpha_error,
                        "alpha": fit_ao.alpha, "alpha_error": fit_ao.alpha_error,
                        "chi2_red": fit_ao.chi_squared_red,
                        "v": np.nan, "v_error": np.nan,
                        "offset": fit_ao.offset, "offset_error": fit_ao.offset_error,
                    })
                except (ValueError, RuntimeError):
                    pass

                # ── Anomalous+drift fit (4D_α τ^α + v²τ²) ─────
                if vstats is not None:
                    try:
                        Da, Da_err, a, a_err, v, v_err, _pcov, \
                            tau_f, msd_f, msd_p, chi2, _rss, _sig = \
                            _fit_anomalous_drift_at_fraction(
                                tamsd.tau, tamsd.msd,
                                tamsd.n_max, tamsd.dt,
                                1.0,  # already capped
                                1e-2, (1e-6, 1e2),
                                1.0, (0.01, 2.0),
                                vstats.v_initial, vstats.v_bounds,
                                tamsd.msd_sem,
                            )
                        rows.append({
                            **base_row,
                            "model": "anomalous_drift",
                            "D_alpha": Da, "D_alpha_error": Da_err,
                            "alpha": a, "alpha_error": a_err,
                            "chi2_red": chi2,
                            "v": v, "v_error": v_err,
                        })
                    except (ValueError, RuntimeError):
                        pass

        print(f"  Collected {sum(1 for r in rows if r['file'] == stem)} fit results.")

    df = pd.DataFrame(rows)
    print(f"\n[taMSD] Total per-track fits: {len(df)}")
    return df


# ── Core analysis — per-file eaMSD ────────────────────────────────────

def extract_per_file_D(csv_files):
    """Fit eaMSD per file → one (D_α, α) per CSV for each model×fraction.

    Returns a DataFrame with one row per (file, model, fraction).
    """
    rows = []

    for csv_path in sorted(csv_files):
        stem = csv_path.stem
        print(f"\n{'='*60}")
        print(f"[eaMSD] Processing: {csv_path.name}")
        trajectories = read_trajectories_from_csv(str(csv_path))
        if not trajectories:
            print("  No trajectories — skipping.")
            continue

        try:
            vstats = analyze_velocities(trajectories, min_points=MIN_TRACK_POINTS)
        except ValueError:
            vstats = None

        for frac in FIT_FRACTIONS:
            pct = int(frac * 100)

            eamsd = calculate_ensemble_msd(trajectories, max_lag_fraction=frac)
            if eamsd.tau.size < 4:
                continue

            base_row = dict(
                file=stem, n_tracks=eamsd.total_trajectories,
                fraction=frac, dt=eamsd.dt,
            )

            # ── Anomalous fit (4D_α τ^α) ──────────────────────
            try:
                Da, Da_err, a, a_err, _pcov, tau_f, msd_f, msd_p, chi2, _rss, sig = \
                    _fit_anomalous_at_fraction(
                        eamsd.tau, eamsd.msd, eamsd.n_max, eamsd.dt,
                        1.0,  # data already capped
                        1e-2, (1e-6, 1e2),
                        1.0, (0.01, 2.0),
                        eamsd.msd_sem,
                    )
                rows.append({
                    **base_row,
                    "model": "anomalous",
                    "D_alpha": Da, "D_alpha_error": Da_err,
                    "alpha": a, "alpha_error": a_err,
                    "chi2_red": chi2,
                    "v": np.nan, "v_error": np.nan,
                })

                tag = f"f{pct:03d}"
                txt = "\n".join([
                    r"$D_\alpha = (%.2e \pm %.1e)\ \mu m^2/s^\alpha$" % (Da, Da_err),
                    r"$\alpha = %.4f \pm %.4f$" % (a, a_err),
                    r"$\chi^2_\nu = %.4f$" % chi2,
                ])
                _save_fit_plot(
                    tau_f, msd_f, msd_p, txt,
                    r"Anomalous: MSD = 4$D_\alpha\tau^\alpha$",
                    LINEAR_DIR / f"{stem}_eamsd_anomalous_{tag}.svg",
                    sig, data_color="C0",
                )
            except (ValueError, RuntimeError) as e:
                print(f"    anomalous {pct}%: FAILED — {e}")

            # ── Anomalous+offset fit (4D_α τ^α + c) ───────────
            try:
                fit_ao = fit_msd_anomalous_offset(
                    eamsd.tau, eamsd.msd, eamsd.n_max, eamsd.dt,
                    fit_fraction=1.0,
                    msd_sigma=eamsd.msd_sem,
                )
                rows.append({
                    **base_row,
                    "model": "anomalous_offset",
                    "D_alpha": fit_ao.D_alpha, "D_alpha_error": fit_ao.D_alpha_error,
                    "alpha": fit_ao.alpha, "alpha_error": fit_ao.alpha_error,
                    "chi2_red": fit_ao.chi_squared_red,
                    "v": np.nan, "v_error": np.nan,
                    "offset": fit_ao.offset, "offset_error": fit_ao.offset_error,
                })

                tag = f"f{pct:03d}"
                txt_lines = [
                    r"$D_\alpha = (%.2e \pm %.1e)\ \mu m^2/s^\alpha$" % (fit_ao.D_alpha, fit_ao.D_alpha_error),
                    r"$\alpha = %.4f \pm %.4f$" % (fit_ao.alpha, fit_ao.alpha_error),
                    r"$c = (%.2e \pm %.1e)\ \mu m^2$" % (fit_ao.offset, fit_ao.offset_error),
                    r"$\chi^2_\nu = %.4f$" % fit_ao.chi_squared_red,
                ]
                _save_fit_plot(
                    fit_ao.tau_fit, fit_ao.msd_fit, fit_ao.msd_predicted,
                    "\n".join(txt_lines),
                    r"Anomalous+offset: MSD = 4$D_\alpha\tau^\alpha$ + $c$",
                    LINEAR_OFFSET_DIR / f"{stem}_eamsd_anomalous_offset_{tag}.svg",
                    fit_ao.msd_sigma_fit, data_color="C0",
                )
            except (ValueError, RuntimeError) as e:
                print(f"    anomalous_offset {pct}%: FAILED — {e}")

            # ── Anomalous+drift fit (4D_α τ^α + v²τ²) ─────────
            if vstats is not None:
                try:
                    Da, Da_err, a, a_err, v, v_err, _pcov, \
                        tau_f, msd_f, msd_p, chi2, _rss, sig = \
                        _fit_anomalous_drift_at_fraction(
                            eamsd.tau, eamsd.msd,
                            eamsd.n_max, eamsd.dt,
                            1.0,  # data already capped
                            1e-2, (1e-6, 1e2),
                            1.0, (0.01, 2.0),
                            vstats.v_initial, vstats.v_bounds,
                            eamsd.msd_sem,
                        )
                    rows.append({
                        **base_row,
                        "model": "anomalous_drift",
                        "D_alpha": Da, "D_alpha_error": Da_err,
                        "alpha": a, "alpha_error": a_err,
                        "chi2_red": chi2,
                        "v": v, "v_error": v_err,
                    })

                    tag = f"f{pct:03d}"
                    txt = "\n".join([
                        r"$D_\alpha = (%.2e \pm %.1e)\ \mu m^2/s^\alpha$" % (Da, Da_err),
                        r"$\alpha = %.4f \pm %.4f$" % (a, a_err),
                        r"$v = (%.2e \pm %.1e)\ \mu m/s$" % (v, v_err),
                        r"$\chi^2_\nu = %.4f$" % chi2,
                    ])
                    _save_fit_plot(
                        tau_f, msd_f, msd_p, txt,
                        r"Anomalous+drift: MSD = 4$D_\alpha\tau^\alpha$ + $v^2\tau^2$",
                        NONLINEAR_DIR / f"{stem}_eamsd_anomalous_drift_{tag}.svg",
                        sig, data_color="C0",
                    )
                except (ValueError, RuntimeError) as e:
                    print(f"    anomalous_drift {pct}%: FAILED — {e}")

    df = pd.DataFrame(rows)
    print(f"\n[eaMSD] Total per-file fits: {len(df)}")
    return df


# ── Core analysis — per-file ⟨taMSD⟩ ──────────────────────────────────

def extract_ensemble_tamsd_D(csv_files):
    """Fit ensemble-averaged taMSD (⟨taMSD⟩) per file → one (D_α, α) per CSV.

    Returns a DataFrame with one row per (file, model, fraction).
    """
    rows = []

    for csv_path in sorted(csv_files):
        stem = csv_path.stem
        print(f"\n{'='*60}")
        print(f"[<taMSD>] Processing: {csv_path.name}")
        trajectories = read_trajectories_from_csv(str(csv_path))
        if not trajectories:
            print("  No trajectories — skipping.")
            continue

        global_dt = estimate_global_time_step(trajectories)

        try:
            vstats = analyze_velocities(trajectories, min_points=MIN_TRACK_POINTS)
        except ValueError:
            vstats = None

        for frac in FIT_FRACTIONS:
            pct = int(frac * 100)

            tau, msd_mean, msd_sem, _ = compute_ensemble_tamsd(
                trajectories, max_lag_fraction=frac, global_dt=global_dt,
            )
            if tau.size < 4:
                continue

            n_max = tau.size
            base_row = dict(
                file=stem, n_tracks=len(trajectories),
                fraction=frac, dt=global_dt,
            )

            # ── Anomalous fit ──────────────────────────────────
            try:
                Da, Da_err, a, a_err, _pcov, tau_f, msd_f, msd_p, chi2, _rss, sig = \
                    _fit_anomalous_at_fraction(
                        tau, msd_mean, n_max, global_dt,
                        1.0,
                        1e-2, (1e-6, 1e2),
                        1.0, (0.01, 2.0),
                        msd_sem,
                    )
                rows.append({
                    **base_row,
                    "model": "anomalous",
                    "D_alpha": Da, "D_alpha_error": Da_err,
                    "alpha": a, "alpha_error": a_err,
                    "chi2_red": chi2,
                    "v": np.nan, "v_error": np.nan,
                })

                tag = f"f{pct:03d}"
                txt = "\n".join([
                    r"$D_\alpha = (%.2e \pm %.1e)\ \mu m^2/s^\alpha$" % (Da, Da_err),
                    r"$\alpha = %.4f \pm %.4f$" % (a, a_err),
                    r"$\chi^2_\nu = %.4f$" % chi2,
                ])
                _save_fit_plot(
                    tau_f, msd_f, msd_p, txt,
                    r"Anomalous: MSD = 4$D_\alpha\tau^\alpha$",
                    LINEAR_DIR / f"{stem}_ens_tamsd_anomalous_{tag}.svg",
                    sig, data_color="C2",
                )
            except (ValueError, RuntimeError) as e:
                print(f"    anomalous {pct}%: FAILED — {e}")

            # ── Anomalous+offset fit ───────────────────────────
            try:
                fit_ao = fit_msd_anomalous_offset(
                    tau, msd_mean, n_max, global_dt,
                    fit_fraction=1.0,
                    msd_sigma=msd_sem,
                )
                rows.append({
                    **base_row,
                    "model": "anomalous_offset",
                    "D_alpha": fit_ao.D_alpha, "D_alpha_error": fit_ao.D_alpha_error,
                    "alpha": fit_ao.alpha, "alpha_error": fit_ao.alpha_error,
                    "chi2_red": fit_ao.chi_squared_red,
                    "v": np.nan, "v_error": np.nan,
                    "offset": fit_ao.offset, "offset_error": fit_ao.offset_error,
                })

                tag = f"f{pct:03d}"
                txt_lines = [
                    r"$D_\alpha = (%.2e \pm %.1e)\ \mu m^2/s^\alpha$" % (fit_ao.D_alpha, fit_ao.D_alpha_error),
                    r"$\alpha = %.4f \pm %.4f$" % (fit_ao.alpha, fit_ao.alpha_error),
                    r"$c = (%.2e \pm %.1e)\ \mu m^2$" % (fit_ao.offset, fit_ao.offset_error),
                    r"$\chi^2_\nu = %.4f$" % fit_ao.chi_squared_red,
                ]
                _save_fit_plot(
                    fit_ao.tau_fit, fit_ao.msd_fit, fit_ao.msd_predicted,
                    "\n".join(txt_lines),
                    r"Anomalous+offset: MSD = 4$D_\alpha\tau^\alpha$ + $c$",
                    LINEAR_OFFSET_DIR / f"{stem}_ens_tamsd_anomalous_offset_{tag}.svg",
                    fit_ao.msd_sigma_fit, data_color="C2",
                )
            except (ValueError, RuntimeError) as e:
                print(f"    anomalous_offset {pct}%: FAILED — {e}")

            # ── Anomalous+drift fit ────────────────────────────
            if vstats is not None:
                try:
                    Da, Da_err, a, a_err, v, v_err, _pcov, \
                        tau_f, msd_f, msd_p, chi2, _rss, sig = \
                        _fit_anomalous_drift_at_fraction(
                            tau, msd_mean,
                            n_max, global_dt,
                            1.0,
                            1e-2, (1e-6, 1e2),
                            1.0, (0.01, 2.0),
                            vstats.v_initial, vstats.v_bounds,
                            msd_sem,
                        )
                    rows.append({
                        **base_row,
                        "model": "anomalous_drift",
                        "D_alpha": Da, "D_alpha_error": Da_err,
                        "alpha": a, "alpha_error": a_err,
                        "chi2_red": chi2,
                        "v": v, "v_error": v_err,
                    })

                    tag = f"f{pct:03d}"
                    txt = "\n".join([
                        r"$D_\alpha = (%.2e \pm %.1e)\ \mu m^2/s^\alpha$" % (Da, Da_err),
                        r"$\alpha = %.4f \pm %.4f$" % (a, a_err),
                        r"$v = (%.2e \pm %.1e)\ \mu m/s$" % (v, v_err),
                        r"$\chi^2_\nu = %.4f$" % chi2,
                    ])
                    _save_fit_plot(
                        tau_f, msd_f, msd_p, txt,
                        r"Anomalous+drift: MSD = 4$D_\alpha\tau^\alpha$ + $v^2\tau^2$",
                        NONLINEAR_DIR / f"{stem}_ens_tamsd_anomalous_drift_{tag}.svg",
                        sig, data_color="C2",
                    )
                except (ValueError, RuntimeError) as e:
                    print(f"    anomalous_drift {pct}%: FAILED — {e}")

    df = pd.DataFrame(rows)
    print(f"\n[<taMSD>] Total per-file fits: {len(df)}")
    return df


# ── Statistics helpers ─────────────────────────────────────────────────

def compute_statistics(values):
    """Return dict of descriptive statistics for an array of values."""
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    n = len(arr)
    if n == 0:
        return dict(n=0, mean=np.nan, median=np.nan, std=np.nan, sem=np.nan)
    return dict(
        n=n,
        mean=float(np.mean(arr)),
        median=float(np.median(arr)),
        std=float(np.std(arr, ddof=1)),
        sem=float(np.std(arr, ddof=1) / np.sqrt(n)),
    )


def one_sample_ttest(values, reference):
    """One-sample t-test: H0: mean(values) = reference.

    Returns (t_statistic, p_value) or (nan, nan) if insufficient data.
    """
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) < 2:
        return float("nan"), float("nan")
    t_stat, p_val = stats.ttest_1samp(arr, reference)
    return float(t_stat), float(p_val)


# ── Full analysis for a group ─────────────────────────────────────────

def analyze_group(group_name, csv_files, group_suffix=""):
    """Run full anomalous diffusion analysis for one dataset group.

    Returns summary_rows list.
    """
    if not csv_files:
        print(f"  No files for group {group_name}")
        return []

    print(f"\n{'#'*70}")
    print(f"  GROUP: {group_name}  ({len(csv_files)} files)")
    print(f"{'#'*70}")

    sfx = f"_{group_suffix}" if group_suffix else ""

    # Phase B: per-track taMSD
    df_tamsd = extract_per_track_D(csv_files)

    # Phase C: per-file eaMSD
    df_eamsd = extract_per_file_D(csv_files)

    # Phase C2: per-file ⟨taMSD⟩
    df_ens_tamsd = extract_ensemble_tamsd_D(csv_files)

    # Phase D+E: histograms, statistics
    summary_rows = []

    for model in ("anomalous", "anomalous_offset", "anomalous_drift"):
        for frac in FIT_FRACTIONS:
            pct = int(frac * 100)
            tag = f"{model}_f{pct:03d}{sfx}"

            # ── taMSD distribution ─────────────────────────────
            mask_ta = (df_tamsd["model"] == model) & (df_tamsd["fraction"] == frac)
            Da_ta = df_tamsd.loc[mask_ta, "D_alpha"].values if not df_tamsd.empty else np.array([])
            a_ta = df_tamsd.loc[mask_ta, "alpha"].values if not df_tamsd.empty else np.array([])

            st_Da_ta = compute_statistics(Da_ta)
            st_a_ta = compute_statistics(a_ta)
            t_stat_a, p_val_a = one_sample_ttest(a_ta, 1.0)  # H0: α = 1 (normal diffusion)

            print(f"\n{'-'*60}")
            print(f"taMSD | {model} | {pct}% lag fraction | {group_name}")
            print(f"  N tracks      = {st_Da_ta['n']}")
            print(f"  Mean(D_α)     = {st_Da_ta['mean']:.4e} µm²/s^α")
            print(f"  Median(D_α)   = {st_Da_ta['median']:.4e}")
            print(f"  Mean(α)       = {st_a_ta['mean']:.4f}")
            print(f"  Median(α)     = {st_a_ta['median']:.4f}")
            print(f"  t-test α vs 1: t = {t_stat_a:.3f}, p = {p_val_a:.4e}")

            if st_Da_ta["n"] > 0:
                Da_fin = Da_ta[np.isfinite(Da_ta)]
                _plot_histogram(
                    Da_fin,
                    xlabel=r"$D_\alpha$ [$\mu$m$^2$/s$^\alpha$]",
                    output_path=STATS_DIR / f"tamsd_Dalpha_histogram_{tag}.svg",
                    mean_val=st_Da_ta["mean"], median_val=st_Da_ta["median"],
                )
                a_fin = a_ta[np.isfinite(a_ta)]
                _plot_histogram(
                    a_fin,
                    xlabel=r"$\alpha$",
                    output_path=STATS_DIR / f"tamsd_alpha_histogram_{tag}.svg",
                    mean_val=st_a_ta["mean"], median_val=st_a_ta["median"],
                    ref_line=1.0,
                    ref_label=r"$\alpha = 1$ (normal)",
                )

            # ── eaMSD distribution ─────────────────────────────
            mask_ea = (df_eamsd["model"] == model) & (df_eamsd["fraction"] == frac)
            Da_ea = df_eamsd.loc[mask_ea, "D_alpha"].values if not df_eamsd.empty else np.array([])
            a_ea = df_eamsd.loc[mask_ea, "alpha"].values if not df_eamsd.empty else np.array([])

            st_Da_ea = compute_statistics(Da_ea)
            st_a_ea = compute_statistics(a_ea)
            t_stat_a_ea, p_val_a_ea = one_sample_ttest(a_ea, 1.0)

            print(f"\neaMSD | {model} | {pct}% | {group_name}")
            print(f"  N files       = {st_Da_ea['n']}")
            print(f"  Mean(D_α)     = {st_Da_ea['mean']:.4e}")
            print(f"  Mean(α)       = {st_a_ea['mean']:.4f}")

            if st_Da_ea["n"] > 0:
                _plot_histogram(
                    Da_ea[np.isfinite(Da_ea)],
                    xlabel=r"$D_\alpha$ [$\mu$m$^2$/s$^\alpha$]",
                    output_path=STATS_DIR / f"eamsd_Dalpha_histogram_{tag}.svg",
                    mean_val=st_Da_ea["mean"], median_val=st_Da_ea["median"],
                )
                _plot_histogram(
                    a_ea[np.isfinite(a_ea)],
                    xlabel=r"$\alpha$",
                    output_path=STATS_DIR / f"eamsd_alpha_histogram_{tag}.svg",
                    mean_val=st_a_ea["mean"], median_val=st_a_ea["median"],
                    ref_line=1.0,
                    ref_label=r"$\alpha = 1$ (normal)",
                )

            # ── ⟨taMSD⟩ distribution ───────────────────────────
            mask_et = (df_ens_tamsd["model"] == model) & (df_ens_tamsd["fraction"] == frac)
            Da_et = df_ens_tamsd.loc[mask_et, "D_alpha"].values if not df_ens_tamsd.empty else np.array([])
            a_et = df_ens_tamsd.loc[mask_et, "alpha"].values if not df_ens_tamsd.empty else np.array([])

            st_Da_et = compute_statistics(Da_et)
            st_a_et = compute_statistics(a_et)
            t_stat_a_et, p_val_a_et = one_sample_ttest(a_et, 1.0)

            print(f"\n<taMSD> | {model} | {pct}% | {group_name}")
            print(f"  N files       = {st_Da_et['n']}")
            print(f"  Mean(D_α)     = {st_Da_et['mean']:.4e}")
            print(f"  Mean(α)       = {st_a_et['mean']:.4f}")

            if st_Da_et["n"] > 0:
                _plot_histogram(
                    Da_et[np.isfinite(Da_et)],
                    xlabel=r"$D_\alpha$ [$\mu$m$^2$/s$^\alpha$]",
                    output_path=STATS_DIR / f"ens_tamsd_Dalpha_histogram_{tag}.svg",
                    mean_val=st_Da_et["mean"], median_val=st_Da_et["median"],
                )
                _plot_histogram(
                    a_et[np.isfinite(a_et)],
                    xlabel=r"$\alpha$",
                    output_path=STATS_DIR / f"ens_tamsd_alpha_histogram_{tag}.svg",
                    mean_val=st_a_et["mean"], median_val=st_a_et["median"],
                    ref_line=1.0,
                    ref_label=r"$\alpha = 1$ (normal)",
                )

            summary_rows.append(dict(
                group=group_name,
                model=model, fraction=frac,
                # taMSD per-track
                Da_tamsd_mean=st_Da_ta["mean"], Da_tamsd_std=st_Da_ta["std"],
                Da_tamsd_sem=st_Da_ta["sem"],
                alpha_tamsd_mean=st_a_ta["mean"], alpha_tamsd_std=st_a_ta["std"],
                N_tracks=st_Da_ta["n"],
                t_stat_alpha_tamsd=t_stat_a, p_val_alpha_tamsd=p_val_a,
                # eaMSD per-file
                Da_eamsd_mean=st_Da_ea["mean"], Da_eamsd_std=st_Da_ea["std"],
                Da_eamsd_sem=st_Da_ea["sem"],
                alpha_eamsd_mean=st_a_ea["mean"], alpha_eamsd_std=st_a_ea["std"],
                N_files_ea=st_Da_ea["n"],
                t_stat_alpha_eamsd=t_stat_a_ea, p_val_alpha_eamsd=p_val_a_ea,
                # ⟨taMSD⟩ per-file
                Da_ens_tamsd_mean=st_Da_et["mean"], Da_ens_tamsd_std=st_Da_et["std"],
                Da_ens_tamsd_sem=st_Da_et["sem"],
                alpha_ens_tamsd_mean=st_a_et["mean"], alpha_ens_tamsd_std=st_a_et["std"],
                N_files_et=st_Da_et["n"],
                t_stat_alpha_ens_tamsd=t_stat_a_et, p_val_alpha_ens_tamsd=p_val_a_et,
            ))

    # Save per-group CSVs
    if not df_tamsd.empty:
        tamsd_csv = DOC_DIR / f"diffusion_statistics_anomalous_tamsd_per_track{sfx}.csv"
        df_tamsd.to_csv(tamsd_csv, index=False, float_format="%.6e")
        print(f"Per-track CSV saved to {tamsd_csv}")

    if not df_eamsd.empty:
        eamsd_csv = DOC_DIR / f"diffusion_statistics_anomalous_eamsd_per_file{sfx}.csv"
        df_eamsd.to_csv(eamsd_csv, index=False, float_format="%.6e")
        print(f"Per-file eaMSD CSV saved to {eamsd_csv}")

    if not df_ens_tamsd.empty:
        et_csv = DOC_DIR / f"diffusion_statistics_anomalous_ens_tamsd_per_file{sfx}.csv"
        df_ens_tamsd.to_csv(et_csv, index=False, float_format="%.6e")
        print(f"Per-file <taMSD> CSV saved to {et_csv}")

    return summary_rows


# ── Main ───────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description="Anomalous diffusion coefficient statistical analysis",
    )
    parser.add_argument(
        "--group", "-g", type=str, default=None,
        choices=["glic50", "glic200"],
        help="Analyse only one glycerol group. If omitted, both are analysed.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Create output dirs
    for d in (STATS_DIR, LINEAR_DIR, LINEAR_OFFSET_DIR, NONLINEAR_DIR, DOC_DIR):
        d.mkdir(parents=True, exist_ok=True)

    groups = group_csv_files(DATA_DIR)
    all_summary = []

    if args.group is None or args.group == "glic50":
        rows = analyze_group("glic50", groups["glic50"], group_suffix="glic50")
        all_summary.extend(rows)

    if args.group is None or args.group == "glic200":
        rows = analyze_group("glic200", groups["glic200"], group_suffix="glic200")
        all_summary.extend(rows)

    # Save combined summary
    if all_summary:
        df_summary = pd.DataFrame(all_summary)
        csv_out = DOC_DIR / "diffusion_statistics_anomalous_results.csv"
        df_summary.to_csv(csv_out, index=False, float_format="%.6e")
        print(f"\nSummary CSV saved to {csv_out}")

    # ── Final summary table ─────────────────────────────────────
    print(f"\n{'='*140}")
    print("FINAL SUMMARY — ANOMALOUS DIFFUSION")
    print(f"{'='*140}")

    header = (f"{'Group':<8} {'Model':<18} {'Frac':>5} | "
              f"{'<D_α>_ta':>11} {'<α>_ta':>8} {'N':>5} {'p(α=1)':>10} | "
              f"{'<D_α>_ea':>11} {'<α>_ea':>8} {'N':>5} {'p(α=1)':>10} | "
              f"{'<D_α>_⟨ta⟩':>11} {'<α>_⟨ta⟩':>8} {'N':>5} {'p(α=1)':>10}")
    print(header)
    print("-" * len(header))
    for r in all_summary:
        pct = f"{int(r['fraction']*100)}%"
        print(
            f"{r['group']:<8} {r['model']:<18} {pct:>5} | "
            f"{r['Da_tamsd_mean']:>11.4e} {r['alpha_tamsd_mean']:>8.4f} "
            f"{r['N_tracks']:>5} {r['p_val_alpha_tamsd']:>10.4e} | "
            f"{r['Da_eamsd_mean']:>11.4e} {r['alpha_eamsd_mean']:>8.4f} "
            f"{r['N_files_ea']:>5} {r['p_val_alpha_eamsd']:>10.4e} | "
            f"{r['Da_ens_tamsd_mean']:>11.4e} {r['alpha_ens_tamsd_mean']:>8.4f} "
            f"{r['N_files_et']:>5} {r['p_val_alpha_ens_tamsd']:>10.4e}"
        )
    print(f"{'='*140}")
    print("Done.")


if __name__ == "__main__":
    main()
