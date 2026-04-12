"""
Rigorous statistical analysis of the diffusion coefficient D.

Workflow
--------
1. Per-track taMSD fitting  → individual D_i values (pooled across all CSV files)
2. Per-file  eaMSD fitting  → one D per CSV file
3. Histograms of D for both estimators (taMSD, eaMSD)
4. One-sample t-test of mean(D) vs D_theory (Einstein-Stokes)
5. Summary comparison plot D_taMSD vs D_eaMSD vs D_theory

Models:  linear (4Dτ),  nonlinear/drift (4Dτ + v²τ²)
Fit fractions: 10 % and 25 %
Dataset: non-anomalous only (Data/31_10_no_anomalous)

Output (SVG for paper, CSV for data):
    Results/no_anomalous/diffusion_statistics/   — histograms, comparison plot
    Results/no_anomalous/linear_fits/            — per-track linear fit plots
    Results/no_anomalous/nonlinear_fits/         — per-track nonlinear fit plots
    Docu/diffusion_statistics_results.csv        — tabular results

Usage:
    python diffusion_statistics.py
"""

from __future__ import annotations

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
    fit_msd_linear,
    linear_msd_model,
    fit_msd_nonlinear,
    nonlinear_msd_model,
    analyze_velocities,
    _fit_nonlinear_at_fraction,
)
from comp_glycerol_viscosity import diffusion_coefficient

# ── Configuration ──────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "Data" / "31_10_no_anomalous"
STATS_DIR = SCRIPT_DIR / "Results" / "no_anomalous" / "diffusion_statistics"
LINEAR_DIR = SCRIPT_DIR / "Results" / "no_anomalous" / "linear_fits"
NONLINEAR_DIR = SCRIPT_DIR / "Results" / "no_anomalous" / "nonlinear_fits"
DOC_DIR = SCRIPT_DIR / "Docu"

MIN_TRACK_POINTS = 30
FIT_FRACTIONS = [0.10, 0.25]

# Experimental parameters for D_theory (Einstein-Stokes)
T_C = 23.3       # Temperature [°C]
C_M = 0.85       # Glycerol mass fraction
R_PARTICLE = 120e-9  # Particle radius [m]


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


def _plot_histogram(D_values, D_theory, xlabel, output_path,
                    mean_val=None, median_val=None, title=None):
    """Plot histogram of D values with theory, mean, and median lines."""
    fig, ax = plt.subplots(figsize=(8, 5.5))
    n_bins = max(10, int(np.sqrt(len(D_values))))
    ax.hist(D_values, bins=n_bins, color="steelblue", edgecolor="white",
            alpha=0.8, density=True, label="Data")
    ax.axvline(D_theory, color="red", linewidth=2, linestyle="--",
               label=rf"$D_{{theory}}$ = {D_theory:.3e}")
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

# ── Core analysis ──────────────────────────────────────────────────────

def extract_per_track_D(csv_files):
    """Phase B: fit taMSD per track → collect D_i for every model×fraction.

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

        # Global velocity stats for nonlinear fits (per-file)
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
                if tamsd.tau.size < 3:
                    continue

                base_row = dict(
                    file=stem, track_id=str(tid), n_points=track.n_points,
                    fraction=frac, dt=global_dt,
                )

                # ── Linear fit ─────────────────────────────────────
                try:
                    fit_lin = fit_msd_linear(
                        tamsd.tau, tamsd.msd, tamsd.n_max, tamsd.dt,
                        fit_fraction=1.0,  # already capped by max_lag_fraction
                        msd_sigma=tamsd.msd_sem,
                    )
                    rows.append({
                        **base_row,
                        "model": "linear",
                        "D": fit_lin.D,
                        "D_error": fit_lin.D_error,
                        "chi2_red": fit_lin.chi_squared_red,
                        "v": np.nan,
                        "v_error": np.nan,
                    })
                except (ValueError, RuntimeError):
                    pass

                # ── Nonlinear (drift) fit ──────────────────────────
                if vstats is not None:
                    try:
                        (D, D_err, v, v_err, _pcov,
                         tau_f, msd_f, msd_p, chi2, rss,
                         _sig) = _fit_nonlinear_at_fraction(
                            tamsd.tau, tamsd.msd,
                            tamsd.n_max, tamsd.dt,
                            1.0,  # full fraction (data already capped)
                            1e-2,                      # D_initial
                            (9e-4, 1.5e-1),            # D_bounds
                            vstats.v_initial,
                            vstats.v_bounds,
                            tamsd.msd_sem,
                        )
                        rows.append({
                            **base_row,
                            "model": "nonlinear",
                            "D": D,
                            "D_error": D_err,
                            "chi2_red": chi2,
                            "v": v,
                            "v_error": v_err,
                        })
                    except (ValueError, RuntimeError):
                        pass

        print(f"  Collected {sum(1 for r in rows if r['file'] == stem)} fit results.")

    df = pd.DataFrame(rows)
    print(f"\n[taMSD] Total per-track fits: {len(df)}")
    return df


def extract_per_file_D(csv_files):
    """Phase C: fit eaMSD per file → one D per CSV for each model×fraction.

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
            if eamsd.tau.size < 3:
                continue

            base_row = dict(
                file=stem, n_tracks=eamsd.total_trajectories,
                fraction=frac, dt=eamsd.dt,
            )

            # ── Linear fit ─────────────────────────────────────
            try:
                fit_lin = fit_msd_linear(
                    eamsd.tau, eamsd.msd, eamsd.n_max, eamsd.dt,
                    fit_fraction=1.0,  # data already capped
                    msd_sigma=eamsd.msd_sem,
                )
                rows.append({
                    **base_row,
                    "model": "linear",
                    "D": fit_lin.D,
                    "D_error": fit_lin.D_error,
                    "chi2_red": fit_lin.chi_squared_red,
                    "v": np.nan,
                    "v_error": np.nan,
                })

                # Save fit plot
                tag = f"f{pct:03d}"
                txt = (r"$D = (%.2e \pm %.1e)\ \mu m^2/s$" % (fit_lin.D, fit_lin.D_error)
                       + "\n" + r"$\chi^2_\nu = %.4f$" % fit_lin.chi_squared_red)
                _save_fit_plot(
                    fit_lin.tau_fit, fit_lin.msd_fit, fit_lin.msd_predicted,
                    txt, r"Linear: MSD = 4D$\tau$",
                    LINEAR_DIR / f"{stem}_eamsd_linear_{tag}.svg",
                    fit_lin.msd_sigma_fit,
                    #title=f"eaMSD — {stem} — linear {pct}%",
                    data_color="C0",
                )
            except (ValueError, RuntimeError) as e:
                print(f"    linear {pct}%: FAILED — {e}")

            # ── Nonlinear (drift) fit ──────────────────────────
            if vstats is not None:
                try:
                    (D, D_err, v, v_err, _pcov,
                     tau_f, msd_f, msd_p, chi2, rss,
                     sig) = _fit_nonlinear_at_fraction(
                        eamsd.tau, eamsd.msd,
                        eamsd.n_max, eamsd.dt,
                        1.0,  # full fraction (data already capped)
                        1e-2,                      # D_initial
                        (9e-4, 1.5e-1),            # D_bounds
                        vstats.v_initial,
                        vstats.v_bounds,
                        eamsd.msd_sem,
                    )
                    rows.append({
                        **base_row,
                        "model": "nonlinear",
                        "D": D,
                        "D_error": D_err,
                        "chi2_red": chi2,
                        "v": v,
                        "v_error": v_err,
                    })

                    tag = f"f{pct:03d}"
                    txt = "\n".join([
                        r"$D = (%.2e \pm %.1e)\ \mu m^2/s$" % (D, D_err),
                        r"$v = (%.2e \pm %.1e)\ \mu m/s$" % (v, v_err),
                        r"$\chi^2_\nu = %.4f$" % chi2,
                    ])
                    _save_fit_plot(
                        tau_f, msd_f, msd_p, txt,
                        r"Nonlinear: MSD = 4D$\tau$ + $v^2\tau^2$",
                        NONLINEAR_DIR / f"{stem}_eamsd_nonlinear_{tag}.svg",
                        sig,
                        data_color="C0",
                    )
                except (ValueError, RuntimeError) as e:
                    print(f"    nonlinear {pct}%: FAILED — {e}")

    df = pd.DataFrame(rows)
    print(f"\n[eaMSD] Total per-file fits: {len(df)}")
    return df


def extract_ensemble_tamsd_D(csv_files):
    """Fit ensemble-averaged taMSD (<taMSD>) per file -> one D per CSV.

    Returns a DataFrame with one row per (file, model, fraction).
    """
    rows = []

    for csv_path in sorted(csv_files):
        stem = csv_path.stem
        print(f"\n{'='*60}")
        print(f"[<taMSD>] Processing: {csv_path.name}")
        trajectories = read_trajectories_from_csv(str(csv_path))
        if not trajectories:
            print("  No trajectories -- skipping.")
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
            if tau.size < 3:
                continue

            n_max = tau.size
            base_row = dict(
                file=stem, n_tracks=len(trajectories),
                fraction=frac, dt=global_dt,
            )

            # -- Linear fit -----------------------------------------
            try:
                fit_lin = fit_msd_linear(
                    tau, msd_mean, n_max, global_dt,
                    fit_fraction=1.0,
                    msd_sigma=msd_sem,
                )
                rows.append({
                    **base_row,
                    "model": "linear",
                    "D": fit_lin.D,
                    "D_error": fit_lin.D_error,
                    "chi2_red": fit_lin.chi_squared_red,
                    "v": np.nan,
                    "v_error": np.nan,
                })

                tag = f"f{pct:03d}"
                txt = (r"$D = (%.2e \pm %.1e)\ \mu m^2/s$" % (fit_lin.D, fit_lin.D_error)
                       + "\n" + r"$\chi^2_\nu = %.4f$" % fit_lin.chi_squared_red)
                _save_fit_plot(
                    fit_lin.tau_fit, fit_lin.msd_fit, fit_lin.msd_predicted,
                    txt, r"Linear: MSD = 4D$\tau$",
                    LINEAR_DIR / f"{stem}_ens_tamsd_linear_{tag}.svg",
                    fit_lin.msd_sigma_fit,
                    data_color="C2",
                )
            except (ValueError, RuntimeError) as e:
                print(f"    linear {pct}%: FAILED -- {e}")

            # -- Nonlinear (drift) fit ------------------------------
            if vstats is not None:
                try:
                    (D, D_err, v, v_err, _pcov,
                     tau_f, msd_f, msd_p, chi2, rss,
                     sig) = _fit_nonlinear_at_fraction(
                        tau, msd_mean,
                        n_max, global_dt,
                        1.0,
                        1e-2,
                        (9e-4, 1.5e-1),
                        vstats.v_initial,
                        vstats.v_bounds,
                        msd_sem,
                    )
                    rows.append({
                        **base_row,
                        "model": "nonlinear",
                        "D": D,
                        "D_error": D_err,
                        "chi2_red": chi2,
                        "v": v,
                        "v_error": v_err,
                    })

                    tag = f"f{pct:03d}"
                    txt = "\n".join([
                        r"$D = (%.2e \pm %.1e)\ \mu m^2/s$" % (D, D_err),
                        r"$v = (%.2e \pm %.1e)\ \mu m/s$" % (v, v_err),
                        r"$\chi^2_\nu = %.4f$" % chi2,
                    ])
                    _save_fit_plot(
                        tau_f, msd_f, msd_p, txt,
                        r"Nonlinear: MSD = 4D$\tau$ + $v^2\tau^2$",
                        NONLINEAR_DIR / f"{stem}_ens_tamsd_nonlinear_{tag}.svg",
                        sig,
                        data_color="C2",
                    )
                except (ValueError, RuntimeError) as e:
                    print(f"    nonlinear {pct}%: FAILED -- {e}")

    df = pd.DataFrame(rows)
    print(f"\n[<taMSD>] Total per-file fits: {len(df)}")
    return df


def compute_statistics(D_values):
    """Return dict of descriptive statistics for an array of D values."""
    arr = np.asarray(D_values, dtype=float)
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


def one_sample_ttest(D_values, D_theory):
    """One-sample t-test: H0: mean(D_i) = D_theory.

    Returns (t_statistic, p_value) or (nan, nan) if insufficient data.
    """
    arr = np.asarray(D_values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) < 2:
        return float("nan"), float("nan")
    t_stat, p_val = stats.ttest_1samp(arr, D_theory)
    return float(t_stat), float(p_val)


# ── Main ───────────────────────────────────────────────────────────────

def main():
    # Create output dirs
    for d in (STATS_DIR, LINEAR_DIR, NONLINEAR_DIR, DOC_DIR):
        d.mkdir(parents=True, exist_ok=True)

    # D_theory in µm²/s
    D_theory_m2s = diffusion_coefficient(T_C, C_M, R_PARTICLE)
    D_theory = D_theory_m2s * 1e12  # m²/s → µm²/s
    print(f"D_theory (Einstein-Stokes) = {D_theory:.6e} um^2/s")
    print(f"  ({D_theory_m2s:.6e} m^2/s)")

    # Gather CSV files
    csv_files = sorted(DATA_DIR.glob("*.csv"))
    if not csv_files:
        print(f"No CSV files in {DATA_DIR}")
        return
    print(f"\nFound {len(csv_files)} CSV files in {DATA_DIR.name}/")

    # Phase B: per-track taMSD D_i
    df_tamsd = extract_per_track_D(csv_files)

    # Phase C: per-file eaMSD D
    df_eamsd = extract_per_file_D(csv_files)

    # Phase C2: per-file ensemble-averaged taMSD D
    df_ens_tamsd = extract_ensemble_tamsd_D(csv_files)

    # Phase D+E: histograms, statistics, and t-tests
    summary_rows = []

    for model in ("linear", "nonlinear"):
        for frac in FIT_FRACTIONS:
            pct = int(frac * 100)
            tag = f"{model}_f{pct:03d}"

            # ── taMSD distribution ─────────────────────────────
            mask_ta = (df_tamsd["model"] == model) & (df_tamsd["fraction"] == frac)
            D_ta = df_tamsd.loc[mask_ta, "D"].values

            st_ta = compute_statistics(D_ta)
            t_stat, p_val = one_sample_ttest(D_ta, D_theory)

            print(f"\n{'-'*60}")
            print(f"taMSD | {model} | {pct}% lag fraction")
            print(f"  N tracks    = {st_ta['n']}")
            print(f"  Mean(D)     = {st_ta['mean']:.4e} um^2/s")
            print(f"  Median(D)   = {st_ta['median']:.4e} um^2/s")
            print(f"  Std(D)      = {st_ta['std']:.4e} um^2/s")
            print(f"  SEM(D)      = {st_ta['sem']:.4e} um^2/s")
            print(f"  D_theory    = {D_theory:.4e} um^2/s")
            print(f"  t-test: t = {t_stat:.3f}, p = {p_val:.4e}")

            if st_ta["n"] > 0:
                _plot_histogram(
                    D_ta[np.isfinite(D_ta)], D_theory,
                    xlabel=r"$D$ [$\mu$m$^2$/s]",
                    output_path=STATS_DIR / f"tamsd_D_histogram_{tag}.svg",
                    mean_val=st_ta["mean"],
                    median_val=st_ta["median"],
                )

            # ── eaMSD distribution ─────────────────────────────
            mask_ea = (df_eamsd["model"] == model) & (df_eamsd["fraction"] == frac)
            D_ea = df_eamsd.loc[mask_ea, "D"].values

            st_ea = compute_statistics(D_ea)
            t_stat_ea, p_val_ea = one_sample_ttest(D_ea, D_theory)

            print(f"\neaMSD | {model} | {pct}% lag fraction")
            print(f"  N files     = {st_ea['n']}")
            print(f"  Mean(D)     = {st_ea['mean']:.4e} um^2/s")
            print(f"  Median(D)   = {st_ea['median']:.4e} um^2/s")
            print(f"  Std(D)      = {st_ea['std']:.4e} um^2/s")
            print(f"  SEM(D)      = {st_ea['sem']:.4e} um^2/s")
            print(f"  D_theory    = {D_theory:.4e} um^2/s")
            print(f"  t-test: t = {t_stat_ea:.3f}, p = {p_val_ea:.4e}")

            if st_ea["n"] > 0:
                _plot_histogram(
                    D_ea[np.isfinite(D_ea)], D_theory,
                    xlabel=r"$D$ [$\mu$m$^2$/s]",
                    output_path=STATS_DIR / f"eamsd_D_histogram_{tag}.svg",
                    mean_val=st_ea["mean"],
                    median_val=st_ea["median"],
                )

            # ── <taMSD> distribution ───────────────────────────
            mask_et = (df_ens_tamsd["model"] == model) & (df_ens_tamsd["fraction"] == frac)
            D_et = df_ens_tamsd.loc[mask_et, "D"].values

            st_et = compute_statistics(D_et)
            t_stat_et, p_val_et = one_sample_ttest(D_et, D_theory)

            print(f"\n<taMSD> | {model} | {pct}% lag fraction")
            print(f"  N files     = {st_et['n']}")
            print(f"  Mean(D)     = {st_et['mean']:.4e} um^2/s")
            print(f"  Median(D)   = {st_et['median']:.4e} um^2/s")
            print(f"  Std(D)      = {st_et['std']:.4e} um^2/s")
            print(f"  SEM(D)      = {st_et['sem']:.4e} um^2/s")
            print(f"  D_theory    = {D_theory:.4e} um^2/s")
            print(f"  t-test: t = {t_stat_et:.3f}, p = {p_val_et:.4e}")

            if st_et["n"] > 0:
                _plot_histogram(
                    D_et[np.isfinite(D_et)], D_theory,
                    xlabel=r"$D$ [$\mu$m$^2$/s]",
                    output_path=STATS_DIR / f"ens_tamsd_D_histogram_{tag}.svg",
                    mean_val=st_et["mean"],
                    median_val=st_et["median"],
                )

            # ── Comparisons ────────────────────────────────────
            if st_ta["n"] > 0 and st_ea["n"] > 0:
                rel_diff = abs(st_ta["mean"] - st_ea["mean"]) / D_theory
                print(f"\n  |<D>_taMSD - <D>_eaMSD| / D_theory = {rel_diff:.4f}")
            if st_ta["n"] > 0 and st_et["n"] > 0:
                rel_diff2 = abs(st_ta["mean"] - st_et["mean"]) / D_theory
                print(f"  |<D>_taMSD - <D>_<taMSD>| / D_theory = {rel_diff2:.4f}")

            summary_rows.append(dict(
                model=model, fraction=frac,
                D_tamsd_mean=st_ta["mean"], D_tamsd_std=st_ta["std"],
                D_tamsd_sem=st_ta["sem"], D_tamsd_median=st_ta["median"],
                N_tracks=st_ta["n"],
                t_stat_tamsd=t_stat, p_val_tamsd=p_val,
                D_eamsd_mean=st_ea["mean"], D_eamsd_std=st_ea["std"],
                D_eamsd_sem=st_ea["sem"], D_eamsd_median=st_ea["median"],
                N_files_ea=st_ea["n"],
                t_stat_eamsd=t_stat_ea, p_val_eamsd=p_val_ea,
                D_ens_tamsd_mean=st_et["mean"], D_ens_tamsd_std=st_et["std"],
                D_ens_tamsd_sem=st_et["sem"], D_ens_tamsd_median=st_et["median"],
                N_files_et=st_et["n"],
                t_stat_ens_tamsd=t_stat_et, p_val_ens_tamsd=p_val_et,
                D_theory=D_theory,
            ))

    # Phase F: comparison plot
    valid_rows = [r for r in summary_rows
                  if np.isfinite(r["D_tamsd_mean"]) or np.isfinite(r["D_eamsd_mean"])]

    # Save summary CSV
    df_summary = pd.DataFrame(summary_rows)
    csv_out = DOC_DIR / "diffusion_statistics_results.csv"
    df_summary.to_csv(csv_out, index=False, float_format="%.6e")
    print(f"Summary CSV saved to {csv_out}")

    # Save per-track CSV
    if not df_tamsd.empty:
        tamsd_csv = DOC_DIR / "diffusion_statistics_tamsd_per_track.csv"
        df_tamsd.to_csv(tamsd_csv, index=False, float_format="%.6e")
        print(f"Per-track CSV saved to {tamsd_csv}")

    # Save per-file CSV
    if not df_eamsd.empty:
        eamsd_csv = DOC_DIR / "diffusion_statistics_eamsd_per_file.csv"
        df_eamsd.to_csv(eamsd_csv, index=False, float_format="%.6e")
        print(f"Per-file eaMSD CSV saved to {eamsd_csv}")

    # Save ensemble-taMSD per-file CSV
    if not df_ens_tamsd.empty:
        ens_tamsd_csv = DOC_DIR / "diffusion_statistics_ens_tamsd_per_file.csv"
        df_ens_tamsd.to_csv(ens_tamsd_csv, index=False, float_format="%.6e")
        print(f"Per-file <taMSD> CSV saved to {ens_tamsd_csv}")

    # ── Final summary table ─────────────────────────────────────
    print(f"\n{'='*120}")
    print("FINAL SUMMARY")
    print(f"{'='*120}")
    print(f"D_theory = {D_theory:.6e} um^2/s\n")

    header = (f"{'Model':<12} {'Frac':>5} | "
              f"{'<D>_taMSD':>11} {'+/-SEM':>11} {'N':>5} {'p-val':>10} | "
              f"{'<D>_eaMSD':>11} {'+/-SEM':>11} {'N':>5} {'p-val':>10} | "
              f"{'<D>_<taMSD>':>11} {'+/-SEM':>11} {'N':>5} {'p-val':>10}")
    print(header)
    print("-" * len(header))
    for r in summary_rows:
        pct = f"{int(r['fraction']*100)}%"
        print(
            f"{r['model']:<12} {pct:>5} | "
            f"{r['D_tamsd_mean']:>11.4e} {r['D_tamsd_sem']:>11.4e} "
            f"{r['N_tracks']:>5} {r['p_val_tamsd']:>10.4e} | "
            f"{r['D_eamsd_mean']:>11.4e} {r['D_eamsd_sem']:>11.4e} "
            f"{r['N_files_ea']:>5} {r['p_val_eamsd']:>10.4e} | "
            f"{r['D_ens_tamsd_mean']:>11.4e} {r['D_ens_tamsd_sem']:>11.4e} "
            f"{r['N_files_et']:>5} {r['p_val_ens_tamsd']:>10.4e}"
        )
    print(f"{'='*120}")
    print("Done.")


if __name__ == "__main__":
    main()
