"""
Batch MSD analysis — plots and fits for all datasets.

For every CSV in each dataset:
  - EA-MSD and TA-MSD plots at multiple lag fractions
  - Model-specific fits with SEM error bars

Output:
    Results/<label>/eamsd/   — EA-MSD SVGs
    Results/<label>/tamsd/   — TA-MSD SVGs
    Results/<label>/fits/    — Fit SVGs (linear_offset_fits/ and nonlinear_fits/)
    Docu/                    — Summary CSV + Markdown tables

Usage:
    python run_all.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from check_ergodicity import compute_ensemble_tamsd
from data_reader import read_trajectories_from_csv
from msd_analyzer import (
    calculate_ensemble_msd,
    compute_ensemble_drift,
)
from msd_fitting import (
    fit_msd_linear,
    fit_msd_linear_offset,
    fit_msd_nonlinear,
    analyze_velocities,
    fit_msd_anomalous,
    fit_msd_anomalous_drift,
    fit_msd_anomalous_offset,
)
from plot_msd import plot_and_save

SCRIPT_DIR = Path(__file__).parent
DOC_DIR = SCRIPT_DIR / "Docu"

DATASETS = {
    "no_anomalous": SCRIPT_DIR / "Data" / "31_10_no_anomalous",
    "anomalous":    SCRIPT_DIR / "Data" / "14_11_anomalous",
}

# Which fit models to run per dataset
DATASET_MODELS = {
    "no_anomalous": ["linear_offset", "nonlinear"],
    "anomalous":    ["anomalous_offset"],
}

LAG_FRACTIONS = [0.10, 0.25]


# ---------------------------------------------------------------------------
# Fit-plot helper (same style as fit_msd.py)
# ---------------------------------------------------------------------------

def _msd_ylabel(msd_kind: str) -> str:
    labels = {
        "eaMSD": r"eaMSD [$\mu$m$^2$]",
        "taMSD": r"taMSD [$\mu$m$^2$]",
        "<taMSD>": r"$\langle$taMSD$\rangle$ [$\mu$m$^2$]",
    }
    return labels.get(msd_kind, r"MSD [$\mu$m$^2$]")


def _msd_data_label(msd_kind: str) -> str:
    labels = {
        "eaMSD": "eaMSD Data",
        "taMSD": "taMSD Data",
        "<taMSD>": "taMSD Data",
    }
    return labels.get(msd_kind, "MSD Data")


def _save_fit_plot(tau_fit, msd_fit, msd_predicted, textstr, fit_label,
                   output_path, msd_sigma=None, msd_kind="MSD",
                   text_below_legend=False):
    fig, ax = plt.subplots(figsize=(8, 6))
    data_label = _msd_data_label(msd_kind)
    if msd_sigma is not None:
        ax.errorbar(tau_fit, msd_fit, yerr=msd_sigma, fmt="o", color="C0",
                    markersize=8, alpha=0.7, capsize=4, capthick=1.2,
                    elinewidth=1.2, label=data_label, zorder=2)
    else:
        ax.plot(tau_fit, msd_fit, "o", color="C0", markersize=8, alpha=0.7,
                label=data_label, zorder=2)
    ax.plot(tau_fit, msd_predicted, "-", color="C3", linewidth=2.5,
            label=fit_label, zorder=3)
    ax.set_xlabel(r"Time Lag $\tau$ [s]", fontsize=12)
    ax.set_ylabel(_msd_ylabel(msd_kind), fontsize=12)
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(loc="upper left", fontsize=10, framealpha=0.9)
    props = dict(boxstyle="round", facecolor="white", alpha=0.95,
                 edgecolor="black", linewidth=1.2)
    ax.text(0.02, 0.87, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment="top", horizontalalignment="left", bbox=props)
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Processing
# ---------------------------------------------------------------------------

def run_msd_plots(csv_file: Path, label: str) -> None:
    """Generate EA-MSD and TA-MSD plots at all lag fractions."""
    trajectories = read_trajectories_from_csv(str(csv_file))
    if not trajectories:
        print(f"    No trajectories in {csv_file.name}")
        return

    stem = csv_file.stem
    eamsd_dir = SCRIPT_DIR / "Results" / label / "eamsd"
    tamsd_dir = SCRIPT_DIR / "Results" / label / "tamsd"
    eamsd_dir.mkdir(parents=True, exist_ok=True)
    tamsd_dir.mkdir(parents=True, exist_ok=True)

    # Use drift correction for anomalous datasets
    use_drift_corr = (label == "anomalous")

    # For the anomalous dataset also save raw (uncorrected) plots for comparison
    if use_drift_corr:
        eamsd_raw_dir = SCRIPT_DIR / "Results" / label / "eamsd_raw"
        tamsd_raw_dir = SCRIPT_DIR / "Results" / label / "tamsd_raw"
        eamsd_raw_dir.mkdir(parents=True, exist_ok=True)
        tamsd_raw_dir.mkdir(parents=True, exist_ok=True)
    else:
        eamsd_raw_dir = None
        tamsd_raw_dir = None

    global_dt = None  # inferred lazily by compute_ensemble_tamsd

    for frac in LAG_FRACTIONS:
        pct = int(frac * 100)
        tag = f"f{pct:03d}"

        # EA-MSD
        try:
            ea = calculate_ensemble_msd(
                trajectories, max_lag_fraction=frac,
                drift_corrected=use_drift_corr,
            )
            if ea.tau.size > 0:
                out = eamsd_dir / f"{stem}_eamsd_{tag}.svg"
                plot_and_save(ea.tau, ea.msd, out, msd_sem=ea.msd_sem,
                              label="EA-MSD", color="C0")
                print(f"    eamsd {pct:3d}% OK")
        except Exception as e:
            print(f"    eamsd {pct:3d}% FAILED: {e}")

        if eamsd_raw_dir is not None:
            try:
                ea_raw = calculate_ensemble_msd(
                    trajectories, max_lag_fraction=frac,
                    drift_corrected=False,
                )
                if ea_raw.tau.size > 0:
                    out = eamsd_raw_dir / f"{stem}_eamsd_{tag}.svg"
                    plot_and_save(ea_raw.tau, ea_raw.msd, out, msd_sem=ea_raw.msd_sem,
                                  label="EA-MSD (raw)", color="C0")
                    print(f"    eamsd_raw {pct:3d}% OK")
            except Exception as e:
                print(f"    eamsd_raw {pct:3d}% FAILED: {e}")

        # TA-MSD — ensemble average over all tracks (SEM = std / sqrt(N_tracks))
        try:
            tau_ta, msd_ta, sem_ta, _ = compute_ensemble_tamsd(
                trajectories, max_lag_fraction=frac,
                drift_corrected=use_drift_corr,
            )
            if tau_ta.size > 0:
                out = tamsd_dir / f"{stem}_tamsd_{tag}.svg"
                plot_and_save(tau_ta, msd_ta, out, msd_sem=sem_ta,
                              label="⟨TA-MSD⟩", color="C1")
                print(f"    tamsd {pct:3d}% OK")
        except Exception as e:
            print(f"    tamsd {pct:3d}% FAILED: {e}")

        if tamsd_raw_dir is not None:
            try:
                tau_ta_raw, msd_ta_raw, sem_ta_raw, _ = compute_ensemble_tamsd(
                    trajectories, max_lag_fraction=frac,
                    drift_corrected=False,
                )
                if tau_ta_raw.size > 0:
                    out = tamsd_raw_dir / f"{stem}_tamsd_{tag}.svg"
                    plot_and_save(tau_ta_raw, msd_ta_raw, out, msd_sem=sem_ta_raw,
                                  label="⟨TA-MSD⟩ (raw)", color="C1")
                    print(f"    tamsd_raw {pct:3d}% OK")
            except Exception as e:
                print(f"    tamsd_raw {pct:3d}% FAILED: {e}")


def run_fits(csv_file: Path, label: str) -> list[dict]:
    """Run all configured fits for a dataset and return summary rows."""
    models = DATASET_MODELS.get(label, [])
    trajectories = read_trajectories_from_csv(str(csv_file))
    if not trajectories:
        print(f"    No trajectories in {csv_file.name}")
        return []

    # Use drift-corrected eaMSD for anomalous datasets
    use_drift_corr = (label == "anomalous")
    msd_result = calculate_ensemble_msd(trajectories, drift_corrected=use_drift_corr)
    if msd_result.tau.size == 0:
        print(f"    Empty MSD for {csv_file.name}")
        return []

    stem = csv_file.stem
    linear_dir = SCRIPT_DIR / "Results" / label / "linear_fits"
    linear_offset_dir = SCRIPT_DIR / "Results" / label / "linear_offset_fits"
    nonlinear_dir = SCRIPT_DIR / "Results" / label / "nonlinear_fits"
    linear_dir.mkdir(parents=True, exist_ok=True)
    linear_offset_dir.mkdir(parents=True, exist_ok=True)
    nonlinear_dir.mkdir(parents=True, exist_ok=True)

    summary = []

    for model in models:
        try:
            row = _run_single_fit(
                model, msd_result, trajectories, stem,
                linear_dir, linear_offset_dir, nonlinear_dir, label, csv_file.name,
            )
            summary.append(row)
            print(f"    {model} OK")
        except Exception as e:
            print(f"    {model} FAILED: {e}")

    return summary


def _run_single_fit(model, msd_result, trajectories, stem,
                    linear_dir, linear_offset_dir, nonlinear_dir,
                    dataset_label, filename):
    """Run one fit model and return a summary dict."""
    tau, msd = msd_result.tau, msd_result.msd
    n_max, dt = msd_result.n_max, msd_result.dt
    sigma = msd_result.msd_sem

    if model == "linear":
        fit = fit_msd_linear(tau, msd, n_max, dt, msd_sigma=sigma)
        txt = (r"$D = (%.2e \pm %.1e)\ \mu m^2/s$" % (fit.D, fit.D_error)
               + "\n" + r"$\chi^2_\nu = %.4f$" % fit.chi_squared_red)
        out = linear_dir / f"{stem}_linear_fit.svg"
        _save_fit_plot(fit.tau_fit, fit.msd_fit, fit.msd_predicted, txt,
                       r"Linear: MSD = 4D$\tau$", out, fit.msd_sigma_fit,
                       msd_kind="eaMSD")
        return dict(Dataset=dataset_label, File=filename,
                    Model="Linear (4Dτ)",
                    D=f"{fit.D:.4e}", D_err=f"{fit.D_error:.2e}",
                    alpha="N/A", alpha_err="N/A",
                    v="N/A", v_err="N/A",
                    chi2_red=f"{fit.chi_squared_red:.4f}")

    elif model == "linear_offset":
        fit = fit_msd_linear_offset(tau, msd, n_max, dt, msd_sigma=sigma)
        txt = "\n".join([
            r"$D = (%.2e \pm %.1e)\ \mu m^2/s$" % (fit.D, fit.D_error),
            r"$c = (%.2e \pm %.1e)\ \mu m^2$" % (fit.offset, fit.offset_error),
            r"$\chi^2_\nu = %.4f$" % fit.chi_squared_red,
        ])
        out = linear_offset_dir / f"{stem}_linear_offset_fit.svg"
        _save_fit_plot(
            fit.tau_fit, fit.msd_fit, fit.msd_predicted, txt,
            r"Linear+offset: MSD = 4D$\tau$ + $c$", out,
            fit.msd_sigma_fit, msd_kind="eaMSD", text_below_legend=True,
        )
        return dict(Dataset=dataset_label, File=filename,
                    Model="Linear+Offset (4Dτ + c)",
                    D=f"{fit.D:.4e}", D_err=f"{fit.D_error:.2e}",
                    alpha="N/A", alpha_err="N/A",
                    v="N/A", v_err="N/A",
                    chi2_red=f"{fit.chi_squared_red:.4f}")

    elif model == "nonlinear":
        vstats = analyze_velocities(trajectories)
        fit = fit_msd_nonlinear(tau, msd, n_max, dt,
                                velocity_stats=vstats, msd_sigma=sigma)
        txt = "\n".join([
            r"$D = (%.2e \pm %.1e)\ \mu m^2/s$" % (fit.D, fit.D_error),
            r"$v = (%.2e \pm %.1e)\ \mu m/s$" % (fit.v, fit.v_error),
            r"$\chi^2_\nu = %.4f$" % fit.chi_squared_red,
        ])
        out = nonlinear_dir / f"{stem}_nonlinear_fit.svg"
        _save_fit_plot(fit.tau_fit, fit.msd_fit, fit.msd_predicted, txt,
                       r"Nonlinear: MSD = 4D$\tau$ + $v^2\tau^2$", out,
                       fit.msd_sigma_fit, msd_kind="eaMSD")
        return dict(Dataset=dataset_label, File=filename,
                    Model="Nonlinear (4Dτ + v²τ²)",
                    D=f"{fit.D:.4e}", D_err=f"{fit.D_error:.2e}",
                    alpha="N/A", alpha_err="N/A",
                    v=f"{fit.v:.4e}", v_err=f"{fit.v_error:.2e}",
                    chi2_red=f"{fit.chi_squared_red:.4f}")

    elif model == "anomalous":
        fit = fit_msd_anomalous(tau, msd, n_max, dt, msd_sigma=sigma)
        txt = "\n".join([
            r"$D_\alpha = (%.2e \pm %.1e)\ \mu m^2/s^\alpha$" % (fit.D_alpha, fit.D_alpha_error),
            r"$\alpha = %.4f \pm %.4f$" % (fit.alpha, fit.alpha_error),
            r"$\chi^2_\nu = %.4f$" % fit.chi_squared_red,
        ])
        out = linear_dir / f"{stem}_anomalous_fit.svg"
        _save_fit_plot(fit.tau_fit, fit.msd_fit, fit.msd_predicted, txt,
                       r"Anomalous: MSD = 4$D_\alpha \tau^\alpha$", out,
                       msd_kind="eaMSD")
        return dict(Dataset=dataset_label, File=filename,
                    Model="Anomalous (4D_α τ^α)",
                    D=f"{fit.D_alpha:.4e}", D_err=f"{fit.D_alpha_error:.2e}",
                    alpha=f"{fit.alpha:.4f}", alpha_err=f"{fit.alpha_error:.4f}",
                    v="N/A", v_err="N/A",
                    chi2_red=f"{fit.chi_squared_red:.4f}")

    elif model == "anomalous_offset":
        fit = fit_msd_anomalous_offset(tau, msd, n_max, dt, msd_sigma=sigma)
        txt = "\n".join([
            r"$D_\alpha = (%.2e \pm %.1e)\ \mu m^2/s^\alpha$" % (fit.D_alpha, fit.D_alpha_error),
            r"$\alpha = %.4f \pm %.4f$" % (fit.alpha, fit.alpha_error),
            r"$c = (%.2e \pm %.1e)\ \mu m^2$" % (fit.offset, fit.offset_error),
            r"$\chi^2_\nu = %.4f$" % fit.chi_squared_red,
        ])
        out = linear_dir / f"{stem}_anomalous_offset_fit.svg"
        _save_fit_plot(fit.tau_fit, fit.msd_fit, fit.msd_predicted, txt,
                       r"Drift-corr. anom.+offset: 4$D_\alpha \tau^\alpha$ + $c$", out,
                       fit.msd_sigma_fit, msd_kind="eaMSD")
        return dict(Dataset=dataset_label, File=filename,
                    Model="AnomalousOffset (4D_α τ^α + c)",
                    D=f"{fit.D_alpha:.4e}", D_err=f"{fit.D_alpha_error:.2e}",
                    alpha=f"{fit.alpha:.4f}", alpha_err=f"{fit.alpha_error:.4f}",
                    v="N/A", v_err="N/A",
                    chi2_red=f"{fit.chi_squared_red:.4f}")

    raise ValueError(f"Unknown model: {model}")


def write_summary(summary: list[dict], dataset_label: str) -> None:
    """Save CSV and Markdown summary tables."""
    DOC_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(summary)

    csv_out = DOC_DIR / f"fits_{dataset_label}_results.csv"
    df.to_csv(csv_out, index=False)

    md_out = DOC_DIR / f"fits_{dataset_label}_results.md"
    title = dataset_label.replace("_", " ").title()
    with open(md_out, "w", encoding="utf-8") as f:
        f.write(f"# MSD Fit Results — {title}\n\n")
        f.write(f"**Date:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("| File | Model | D | D_err | alpha | alpha_err | v | v_err | chi2_red |\n")
        f.write("|---|---|---|---|---|---|---|---|---|\n")
        for _, row in df.iterrows():
            f.write(f"| {row['File']} | {row['Model']} | {row['D']} | {row['D_err']} "
                    f"| {row['alpha']} | {row['alpha_err']} "
                    f"| {row['v']} | {row['v_err']} | {row['chi2_red']} |\n")

    print(f"  Summary: {csv_out.name}, {md_out.name}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    all_summary: list[dict] = []

    for label, data_dir in DATASETS.items():
        files = sorted(data_dir.rglob("*.csv"))
        if not files:
            print(f"No CSV files in {data_dir}")
            continue

        print(f"\n{'=' * 60}")
        print(f"Dataset: {label}  ({len(files)} files)")
        print(f"{'=' * 60}")

        dataset_summary: list[dict] = []

        for i, csv_file in enumerate(files, 1):
            print(f"\n  [{i:02d}/{len(files):02d}] {csv_file.name}")
            run_msd_plots(csv_file, label)
            rows = run_fits(csv_file, label)
            dataset_summary.extend(rows)

        if dataset_summary:
            write_summary(dataset_summary, label)
            all_summary.extend(dataset_summary)

    if all_summary:
        write_summary(all_summary, "all_datasets")

    print("\nAll done!")


if __name__ == "__main__":
    main()
