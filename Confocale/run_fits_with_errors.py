"""
Batch fitting script — runs MSD fits on ALL CSV files
in both Data/31_10_no_anomalous/ and Data/14_11_anomalous/.

For no_anomalous: standard models
    linear_fits/    — MSD = 4Dτ           (plot_msd_fit.py)
    nonlinear_fits/ — MSD = 4Dτ + v²τ²   (msd_fitting_nonlinear.py)

For anomalous: anomalous models
    linear_fits/    — MSD = 4D_α τ^α              (msd_fitting_anomalous.py)
    nonlinear_fits/ — MSD = 4D_α τ^α + v²τ²       (msd_fitting_anomalous_drift.py)

Output folder structure:
    Results/
        no_anomalous/
            linear_fits/
            nonlinear_fits/
        anomalous/
            linear_fits/
            nonlinear_fits/

Summary tables saved to Docu/:
    fits_no_anomalous_results.csv / .md
    fits_anomalous_results.csv   / .md

Usage:
    python run_fits_with_errors.py
"""

import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).parent
DOC_DIR = SCRIPT_DIR / "Docu"

DATASETS = {
    "no_anomalous": SCRIPT_DIR / "Data" / "31_10_no_anomalous",
    "anomalous":    SCRIPT_DIR / "Data" / "14_11_anomalous",
}

# Datasets that require the anomalous diffusion models (MSD = 4D_α τ^α)
ANOMALOUS_DATASETS = {"anomalous"}


def run(cmd: list, label: str) -> subprocess.CompletedProcess:
    print(f"  [{label}]", end=" ", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True, check=False, cwd=SCRIPT_DIR)
    if result.returncode != 0:
        print(f"FAILED\n    {result.stderr.strip()[:300]}")
    else:
        print("OK")
    return result


def process_dataset(dataset_label: str, data_dir: Path) -> list:
    """Run fits on every CSV in data_dir.

    For standard datasets: linear + nonlinear (drift) fits.
    For anomalous datasets: anomalous (τ^α) + anomalous+drift fits.

    Plots are saved under Results/<dataset_label>/linear_fits/ and
    Results/<dataset_label>/nonlinear_fits/.

    Returns a list of summary-row dicts.
    """
    is_anomalous = dataset_label in ANOMALOUS_DATASETS

    linear_dir   = SCRIPT_DIR / "Results" / dataset_label / "linear_fits"
    nonlinear_dir = SCRIPT_DIR / "Results" / dataset_label / "nonlinear_fits"
    linear_dir.mkdir(parents=True, exist_ok=True)
    nonlinear_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(data_dir.glob("*.csv"))
    if not files:
        print(f"  No CSV files found in {data_dir}")
        return []

    linear_model   = "Anomalous (4D_α τ^α)"       if is_anomalous else "Linear (4Dτ)"
    nonlinear_model = "AnomalousDrift (4D_α τ^α + v²τ²)" if is_anomalous else "Nonlinear (4Dτ + v²τ²)"

    print(f"\n{'='*60}")
    print(f"Dataset: {dataset_label}  ({len(files)} files)")
    print(f"  Data dir  : {data_dir}")
    print(f"  Fit model : {'ANOMALOUS' if is_anomalous else 'STANDARD'}")
    print(f"  Linear dir: {linear_dir}")
    print(f"  Nonlin dir: {nonlinear_dir}")
    print(f"{'='*60}")

    summary = []

    for i, csv_file in enumerate(files, 1):
        stem = csv_file.stem
        print(f"  [{i:02d}/{len(files):02d}] {csv_file.name}")

        if is_anomalous:
            # ---- Anomalous fit: MSD = 4D_α τ^α ----------------------------
            res_lin = run(
                [sys.executable, "msd_fitting_anomalous.py", str(csv_file),
                 "--output-dir", str(linear_dir)],
                "anomalous fit",
            )

            d_lin = d_lin_err = alpha_lin = alpha_lin_err = chi_lin = "N/A"
            if res_lin and res_lin.stdout:
                md = re.search(
                    r"Diffusion Coeff \(D_alpha\):\s+\(([\d\.e\+\-]+) \+/- ([\d\.e\+\-]+)\)",
                    res_lin.stdout,
                )
                ma = re.search(
                    r"Anomalous exponent \(alpha\):\s+\(([\d\.e\+\-]+) \+/- ([\d\.e\+\-]+)\)",
                    res_lin.stdout,
                )
                mc = re.search(r"chi\^2_red:\s+([\d\.e\+\-]+|nan)", res_lin.stdout)
                if md:
                    d_lin, d_lin_err = md.group(1), md.group(2)
                if ma:
                    alpha_lin, alpha_lin_err = ma.group(1), ma.group(2)
                if mc:
                    chi_lin = mc.group(1)

            summary.append({
                "Dataset": dataset_label, "File": csv_file.name, "Model": linear_model,
                "D": d_lin, "D_err": d_lin_err,
                "alpha": alpha_lin, "alpha_err": alpha_lin_err,
                "v": "N/A", "v_err": "N/A",
                "chi2_red": chi_lin,
            })

            # ---- Anomalous+Drift fit: MSD = 4D_α τ^α + v²τ² ---------------
            res_nl = run(
                [sys.executable, "msd_fitting_anomalous_drift.py", str(csv_file),
                 "--output-dir", str(nonlinear_dir)],
                "anomalous+drift fit",
            )

            d_nl = d_nl_err = alpha_nl = alpha_nl_err = v_nl = v_nl_err = chi_nl = "N/A"
            if res_nl and res_nl.stdout:
                md = re.search(
                    r"Diffusion Coeff \(D_alpha\):\s+\(([\d\.e\+\-]+) \+/- ([\d\.e\+\-]+)\)",
                    res_nl.stdout,
                )
                ma = re.search(
                    r"Anomalous exponent \(alpha\):\s+\(([\d\.e\+\-]+) \+/- ([\d\.e\+\-]+)\)",
                    res_nl.stdout,
                )
                mv = re.search(
                    r"Drift Velocity \(v\):\s+\(([\d\.e\+\-]+) \+/- ([\d\.e\+\-]+)\)",
                    res_nl.stdout,
                )
                mc = re.search(r"chi\^2_red:\s+([\d\.e\+\-]+|nan)", res_nl.stdout)
                if md:
                    d_nl, d_nl_err = md.group(1), md.group(2)
                if ma:
                    alpha_nl, alpha_nl_err = ma.group(1), ma.group(2)
                if mv:
                    v_nl, v_nl_err = mv.group(1), mv.group(2)
                if mc:
                    chi_nl = mc.group(1)

            summary.append({
                "Dataset": dataset_label, "File": csv_file.name, "Model": nonlinear_model,
                "D": d_nl, "D_err": d_nl_err,
                "alpha": alpha_nl, "alpha_err": alpha_nl_err,
                "v": v_nl, "v_err": v_nl_err,
                "chi2_red": chi_nl,
            })

        else:
            # ---- Standard linear fit: MSD = 4Dτ ----------------------------
            linear_out = linear_dir / f"{stem}_linear_fit.svg"
            res_lin = run(
                [sys.executable, "plot_msd_fit.py", str(csv_file), "--output", str(linear_out)],
                "linear fit",
            )

            d_lin = d_lin_err = chi_lin = "N/A"
            if res_lin and res_lin.stdout:
                m = re.search(
                    r"Diffusion Coefficient \(D\): \(([\d\.e\+\-]+) \+/- ([\d\.e\+\-]+)\)",
                    res_lin.stdout,
                )
                c = re.search(r"chi\^2_red:\s+([\d\.e\+\-]+|nan)", res_lin.stdout)
                if m:
                    d_lin, d_lin_err = m.group(1), m.group(2)
                if c:
                    chi_lin = c.group(1)

            summary.append({
                "Dataset": dataset_label, "File": csv_file.name, "Model": linear_model,
                "D": d_lin, "D_err": d_lin_err,
                "alpha": "N/A", "alpha_err": "N/A",
                "v": "N/A", "v_err": "N/A",
                "chi2_red": chi_lin,
            })

            # ---- Standard nonlinear fit: MSD = 4Dτ + v²τ² -----------------
            res_nl = run(
                [sys.executable, "msd_fitting_nonlinear.py", str(csv_file),
                 "--output-dir", str(nonlinear_dir)],
                "nonlinear fit",
            )

            d_nl = d_nl_err = v_nl = v_nl_err = chi_nl = "N/A"
            if res_nl and res_nl.stdout:
                md = re.search(
                    r"Diffusion Coeff \(D\):\s+\(([\d\.e\+\-]+) \+/- ([\d\.e\+\-]+)\)",
                    res_nl.stdout,
                )
                mv = re.search(
                    r"Drift Velocity \(v\):\s+\(([\d\.e\+\-]+) \+/- ([\d\.e\+\-]+)\)",
                    res_nl.stdout,
                )
                mc = re.search(r"chi\^2_red:\s+([\d\.e\+\-]+|nan)", res_nl.stdout)
                if md:
                    d_nl, d_nl_err = md.group(1), md.group(2)
                if mv:
                    v_nl, v_nl_err = mv.group(1), mv.group(2)
                if mc:
                    chi_nl = mc.group(1)

            summary.append({
                "Dataset": dataset_label, "File": csv_file.name, "Model": nonlinear_model,
                "D": d_nl, "D_err": d_nl_err,
                "alpha": "N/A", "alpha_err": "N/A",
                "v": v_nl, "v_err": v_nl_err,
                "chi2_red": chi_nl,
            })

    return summary


def write_summary(summary: list, dataset_label: str) -> None:
    """Save CSV and Markdown summary tables for one dataset."""
    DOC_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(summary)

    # Ensure alpha columns exist (may be absent for older standard-only runs)
    for col in ("alpha", "alpha_err"):
        if col not in df.columns:
            df[col] = "N/A"

    csv_out = DOC_DIR / f"fits_{dataset_label}_results.csv"
    df.to_csv(csv_out, index=False)

    md_out = DOC_DIR / f"fits_{dataset_label}_results.md"
    title = dataset_label.replace("_", " ").title()
    with open(md_out, "w", encoding="utf-8") as f:
        f.write(f"# MSD Fit Results — {title} (SEM-weighted, chi2_red metric)\n\n")
        f.write(f"**Date:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("| File | Model | D (um2/s or um2/s^α) | D_err | alpha | alpha_err | v (um/s) | v_err | chi2_red |\n")
        f.write("|---|---|---|---|---|---|---|---|---|\n")
        for _, row in df.iterrows():
            f.write(
                f"| {row['File']} | {row['Model']} | {row['D']} | {row['D_err']} "
                f"| {row.get('alpha', 'N/A')} | {row.get('alpha_err', 'N/A')} "
                f"| {row['v']} | {row['v_err']} | {row['chi2_red']} |\n"
            )

    print(f"  Summary CSV : {csv_out}")
    print(f"  Summary MD  : {md_out}")


def main() -> None:
    all_summary = []

    for label, data_dir in DATASETS.items():
        summary = process_dataset(label, data_dir)
        if summary:
            write_summary(summary, label)
            all_summary.extend(summary)

    # Write a single combined summary as well
    if all_summary:
        write_summary(all_summary, "all_datasets")

    print("\nAll done!")


if __name__ == "__main__":
    main()
