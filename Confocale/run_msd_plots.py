"""
Batch EA-MSD and TA-MSD plots for both datasets.

For every CSV file in Data/31_10_no_anomalous/ and Data/14_11_anomalous/ and
for each of the four lag fractions [0.10, 0.25, 0.50, 1.0]:
  - Generates one EA-MSD plot with SEM error bars  (eamsd_plot.py)
  - Generates one TA-MSD plot with SEM error bars  (tamsd_plot.py, first/default track)

Output folder structure:
    Results/
        no_anomalous/
            eamsd/     — SVGs named  <stem>_eamsd_f<pct>.svg
            tamsd/     — SVGs named  <stem>_tamsd_f<pct>.svg
        anomalous/
            eamsd/
            tamsd/

Usage:
    python run_msd_plots.py
"""

import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent

DATASETS = {
    "no_anomalous": SCRIPT_DIR / "Data" / "31_10_no_anomalous",
    "anomalous":    SCRIPT_DIR / "Data" / "14_11_anomalous",
}

LAG_FRACTIONS = [0.10, 0.25, 0.50, 1.0]


def run(cmd: list, label: str) -> subprocess.CompletedProcess:
    print(f"  [{label}]", end=" ", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True, check=False, cwd=SCRIPT_DIR)
    if result.returncode != 0:
        print(f"FAILED\n    {result.stderr.strip()[:400]}")
    else:
        print("OK")
    return result


def process_dataset(label: str, data_dir: Path) -> None:
    eamsd_dir = SCRIPT_DIR / "Results" / label / "eamsd"
    tamsd_dir  = SCRIPT_DIR / "Results" / label / "tamsd"
    eamsd_dir.mkdir(parents=True, exist_ok=True)
    tamsd_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(data_dir.glob("*.csv"))
    if not files:
        print(f"  No CSV files found in {data_dir}")
        return

    print(f"\n{'='*60}")
    print(f"Dataset: {label}  ({len(files)} files × {len(LAG_FRACTIONS)} fractions)")
    print(f"{'='*60}")

    for i, csv_file in enumerate(files, 1):
        stem = csv_file.stem
        print(f"  [{i:02d}/{len(files):02d}] {csv_file.name}")

        for frac in LAG_FRACTIONS:
            pct = int(frac * 100)
            frac_label = f"f{pct:03d}"

            # ---- EA-MSD ------------------------------------------------
            eamsd_out = eamsd_dir / f"{stem}_eamsd_{frac_label}.svg"
            run(
                [sys.executable, "eamsd_plot.py", str(csv_file),
                 "--max-lag-fraction", str(frac),
                 "--output", str(eamsd_out)],
                f"eamsd {pct:3d}%",
            )

            # ---- TA-MSD ------------------------------------------------
            tamsd_out = tamsd_dir / f"{stem}_tamsd_{frac_label}.svg"
            run(
                [sys.executable, "tamsd_plot.py", str(csv_file),
                 "--max-lag-fraction", str(frac),
                 "--output", str(tamsd_out)],
                f"tamsd {pct:3d}%",
            )


def main() -> None:
    for label, data_dir in DATASETS.items():
        process_dataset(label, data_dir)
    print("\nAll done!")


if __name__ == "__main__":
    main()
