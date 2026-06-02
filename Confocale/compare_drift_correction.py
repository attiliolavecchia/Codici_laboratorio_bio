"""
Drift-correction impact analysis — raw vs drift-corrected eaMSD and ⟨taMSD⟩.

For each CSV file in the anomalous dataset we compute both MSD estimators
(eaMSD and ensemble-averaged taMSD) **with** and **without** drift correction
and overlay them on the same log-log plot.  This lets us quantify how much the
drift correction actually changed the curves.

Drift-correction methods used
------------------------------
  eaMSD : variance method  — MSD_corr(τ) = Var(Δx_n) + Var(Δy_n)
              (subtracts the squared ensemble-mean displacement |⟨Δr⟩|²)
  ⟨taMSD⟩ : Crocker-Grier ensemble common-mode drift subtraction
              (Crocker & Grier 1996, J. Colloid Interface Sci. 179:298)

Output
------
    Results/anomalous/drift_correction_comparison/
        eaMSD_drift_cmp_<group>_<file>.svg   — per-file eaMSD overlay
        taMSD_drift_cmp_<group>_<file>.svg   — per-file ⟨taMSD⟩ overlay
        eaMSD_drift_cmp_<group>_median.svg   — group-level median overlay
        taMSD_drift_cmp_<group>_median.svg   — group-level median overlay

Usage:
    python compare_drift_correction.py                  # both groups
    python compare_drift_correction.py --group glic50
    python compare_drift_correction.py --group glic200
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from check_ergodicity import compute_ensemble_tamsd
from data_reader import read_trajectories_from_csv, estimate_global_time_step
from msd_analyzer import calculate_ensemble_msd, compute_ensemble_drift

# ── Configuration ──────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "Data" / "14_11_anomalous"
OUT_DIR = SCRIPT_DIR / "Results" / "anomalous" / "drift_correction_comparison"

GLIC50_SUBDIR = "PEG_15_Glicerolo_50_H2O_35"
GLIC200_SUBDIR = "PEG_18_Glicerolo_40_H2O_42"

MSD_LAG_FRACTION = 0.50   # compute MSD up to 50% of the longest trajectory
FIT_FRACTIONS = [0.10, 0.25]  # fractions of lag points shown in each plot

# Visual colours: raw = warm, corrected = cool
COL_RAW = "#e07b54"        # orange-red  → raw (no drift correction)
COL_CORR = "#3b82f6"       # blue        → drift-corrected


# ── Plot helpers ───────────────────────────────────────────────────────

def _clean_axes(ax):
    ax.grid(True, linestyle=":", alpha=0.6)


def _plot_overlay(tau_raw, msd_raw, sem_raw,
                  tau_corr, msd_corr, sem_corr,
                  ylabel, out_path, n_show=None):
    """Save a linear-axes overlay plot (raw vs drift-corrected),
    restricted to the first *n_show* lag points."""
    if n_show is not None:
        n_show = max(4, int(n_show))
        tau_raw  = tau_raw[:n_show]
        msd_raw  = msd_raw[:n_show]
        sem_raw  = sem_raw[:n_show]
        tau_corr = tau_corr[:n_show]
        msd_corr = msd_corr[:n_show]
        sem_corr = sem_corr[:n_show]

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    _clean_axes(ax)

    ax.errorbar(tau_raw, msd_raw, yerr=sem_raw,
                fmt="o-", color=COL_RAW,
                capsize=4, capthick=1.2, elinewidth=1.2,
                label="Raw (no drift correction)")
    ax.errorbar(tau_corr, msd_corr, yerr=sem_corr,
                fmt="s-", color=COL_CORR,
                capsize=4, capthick=1.2, elinewidth=1.2,
                label="Drift-corrected")

    ax.set_xlabel(r"Time Lag ($\tau$) [seconds]")
    ax.set_ylabel(ylabel)
    ax.legend()
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


# ── Per-file computation ───────────────────────────────────────────────

def analyse_file(csv_path: Path):
    """Return eaMSD and ⟨taMSD⟩ curves raw vs corrected for one CSV file."""
    trajectories = read_trajectories_from_csv(str(csv_path))
    if not trajectories:
        return None

    global_dt = estimate_global_time_step(trajectories)

    # ── eaMSD ─────────────────────────────────────────────────────────
    ea_raw = calculate_ensemble_msd(
        trajectories, max_lag_fraction=MSD_LAG_FRACTION,
        global_dt=global_dt, drift_corrected=False,
    )
    ea_corr = calculate_ensemble_msd(
        trajectories, max_lag_fraction=MSD_LAG_FRACTION,
        global_dt=global_dt, drift_corrected=True,
    )
    if ea_raw.tau.size < 4 or ea_corr.tau.size < 4:
        return None

    # ── ⟨taMSD⟩ ──────────────────────────────────────────────────────
    tau_ta_raw, msd_ta_raw, sem_ta_raw, _ = compute_ensemble_tamsd(
        trajectories, max_lag_fraction=MSD_LAG_FRACTION,
        global_dt=global_dt, drift_corrected=False,
    )
    tau_ta_corr, msd_ta_corr, sem_ta_corr, _ = compute_ensemble_tamsd(
        trajectories, max_lag_fraction=MSD_LAG_FRACTION,
        global_dt=global_dt, drift_corrected=True,
    )

    # Common τ-length for eaMSD pair
    n_ea = min(ea_raw.tau.size, ea_corr.tau.size)
    # Common τ-length for taMSD pair
    n_ta = min(tau_ta_raw.size, tau_ta_corr.size)

    return dict(
        stem=csv_path.stem,
        n_tracks=len(trajectories),
        dt=global_dt,
        # eaMSD
        tau=ea_raw.tau[:n_ea],
        ea_raw=ea_raw.msd[:n_ea],
        ea_raw_sem=ea_raw.msd_sem[:n_ea],
        ea_corr=ea_corr.msd[:n_ea],
        ea_corr_sem=ea_corr.msd_sem[:n_ea],
        # ⟨taMSD⟩  (may have different length)
        tau_ta=tau_ta_raw[:n_ta],
        ta_raw=msd_ta_raw[:n_ta],
        ta_raw_sem=sem_ta_raw[:n_ta],
        ta_corr=msd_ta_corr[:n_ta],
        ta_corr_sem=sem_ta_corr[:n_ta],
    )


# ── Per-group routine ──────────────────────────────────────────────────

def run_group(group_name: str, csv_files):
    print(f"\n{'='*60}\n[{group_name}]  {len(csv_files)} file(s)\n{'='*60}")
    results = []

    for csv_path in sorted(csv_files):
        print(f"  -> {csv_path.name}")
        r = analyse_file(csv_path)
        if r is None:
            print("     skipped (no/short trajectories)")
            continue
        results.append(r)

        short = r["stem"]

        for frac in FIT_FRACTIONS:
            pct = int(frac * 100)
            n_show_ea = max(4, int(frac * len(r["tau"])))
            n_show_ta = max(4, int(frac * len(r["tau_ta"])))

            # ── Per-file eaMSD overlay ─────────────────────────────
            _plot_overlay(
                r["tau"], r["ea_raw"], r["ea_raw_sem"],
                r["tau"], r["ea_corr"], r["ea_corr_sem"],
                ylabel=r"eaMSD [$\mu$m$^2$]",
                out_path=OUT_DIR / f"eaMSD_drift_cmp_{group_name}_{short}_f{pct:03d}.svg",
                n_show=n_show_ea,
            )

            # ── Per-file ⟨taMSD⟩ overlay ──────────────────────────
            _plot_overlay(
                r["tau_ta"], r["ta_raw"], r["ta_raw_sem"],
                r["tau_ta"], r["ta_corr"], r["ta_corr_sem"],
                ylabel=r"$\langle$taMSD$\rangle$ [$\mu$m$^2$]",
                out_path=OUT_DIR / f"taMSD_drift_cmp_{group_name}_{short}_f{pct:03d}.svg",
                n_show=n_show_ta,
            )

        print(f"     eaMSD: {r['ea_raw'].max():.3e} -> {r['ea_corr'].max():.3e} (max)"
              f"  taMSD: {r['ta_raw'].max():.3e} -> {r['ta_corr'].max():.3e} (max)")

    if not results:
        print(f"[{group_name}] no usable files.")
        return

    print(f"\n[{group_name}] Plots saved to:  {OUT_DIR}")


# ── Entry point ────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Compare MSD before/after drift correction."
    )
    parser.add_argument("--group", choices=["glic50", "glic200"],
                        help="Run only one group (default: both).")
    args = parser.parse_args()

    groups = {
        "glic50": sorted((DATA_DIR / GLIC50_SUBDIR).glob("*.csv")),
        "glic200": sorted((DATA_DIR / GLIC200_SUBDIR).glob("*.csv")),
    }

    targets = [args.group] if args.group else list(groups.keys())
    for g in targets:
        run_group(g, groups[g])

    print("\nDone.")


if __name__ == "__main__":
    main()
