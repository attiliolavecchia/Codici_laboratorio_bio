"""
EA‑MSD and TA‑MSD Plotter

Behavior:
    1) EA‑MSD comparison: scans Data/ for Experiment*_spots_<timestep>minstep.csv files,
         lets you pick FOUR minimal time steps, and overlays their EA‑MSD curves (ensemble, initial‑displacement only).
    2) TA‑MSD single‑trajectory study: lets you choose ONE trajectory from ONE file and
         overlays TA‑MSD computed from FOUR different truncation lengths (variable L of the SAME trajectory),
         as in the reference slide (single track with multiple lengths).

Headers are normalized by data_reader (TrackMate‑style POSITION_X/Y/T and TRACK_ID supported).

Output:
    - eamsd_plots/compare_eamsd_<timesteps>_f<frac>_<timestamp>.svg  (four minimal time steps)
    - tamsd_plots/compare_tamsd_ts<file_timestep>_track<ID>_L<lengths>_f<frac>_<timestamp>.svg  (single trajectory, four lengths)

Notes:
    - Lag extent can be limited via a max lag fraction prompt (0 < f ≤ 1).
    - Filenames are parameter‑encoded so runs with different settings do not overwrite each other.
    - EA‑MSD legend shows number of trajectories per dataset.
    - TA‑MSD legend shows truncation length L (points) for the single trajectory.
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Union

import matplotlib.pyplot as plt
import numpy as np

from data_reader import read_trajectories_from_csv, Trajectory
from msd_analyzer import (
    calculate_ensemble_msd,
    calculate_time_averaged_msd_per_track,
    determine_maximum_lag_steps,
)


# ---------- Utilities ----------

TIMESTEP_PATTERN = re.compile(r"_spots_(\d+)minstep\.csv$", re.IGNORECASE)


def find_data_dir(start: Path) -> Path:
    """Return the Data/ directory (searching typical locations relative to this script)."""
    # Preferred: Confocale/Data alongside this script
    candidate = (start / "Data").resolve()
    if candidate.is_dir():
        return candidate
    # Fallback: Codici/Data relative to repo root
    codici = start.parent
    fallback = (codici / "Data").resolve()
    if fallback.is_dir():
        return fallback
    # Last resort: current working directory/Data
    cwd_fallback = (Path.cwd() / "Data").resolve()
    return cwd_fallback


@dataclass(frozen=True)
class DataFile:
    path: Path
    timestep_min: int
    n_trajectories: int


def discover_experiment_files(data_dir: Path) -> List[DataFile]:
    files: List[DataFile] = []
    for p in sorted(data_dir.glob("*.csv")):
        m = TIMESTEP_PATTERN.search(p.name)
        if not m:
            continue
        try:
            minutes = int(m.group(1))
        except Exception:
            continue
        # Count trajectories for legend labeling and user verification
        try:
            tr = read_trajectories_from_csv(str(p))
            n_tr = len(tr)
        except Exception:
            n_tr = 0
        files.append(DataFile(path=p, timestep_min=minutes, n_trajectories=n_tr))
    # Sort by timestep
    files.sort(key=lambda df: df.timestep_min)
    return files


def prompt_for_four_timesteps(files: List[DataFile]) -> List[DataFile]:
    """Interactively select four distinct timestep files.

    Also supports non-interactive selection via environment variable
    COMPARE_MSD_CHOICES="idx1,idx2,idx3,idx4" to facilitate automation.
    """
    if len(files) < 4:
        raise RuntimeError(
            f"Need at least 4 files with pattern '*_spots_<N>minstep.csv' in Data/. Found {len(files)}."
        )

    print("Available minimal time steps (minutes):")
    for idx, df in enumerate(files, start=1):
        print(f"  [{idx}] {df.timestep_min:>4} min  |  {df.n_trajectories:>4} trajectories  —  {df.path.name}")
    import os
    preset = os.environ.get("COMPARE_MSD_CHOICES")
    if preset:
        parts = re.split(r"[\s,;]+", preset.strip())
        idxs = []
        for x in parts:
            if not x:
                continue
            if not x.isdigit():
                raise ValueError(f"Invalid index in COMPARE_MSD_CHOICES: {x}")
            idxs.append(int(x))
        if len(idxs) != 4 or len(set(idxs)) != 4 or any(i < 1 or i > len(files) for i in idxs):
            raise ValueError(
                "Environment COMPARE_MSD_CHOICES must contain four distinct valid indices (1-based)."
            )
        print(f"Using preset selection from COMPARE_MSD_CHOICES: {preset}")
        return [files[i - 1] for i in idxs]

    print("\nSelect FOUR distinct entries by index (comma-separated), e.g., 1,3,5,7")
    while True:
        raw = input("Your choice: ").strip()
        parts = re.split(r"[\s,;]+", raw)
        try:
            idxs = [int(x) for x in parts if x]
        except ValueError:
            print("Please enter only integers separated by commas/space.")
            continue
        if len(idxs) != 4:
            print("Please enter exactly FOUR indices.")
            continue
        if len(set(idxs)) != 4:
            print("Indices must be distinct.")
            continue
        # Range check
        if any(i < 1 or i > len(files) for i in idxs):
            print("One or more indices are out of range.")
            continue
        return [files[i - 1] for i in idxs]


# ---------- Computation ----------

@dataclass
class MSDCurves:
    label: str
    e_tau_s: List[float]
    e_msd: List[float]
    t_tau_s: List[float]
    t_msd: List[float]

def compute_curves_for_file(file: DataFile, *, max_lag_fraction: float | None) -> MSDCurves:
    trajectories = read_trajectories_from_csv(str(file.path))
    if not trajectories:
        raise RuntimeError(f"No trajectories found in {file.path}")

    # EA-MSD across all trajectories (initial-displacement only)
    ea = calculate_ensemble_msd(trajectories, max_lag_fraction=max_lag_fraction)

    # Label shows number of trajectories (user request)
    label = f"{file.n_trajectories} trajectories"
    return MSDCurves(
        label=label,
        e_tau_s=list(map(float, ea.tau)),
        e_msd=list(map(float, ea.msd)),
        t_tau_s=[],
        t_msd=[],
    )


def compute_tamsd_for_length(track: Trajectory, L: int, *, max_lag_fraction: float | None) -> MSDCurves:
    """Compute TA‑MSD for a truncated view of ``track`` with first ``L`` points.

    L must satisfy 2 ≤ L ≤ track.n_points. Legend label will be "L = <points>".
    """
    from dataclasses import replace
    if L < 2 or L > track.n_points:
        raise ValueError("L must be between 2 and the trajectory length")
    if L < track.n_points:
        track = replace(track, time=track.time[:L], x=track.x[:L], y=track.y[:L])
    ta = calculate_time_averaged_msd_per_track(track, max_lag_fraction=max_lag_fraction)
    return MSDCurves(
        label=f"L = {L} points",
        e_tau_s=[],
        e_msd=[],
        t_tau_s=list(map(float, ta.tau)),
        t_msd=list(map(float, ta.msd)),
    )


# ---------- Plotting ----------

COLORS = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"]

# Output subdirectories (created automatically inside the script directory)
EAMSD_DIRNAME = "eamsd_plots"
TAMSD_DIRNAME = "tamsd_plots"


def plot_overlaid_eamsd(curves: List[MSDCurves], out_path: Path) -> None:
    plt.figure(figsize=(8.0, 5.2))
    for i, c in enumerate(curves):
        color = COLORS[i % len(COLORS)]
        plt.plot(c.e_tau_s, c.e_msd, marker="o", ms=3.5, lw=1.5, color=color, label=c.label)
    plt.xlabel("Time Lag (τ) [seconds]")
    plt.ylabel("EA-MSD")
    #plt.title("Ensemble-Averaged MSD")
    plt.grid(True, linestyle=":", alpha=0.5)
    plt.legend(title="Trajectories")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def plot_overlaid_tamsd(curves: List[MSDCurves], out_path: Path) -> None:
    plt.figure(figsize=(8.0, 5.2))
    for i, c in enumerate(curves):
        color = COLORS[i % len(COLORS)]
        plt.plot(c.t_tau_s, c.t_msd, marker="s", ms=3.5, lw=1.5, color=color, label=c.label)
    plt.xlabel("Time Lag (τ) [seconds]")
    plt.ylabel("TA-MSD")
    #plt.title("Time-Averaged MSD")
    plt.grid(True, linestyle=":", alpha=0.5)
    plt.legend(title="Truncation length")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def _format_fraction(f: float | None) -> str:
    """Return a compact string for the max lag fraction.

    Examples: None -> "full", 0.5 -> "0.5", 0.25 -> "0.25", 1.0 -> "1".
    """
    if f is None:
        return "full"
    s = f"{f:.3f}"  # up to 3 decimals
    s = s.rstrip("0").rstrip(".")
    return s


def _timestamp() -> str:
    return datetime.now().strftime("%Y%m%d-%H%M%S")
def _prompt_max_lag_fraction() -> float | None:
    print("\nOptional: choose lag extent as a fraction of available points (0 < f <= 1).")
    print("Press Enter for full range (use N_max - 1).")
    while True:
        s = input("Max lag fraction [blank = full]: ").strip()
        if s == "":
            return None
        try:
            f = float(s)
        except ValueError:
            print("Please enter a number between 0 and 1, or leave blank.")
            continue
        if not (0 < f <= 1):
            print("Fraction must be in (0, 1].")
            continue
        return f


# ---------- Main ----------

def main(argv: List[str]) -> int:
    here = Path(__file__).resolve().parent
    data_dir = find_data_dir(here)

    if not data_dir.is_dir():
        print(f"Data directory not found: {data_dir}")
        return 2

    files = discover_experiment_files(data_dir)
    if not files:
        print("No matching CSV files found (pattern '*_spots_<N>minstep.csv').")
        return 3

    try:
        picked = prompt_for_four_timesteps(files)
    except KeyboardInterrupt:
        print("\nCancelled.")
        return 130
    except Exception as e:
        print(f"Error: {e}")
        return 4

    # Prompt for maximum lag fraction (applied to both EA and TA computations)
    try:
        max_lag_fraction = _prompt_max_lag_fraction()
    except KeyboardInterrupt:
        print("\nCancelled.")
        return 130

    # --- EA-MSD: four minimal time steps ---
    curves: List[MSDCurves] = []
    for df in picked:
        print(f"\nProcessing: {df.path.name}  (timestep ≈ {df.timestep_min} min, {df.n_trajectories} trajectories)")
        curves.append(compute_curves_for_file(df, max_lag_fraction=max_lag_fraction))

    # Enforce SAME time-lag (τ) grid across EA-MSD curves by interpolating to a common τ_grid
    # Choose τ_grid from dt_grid = max(dt_i) up to τ_max_common = min(max τ_i)
    ea_taus = [np.asarray(c.e_tau_s, dtype=float) for c in curves if c.e_tau_s]
    ea_msds = [np.asarray(c.e_msd, dtype=float) for c in curves if c.e_msd]
    if ea_taus:
        tau_max_common = float(min(t[-1] for t in ea_taus))
        # Estimate per-curve dt as first τ entry; guard against empty or zero
        dt_candidates = [float(t[0]) for t in ea_taus if t.size > 0 and t[0] > 0]
        if dt_candidates:
            dt_grid = float(max(dt_candidates))
            if tau_max_common > 0 and dt_grid > 0 and tau_max_common >= dt_grid:
                tau_grid = np.arange(dt_grid, tau_max_common + 1e-12, dt_grid)
                # Ensure at least two points; otherwise skip interpolation
                if tau_grid.size >= 1:
                    for c in curves:
                        t = np.asarray(c.e_tau_s, dtype=float)
                        m = np.asarray(c.e_msd, dtype=float)
                        if t.size >= 2:
                            new_m = np.interp(tau_grid, t, m)
                        elif t.size == 1:
                            new_m = np.full_like(tau_grid, m[0], dtype=float)
                        else:
                            new_m = np.asarray([], dtype=float)
                        c.e_tau_s = list(map(float, tau_grid))
                        c.e_msd = list(map(float, new_m))
                    print(f"Using common EA-MSD lag grid with Δt = {dt_grid:.6g} s and τ_max ≈ {tau_max_common:.6g} s.")

    # --- TA-MSD: one trajectory, four truncations ---
    print("\nTA-MSD (single trajectory with variable lengths)")
    print("Choose the source file for the trajectory (indices as above):")
    while True:
        sel_file = input(f"File index for TA-MSD (1..{len(files)}): ").strip()
        if not sel_file.isdigit():
            print("Enter an integer index.")
            continue
        fi = int(sel_file)
        if not (1 <= fi <= len(files)):
            print("Out of range.")
            continue
        ta_file = files[fi - 1]
        break

    traj_map = read_trajectories_from_csv(str(ta_file.path))
    if not traj_map:
        print("No trajectories found in selected file for TA-MSD.")
        return 5
    # Auto-select the longest trajectory by number of points
    longest_id = max(traj_map.keys(), key=lambda k: traj_map[k].n_points)
    track = traj_map[longest_id]
    N = track.n_points
    print(f"Selected longest trajectory automatically: ID {longest_id} (length: {N} points).")
    print("Enter FOUR distinct truncation lengths L (points), each with 2 ≤ L ≤ N.")
    print("Example: 20, 40, 60, 80")
    while True:
        raw = input("Lengths: ").strip()
        parts = re.split(r"[\s,;]+", raw)
        try:
            lengths = [int(x) for x in parts if x]
        except ValueError:
            print("Please enter integers only.")
            continue
        if len(lengths) != 4:
            print("Please enter exactly FOUR lengths.")
            continue
        if len(set(lengths)) != 4:
            print("Lengths must be distinct.")
            continue
        if any((L < 2 or L > N) for L in lengths):
            print(f"Each L must satisfy 2 ≤ L ≤ {N}.")
            continue
        lengths = sorted(lengths)
        break

    t_curves: List[MSDCurves] = []
    # Ensure SAME time-lag (τ) array for all truncation lengths by limiting to a common K.
    # Correct logic: fraction applies to the BASE (full) trajectory length N, then we also
    # cap by the shortest truncation length (min(L) - 1). This avoids the degenerate case
    # where small L values combined with a small fraction yield K=1.
    eff_fraction = 1.0 if (max_lag_fraction is None) else float(max_lag_fraction)
    K_fractional = determine_maximum_lag_steps(N, eff_fraction)
    K_len_limit = max(1, min(L - 1 for L in lengths))
    K_common = int(min(K_fractional, K_len_limit)) if lengths else 0
    if K_common < 1:
        print("Chosen lengths and fraction yield no valid lag steps (K_common < 1). Aborting TA-MSD plot.")
        return 6

    # Compute each TA-MSD with FULL available lags (ignore fraction here),
    # then trim to the common lag count so all τ are identical
    for L in lengths:
        c = compute_tamsd_for_length(track, L, max_lag_fraction=None)
        c.t_tau_s = c.t_tau_s[:K_common]
        c.t_msd = c.t_msd[:K_common]
        t_curves.append(c)

    # Inform the user about the enforced common τ range
    if t_curves:
        dt_seconds = t_curves[0].t_tau_s[0]  # this is 1 * dt
        # Guard in case of unexpected empty arrays
        if K_common >= 1 and dt_seconds:
            tau_max_seconds = K_common * dt_seconds
            print(
                "Using common lag steps K = "
                f"{K_common} (τ_max ≈ {tau_max_seconds:.6g} s) across all truncations "
                f"[limited by min(L)-1 = {K_len_limit}, fraction on N gave {K_fractional}]."
            )

    # ----- Build unique, parameter-encoded filenames to avoid overwrites -----
    f_str = _format_fraction(max_lag_fraction)
    ts_selected = "-".join(str(df.timestep_min) for df in picked)
    ts_now = _timestamp()
    # Ensure output directories exist
    e_dir = here / EAMSD_DIRNAME
    t_dir = here / TAMSD_DIRNAME
    e_dir.mkdir(parents=True, exist_ok=True)
    t_dir.mkdir(parents=True, exist_ok=True)

    out_e = e_dir / f"compare_eamsd_ts{ts_selected}_f{f_str}_{ts_now}.svg"

    ta_ts = ta_file.timestep_min
    L_str = "-".join(map(str, lengths))
    out_t = t_dir / f"compare_tamsd_ts{ta_ts}_track{longest_id}_L{L_str}_f{f_str}_{ts_now}.svg"
    plot_overlaid_eamsd(curves, out_e)
    plot_overlaid_tamsd(t_curves, out_t)
    print("\nSaved plots:")
    print(f"  EA-MSD: {out_e}")
    print(f"  TA-MSD: {out_t}")
    print(f"  (Directories: '{e_dir.name}', '{t_dir.name}')")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
