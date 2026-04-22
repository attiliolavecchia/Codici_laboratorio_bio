"""
MSD analyzer for 2D particle tracking data.

This module provides two complementary analyses:

1) Ensemble-averaged MSD using only the initial displacement per trajectory ("eaMSD, initial-displacement only").
   No time-averaging over starting indices is performed.

       MSD(τ = n Δt) = (1/M) ∑_{i=1..M} [ (x_i[n] − x_i[0])² + (y_i[n] − y_i[0])² ]

   For each lag n (in index steps), a trajectory contributes exactly one
   squared displacement relative to its initial point if it has at least n+1
   samples; values are then averaged across trajectories available at that lag.

2) Time-averaged MSD (TAMSD) for a single trajectory, which averages over all
   valid windows within that track:

       TAMSD(n) = (1/(N − n)) ∑_{i=0..N−n−1} [ (x[i+n] − x[i])² + (y[i+n] − y[i])² ]

Units and conventions:
    - Time lag τ is expressed in seconds using a global Δt (eaMSD) or a per-track/override Δt (TAMSD),
      both derived from the input Time column (seconds) unless explicitly overridden.
    - MSD values are in squared position units (e.g., micron² if the CSV positions are in micron).
    - Input trajectories are provided by ``data_reader`` which normalizes common TrackMate headers and
      estimates per-trajectory Δt via the median of time differences.

Dependencies: numpy (computations). Data loading utilities are in data_reader.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Mapping, Optional, Tuple, Union

import numpy as np

from data_reader import (
    Trajectory,
    read_trajectories_from_csv,
    estimate_global_time_step,
)


@dataclass(frozen=True)
class MSDResult:
    """Container for the ensemble-averaged MSD vs. lag.

    Attributes:
        tau: 1D array of lag times τ = n Δt in the same units as the input Time.
        msd: 1D array of MSD values (ensemble averaged) corresponding to ``tau``.
        msd_std: 1D array of sample standard deviation of r² values per lag (ddof=1).
        msd_sem: 1D array of standard error of the mean (SEM = msd_std / sqrt(M(n))) per lag.
        tracks_per_lag: Number of trajectories contributing to each lag.
        dt: Global Δt used to convert steps to lag times.
        n_max: Maximum lag in steps used in the calculation.
    total_trajectories: Total number of trajectories (M) loaded.
    longest_trajectory_points: N_max, the number of points in the longest trajectory.
    """

    tau: np.ndarray
    msd: np.ndarray
    msd_std: np.ndarray
    msd_sem: np.ndarray
    tracks_per_lag: np.ndarray
    dt: float
    n_max: int
    total_trajectories: int
    longest_trajectory_points: int

    def as_pairs(self) -> List[Tuple[float, float]]:
        """Return the result as a list of (tau, msd) tuples for convenience."""
        return [(float(t), float(m)) for t, m in zip(self.tau, self.msd)]


def calculate_initial_displacement_msd_per_track(track: Trajectory, maximum_lag_steps: int) -> np.ndarray:
    """Compute per-lag MSD using ONLY initial displacement for a single track.

    For each lag n = 1..maximum_lag_steps, we compute exactly one squared
    displacement relative to the initial point (index 0):

        (x[n] - x[0])^2 + (y[n] - y[0])^2

    Args:
        track: A single trajectory containing x, y, and time arrays.
        maximum_lag_steps: The largest step lag to consider (in index steps).

    Returns:
        A (maximum_lag_steps,) float array where entry n-1 holds the MSD value
        for lag n if n < track length, otherwise NaN.
    """
    x_positions = track.x
    y_positions = track.y
    num_points = track.n_points

    per_lag_values = np.full(maximum_lag_steps, np.nan, dtype=float)

    if num_points < 2:
        return per_lag_values

    x0 = float(x_positions[0])
    y0 = float(y_positions[0])

    # For lag n, compare position at index n with index 0
    max_valid_lag = min(maximum_lag_steps, max(0, num_points - 1))
    for lag_steps in range(1, max_valid_lag + 1):
        dx = float(x_positions[lag_steps]) - x0
        dy = float(y_positions[lag_steps]) - y0
        per_lag_values[lag_steps - 1] = dx * dx + dy * dy

    return per_lag_values


def _collect_initial_displacements_per_track(
    track: Trajectory, maximum_lag_steps: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Collect per-lag Δx and Δy (signed, not squared) for drift-corrected MSD.

    For each lag n = 1..maximum_lag_steps:
        dx[n-1] = x[n] - x[0],   dy[n-1] = y[n] - y[0]

    Returns:
        (dx_array, dy_array) each of shape (maximum_lag_steps,), with NaN
        for lags beyond the track length.
    """
    num_points = track.n_points
    dx_arr = np.full(maximum_lag_steps, np.nan, dtype=float)
    dy_arr = np.full(maximum_lag_steps, np.nan, dtype=float)

    if num_points < 2:
        return dx_arr, dy_arr

    x = np.asarray(track.x, dtype=float)
    y = np.asarray(track.y, dtype=float)
    x0, y0 = x[0], y[0]

    max_valid = min(maximum_lag_steps, num_points - 1)
    dx_arr[:max_valid] = x[1:max_valid + 1] - x0
    dy_arr[:max_valid] = y[1:max_valid + 1] - y0
    return dx_arr, dy_arr


def calculate_ensemble_msd(
    trajectories: Mapping[Union[int, str], Trajectory],
    *,
    max_lag_fraction: Optional[float] = None,
    global_dt: Optional[float] = None,
    drift_corrected: bool = True,
) -> MSDResult:
    """Compute the 2D ensemble-averaged MSD as a function of lag time τ.

     High-level steps:
        1) Determine the maximum lag in index steps as a user-specified fraction of the
            longest trajectory length. If ``max_lag_fraction`` is None, use the full
            supported range up to N_max−1.
        2) For each trajectory, compute per-lag MSD using ONLY initial
            displacement (index 0) as reference — no time-averaging.
        3) Average those per-trajectory values equally across trajectories for each lag (NaNs ignored).
        4) Convert lag steps (n) to τ = n Δt using a global Δt (seconds).

    When ``drift_corrected=True``, the variance method is used instead of the
    raw ensemble mean (Michalet 2010, Phys. Rev. E 82, 041914):

        MSD_corr(n) = Var(Δx_n) + Var(Δy_n)

    This removes the coherent drift component |⟨Δr⟩|² from the MSD.

    Args:
        trajectories: Mapping of Track ID to Trajectory.
        max_lag_fraction: Optional fraction (0 < f ≤ 1] of the longest trajectory's points
            to cap the maximum lag steps (e.g., 0.25 → up to ⌊0.25·N_max⌋; capped at N_max−1).
            If None, use f = 1.0 (i.e., up to N_max−1).
        global_dt: Optional override for Δt in seconds. If None, a robust estimate is used
            (median of per-trajectory dt values > 0).
        drift_corrected: If True, subtract the ensemble-mean displacement squared
            from the MSD at each lag (variance method).

    Returns:
        MSDResult with:
            - tau: (K,) lag times τ = n·Δt in seconds
            - msd: (K,) ensemble-averaged MSD per lag in squared position units
            - tracks_per_lag: (K,) number of contributing trajectories per lag
            - dt: global Δt used (seconds)
            - n_max: maximum lag in steps K
            - total_trajectories: total number of trajectories loaded (M)
            - longest_trajectory_points: length in points of the longest trajectory (N_max)
    """
    if not trajectories:
        return MSDEmptyResult()

    maximum_points = max(t.n_points for t in trajectories.values())
    if maximum_points < 2:
        return MSDEmptyResult()

    # Resolve lag fraction (None ⇒ use full range)
    eff_fraction = 1.0 if (max_lag_fraction is None) else float(max_lag_fraction)
    maximum_lag_steps = determine_maximum_lag_steps(maximum_points, eff_fraction)

    if global_dt is None:
        global_dt = estimate_global_time_step(trajectories)

    tau = build_tau_array(maximum_lag_steps, global_dt)

    if drift_corrected:
        # Variance method: MSD_corr(n) = Var(Δx_n) + Var(Δy_n)
        dx_all: List[np.ndarray] = []
        dy_all: List[np.ndarray] = []
        for track in trajectories.values():
            dx, dy = _collect_initial_displacements_per_track(track, maximum_lag_steps)
            dx_all.append(dx)
            dy_all.append(dy)

        dx_stack = np.vstack(dx_all)  # (M, K)
        dy_stack = np.vstack(dy_all)  # (M, K)

        # Per-lag variance (ignoring NaN for short tracks)
        with np.errstate(invalid='ignore'):
            msd_values = np.nanvar(dx_stack, axis=0, ddof=1) + np.nanvar(dy_stack, axis=0, ddof=1)
            tracks_per_lag = np.sum(np.isfinite(dx_stack), axis=0).astype(int)

            # SEM for the variance estimator: SE(Var) ≈ Var * sqrt(2/(M-1))
            # Combined SE for Var_x + Var_y propagated in quadrature
            var_x = np.nanvar(dx_stack, axis=0, ddof=1)
            var_y = np.nanvar(dy_stack, axis=0, ddof=1)
            M = tracks_per_lag.astype(float)
            se_var_x = var_x * np.sqrt(2.0 / np.maximum(M - 1, 1))
            se_var_y = var_y * np.sqrt(2.0 / np.maximum(M - 1, 1))
            msd_sem = np.sqrt(se_var_x**2 + se_var_y**2)
            msd_std = msd_values * np.sqrt(2.0 / np.maximum(M - 1, 1))
    else:
        per_track_values: List[np.ndarray] = []
        for track in trajectories.values():
            values = calculate_initial_displacement_msd_per_track(track, maximum_lag_steps)
            per_track_values.append(values)

        msd_values, msd_std, tracks_per_lag = average_across_trajectories(per_track_values)
        with np.errstate(invalid='ignore'):
            msd_sem = msd_std / np.sqrt(tracks_per_lag.astype(float))

    return MSDResult(
        tau=tau,
        msd=msd_values,
        msd_std=msd_std,
        msd_sem=msd_sem,
        tracks_per_lag=tracks_per_lag,
        dt=float(global_dt),
        n_max=int(maximum_lag_steps),
        total_trajectories=int(len(trajectories)),
        longest_trajectory_points=int(maximum_points),
    )



@dataclass(frozen=True)
class EnsembleDrift:
    """Common-mode drift trajectory R(t_k) estimated from an ensemble.

    Attributes:
        time: 1D array of time points (seconds), aligned to a global frame grid.
        Rx, Ry: cumulative ensemble drift positions at each time.
        M_per_step: number of tracks contributing to each single-step displacement.
        min_M: minimum M observed across all retained steps.
    """
    time: np.ndarray
    Rx: np.ndarray
    Ry: np.ndarray
    M_per_step: np.ndarray
    min_M: int


def compute_ensemble_drift(
    trajectories: Mapping[Union[int, str], Trajectory],
    *,
    global_dt: Optional[float] = None,
    min_tracks_per_step: int = 5,
) -> EnsembleDrift:
    """Estimate the common-mode drift trajectory from an ensemble of tracks.

    Algorithm (Crocker & Grier 1996, J. Colloid Interface Sci. 179:298, §III-D):
        1) Build a global frame grid t_k = k·Δt covering the union of all tracks.
        2) For each consecutive frame pair (k, k+1), compute the average single-step
           displacement ⟨Δr⟩(t_k) over all tracks present at BOTH frames k and k+1.
        3) The ensemble-drift position is the cumulative sum:
                R(t_k) = Σ_{j<k} ⟨Δr⟩(t_j).
        4) Steps with fewer than ``min_tracks_per_step`` contributing tracks are
           reported via ``min_M`` so callers can warn the user.

    Why this estimator (and not per-track OLS):
        For pure Brownian motion with no drift, the OLS slope of x(t) on a single
        track has variance Var(v̂) ≈ 2D/T (Qian, Sheetz, Elson 1991, Biophys. J.
        60:910), so subtracting it removes part of the genuine diffusion signal
        and biases the TAMSD as ⟨TAMSD_corr⟩ ≈ 4Dτ·(1 − τ/T·h(τ/T)). Averaging
        single-step displacements across M tracks suppresses the diffusive
        contribution by 1/√M while preserving the common drift, so the TAMSD
        bias drops to O(τ/(M·T)) (Vestergaard, Blainey, Flyvbjerg 2014, PRE
        89:022726; Manzo & Garcia-Parajo 2015, Rep. Prog. Phys. 78:124601).

    Args:
        trajectories: Mapping of Track ID to Trajectory. All tracks are assumed
            to share a common Δt (validated against ``global_dt``).
        global_dt: Optional override for Δt (seconds). If None, a robust estimate
            is computed via ``estimate_global_time_step``.
        min_tracks_per_step: Minimum number of tracks that must contribute to a
            given single-step displacement for it to be considered reliable.
            Steps below this threshold are kept (the cumulative sum still includes
            them) but reported via ``min_M`` so callers can decide how to act.

    Returns:
        EnsembleDrift with cumulative R(t_k) for k = 0..K_max−1 (R(t_0) = 0).
    """
    if not trajectories:
        return EnsembleDrift(
            np.asarray([], float), np.asarray([], float), np.asarray([], float),
            np.asarray([], int), 0,
        )

    if global_dt is None:
        global_dt = estimate_global_time_step(trajectories)

    # Build global frame index for each sample of each track
    # frame_k = round((t - t_min_global) / Δt)
    t_min = min(float(np.min(tr.time)) for tr in trajectories.values())
    t_max = max(float(np.max(tr.time)) for tr in trajectories.values())
    K_total = int(round((t_max - t_min) / global_dt)) + 1   # number of frames

    # Per-step accumulators: sum of Δx, Δy, and count of contributing tracks
    sum_dx = np.zeros(K_total - 1, dtype=float)
    sum_dy = np.zeros(K_total - 1, dtype=float)
    cnt    = np.zeros(K_total - 1, dtype=int)

    for tr in trajectories.values():
        t = np.asarray(tr.time, dtype=float)
        x = np.asarray(tr.x, dtype=float)
        y = np.asarray(tr.y, dtype=float)
        if t.size < 2:
            continue
        # Frame indices (0-based, on the global grid)
        frames = np.rint((t - t_min) / global_dt).astype(int)
        # Single-step displacements only between *consecutive* frames
        df = np.diff(frames)            # frame gaps
        consecutive = (df == 1)
        if not consecutive.any():
            continue
        dx = np.diff(x)[consecutive]
        dy = np.diff(y)[consecutive]
        idx = frames[:-1][consecutive]  # the starting frame of each step
        # Accumulate (vectorized)
        np.add.at(sum_dx, idx, dx)
        np.add.at(sum_dy, idx, dy)
        np.add.at(cnt,    idx, 1)

    # Average single-step displacement at each step (NaN where no track contributed)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean_dx = np.where(cnt > 0, sum_dx / np.maximum(cnt, 1), 0.0)
        mean_dy = np.where(cnt > 0, sum_dy / np.maximum(cnt, 1), 0.0)
    # Cumulative drift trajectory; R(t_0) = 0 by construction
    Rx = np.concatenate([[0.0], np.cumsum(mean_dx)])
    Ry = np.concatenate([[0.0], np.cumsum(mean_dy)])
    time_grid = t_min + np.arange(K_total) * global_dt

    # Report the worst-step coverage among steps that had at least one track
    populated = cnt[cnt > 0]
    min_M = int(populated.min()) if populated.size else 0
    if populated.size and min_M < min_tracks_per_step:
        n_thin = int(np.sum(cnt < min_tracks_per_step))
        # warn but do not fail — the caller can decide
        import warnings as _w
        _w.warn(
            f"compute_ensemble_drift: {n_thin}/{cnt.size} steps have "
            f"M < {min_tracks_per_step} tracks (min M = {min_M}). "
            "Drift estimate may be noisy at those lags.",
            RuntimeWarning, stacklevel=2,
        )

    return EnsembleDrift(
        time=time_grid, Rx=Rx, Ry=Ry, M_per_step=cnt, min_M=min_M,
    )


def calculate_time_averaged_msd_per_track(
    track: Trajectory,
    *,
    max_lag_fraction: Optional[float] = None,
    dt_override: Optional[float] = None,
    drift_corrected: bool = True,
    ensemble_drift: Optional[EnsembleDrift] = None,
) -> MSDResult:
    """Compute the time-averaged MSD (TAMSD) for a single trajectory.

    Implements the definition:

        TAMSD(n) = (1 / (N - n)) * sum_{i=0..N-n-1} [ (x[i+n] - x[i])^2 + (y[i+n] - y[i])^2 ]

    for integer lags n = 1..K, where K is determined as a fraction of the
    available points. The result is returned in the same ``MSDResult`` container
    used elsewhere for convenient plotting alongside ensemble MSD.

    Drift correction (when ``drift_corrected=True``):
        The common-mode drift R(t_k) computed from the full ensemble (see
        ``compute_ensemble_drift``) is subtracted from this track's coordinates:
        x_corr(t_k) = x(t_k) − R_x(t_k), y analogously. This preserves the
        track's individual diffusive fluctuations because R is built from a
        Σ Δr / M average over all tracks (Crocker & Grier 1996, J. Colloid
        Interface Sci. 179:298; Vestergaard et al. 2014, PRE 89:022726).

        Per-track OLS detrending is intentionally NOT supported here: it biases
        the TAMSD downward as ⟨TAMSD_corr⟩ ≈ 4Dτ·(1 − τ/T·h(τ/T)) because the
        OLS slope on a single Brownian track has variance 2D/T and absorbs part
        of the genuine diffusion (Qian, Sheetz, Elson 1991, Biophys. J. 60:910).

    Args:
        track: The single trajectory to analyze.
        max_lag_fraction: Optional fraction (0 < f ≤ 1] of the track length to cap
            the maximum lag steps (at most N−1). If None, use the full range (f = 1.0).
        dt_override: If provided, this Δt (seconds) is used to convert steps to seconds.
            Otherwise we use ``track.dt`` if finite and positive; if that is
            not available, we fall back to 1.0 seconds.
        drift_corrected: If True, subtract the ensemble drift before computing TAMSD.
            Requires ``ensemble_drift`` to be provided.
        ensemble_drift: Pre-computed common-mode drift for the file the track
            belongs to (see ``compute_ensemble_drift``). Required when
            ``drift_corrected=True``.

    Returns:
        MSDResult with fields populated for this single-trajectory TAMSD.
        Returns an empty result if the track has fewer than 2 points.

    Raises:
        ValueError: If ``drift_corrected=True`` but ``ensemble_drift`` is not provided.
    """
    N = int(track.n_points)
    if N < 2:
        return MSDEmptyResult()

    eff_fraction = 1.0 if (max_lag_fraction is None) else float(max_lag_fraction)
    K = determine_maximum_lag_steps(N, eff_fraction)

    # Determine Δt to express τ in seconds
    if dt_override is not None and np.isfinite(dt_override) and dt_override > 0:
        dt = float(dt_override)
    else:
        dt = float(track.dt) if (np.isfinite(track.dt) and track.dt > 0) else 1.0

    tau = build_tau_array(K, dt)

    x = np.asarray(track.x, dtype=float)
    y = np.asarray(track.y, dtype=float)

    if drift_corrected:
        if ensemble_drift is None or ensemble_drift.time.size == 0:
            raise ValueError(
                "calculate_time_averaged_msd_per_track: drift_corrected=True "
                "requires a precomputed ensemble_drift (see compute_ensemble_drift). "
                "Per-track OLS detrending is no longer supported because it biases "
                "the TAMSD by a factor (1 − τ/T·h(τ/T)) [Qian-Sheetz-Elson 1991, "
                "Vestergaard et al. 2014]."
            )
        # Ensemble (common-mode) drift subtraction — Crocker & Grier 1996.
        # Map this track's time samples to the global frame grid by index.
        t = np.asarray(track.time, dtype=float)
        t_min = float(ensemble_drift.time[0])
        frames = np.rint((t - t_min) / dt).astype(int)
        frames = np.clip(frames, 0, ensemble_drift.time.size - 1)
        x = x - ensemble_drift.Rx[frames]
        y = y - ensemble_drift.Ry[frames]

    tamsd = np.full(K, np.nan, dtype=float)
    tamsd_std = np.full(K, np.nan, dtype=float)
    tamsd_sem = np.full(K, np.nan, dtype=float)
    # For each lag n, average over all windows i where i+n < N
    for n in range(1, K + 1):
        dx = x[n:] - x[:-n]
        dy = y[n:] - y[:-n]
        vals = dx * dx + dy * dy
        # N - n windows; protected mean
        tamsd[n - 1] = float(np.mean(vals)) if vals.size else float("nan")
        if vals.size >= 2:
            std = float(np.std(vals, ddof=1))
            tamsd_std[n - 1] = std
            # Due to overlapping windows, the effective number of independent samples
            # is roughly the total windows divided by the lag step n.
            n_eff = max(1.0, float(vals.size) / n)
            tamsd_sem[n - 1] = std / np.sqrt(n_eff)

    tracks_per_lag = np.ones(K, dtype=int)

    return MSDResult(
        tau=tau,
        msd=tamsd,
        msd_std=tamsd_std,
        msd_sem=tamsd_sem,
        tracks_per_lag=tracks_per_lag,
        dt=float(dt),
        n_max=int(K),
        total_trajectories=1,
        longest_trajectory_points=N,
    )



def determine_maximum_lag_steps(maximum_points: int, max_lag_fraction: float) -> int:
    """Compute the maximum lag in steps from the longest trajectory.

    We choose: max(1, floor(Nmax * fraction)), but never more than Nmax - 1 to
    ensure each lag has at least one valid displacement window.
    """
    # Clamp fraction into (0, 1] to avoid surprises
    f = float(max_lag_fraction)
    if not np.isfinite(f) or f <= 0:
        f = 1.0
    if f > 1.0:
        f = 1.0
    candidate = max(1, int(np.floor(maximum_points * f)))
    return min(candidate, max(1, maximum_points - 1))


def build_tau_array(maximum_lag_steps: int, global_dt: float) -> np.ndarray:
    """Create the τ array from lag steps and global Δt: τ = n Δt for n = 1..K."""
    return (np.arange(1, maximum_lag_steps + 1, dtype=float)) * float(global_dt)


def average_across_trajectories(per_track_means: List[np.ndarray]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Average per-trajectory MSD values equally across trajectories.

    Args:
        per_track_means: List of arrays (one per trajectory), each shaped (K,)
            containing time-averaged MSD values per lag. Entries may be NaN for
            lags not supported by that trajectory's length.

    Returns:
        Tuple (msd_values, msd_std, tracks_per_lag) where:
          - msd_values: (K,) ensemble-averaged MSD values per lag (nanmean).
          - msd_std: (K,) sample standard deviation of per-trajectory r² values (ddof=1).
          - tracks_per_lag: (K,) count of trajectories contributing a finite value per lag.
    """
    if not per_track_means:
        empty = np.asarray([], dtype=float)
        return empty, empty.copy(), np.asarray([], dtype=int)

    stacked = np.vstack(per_track_means)
    msd_values = np.nanmean(stacked, axis=0)
    msd_std = np.nanstd(stacked, axis=0, ddof=1)
    tracks_per_lag = np.sum(np.isfinite(stacked), axis=0).astype(int)
    return msd_values, msd_std, tracks_per_lag


def MSDEmptyResult() -> MSDResult:
    """Return an empty MSDResult for edge cases (no data or too few points)."""
    return MSDResult(
        np.asarray([], dtype=float),
        np.asarray([], dtype=float),
        np.asarray([], dtype=float),
        np.asarray([], dtype=float),
        np.asarray([], dtype=int),
        float("nan"),
        0,
        0,
        0,
    )


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Compute 2D ensemble-averaged MSD from a CSV file.")
    p.add_argument("csv", help="Path to sferette240nm_spots.csv (or similar)")
    p.add_argument(
        "--max-lag-fraction",
        type=float,
        default=None,
        help="Fraction (0<f<=1] of longest track to cap maximum lag; omit for full range (N_max-1)",
    )
    args = p.parse_args()

    trajectories = read_trajectories_from_csv(args.csv)
    res = calculate_ensemble_msd(trajectories, max_lag_fraction=args.max_lag_fraction)
    print(f"Δt (global) = {res.dt}")
    print(f"n_max (steps) = {res.n_max}, τ_max = {res.tau[-1] if res.tau.size else float('nan')}")
    print(f"Total trajectories (M) = {res.total_trajectories}")
    print(f"Longest trajectory length (points) = {res.longest_trajectory_points}")
    for t, m in zip(res.tau, res.msd):
        print(f"tau={t:.6g}, MSD={m:.6g}")
