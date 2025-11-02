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
        tracks_per_lag: Number of trajectories contributing to each lag.
        dt: Global Δt used to convert steps to lag times.
        n_max: Maximum lag in steps used in the calculation.
    total_trajectories: Total number of trajectories (M) loaded.
    longest_trajectory_points: N_max, the number of points in the longest trajectory.
    """

    tau: np.ndarray
    msd: np.ndarray
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


def calculate_ensemble_msd(
    trajectories: Mapping[Union[int, str], Trajectory],
    *,
    max_lag_fraction: Optional[float] = None,
    global_dt: Optional[float] = None,
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

    Args:
        trajectories: Mapping of Track ID to Trajectory.
        max_lag_fraction: Optional fraction (0 < f ≤ 1] of the longest trajectory's points
            to cap the maximum lag steps (e.g., 0.25 → up to ⌊0.25·N_max⌋; capped at N_max−1).
            If None, use f = 1.0 (i.e., up to N_max−1).
        global_dt: Optional override for Δt in seconds. If None, a robust estimate is used
            (median of per-trajectory dt values > 0).

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

    per_track_values: List[np.ndarray] = []
    for track in trajectories.values():
        values = calculate_initial_displacement_msd_per_track(track, maximum_lag_steps)
        per_track_values.append(values)

    msd_values, tracks_per_lag = average_across_trajectories(per_track_values)

    return MSDResult(
        tau=tau,
        msd=msd_values,
        tracks_per_lag=tracks_per_lag,
        dt=float(global_dt),
        n_max=int(maximum_lag_steps),
        total_trajectories=int(len(trajectories)),
        longest_trajectory_points=int(maximum_points),
    )



def calculate_time_averaged_msd_per_track(
    track: Trajectory,
    *,
    max_lag_fraction: Optional[float] = None,
    dt_override: Optional[float] = None,
) -> MSDResult:
    """Compute the time-averaged MSD (TAMSD) for a single trajectory.

    Implements the definition:

        TAMSD(n) = (1 / (N - n)) * sum_{i=0..N-n-1} [ (x[i+n] - x[i])^2 + (y[i+n] - y[i])^2 ]

    for integer lags n = 1..K, where K is determined as a fraction of the
    available points. The result is returned in the same ``MSDResult`` container
    used elsewhere for convenient plotting alongside ensemble MSD.

    Args:
        track: The single trajectory to analyze.
        max_lag_fraction: Optional fraction (0 < f ≤ 1] of the track length to cap
            the maximum lag steps (at most N−1). If None, use the full range (f = 1.0).
        dt_override: If provided, this Δt (seconds) is used to convert steps to seconds.
            Otherwise we use ``track.dt`` if finite and positive; if that is
            not available, we fall back to 1.0 seconds.

    Returns:
        MSDResult with fields populated for this single-trajectory TAMSD:
            - tau: (K,) array, τ = 1..K multiplied by Δt (seconds)
            - msd: (K,) TAMSD per lag in squared position units
            - tracks_per_lag: array of ones of length K (one contributing track)
            - dt: Δt used (seconds)
            - n_max: K (maximum lag in steps)
            - total_trajectories: 1
            - longest_trajectory_points: N (= track.n_points)
        Returns an empty result if the track has fewer than 2 points.
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

    tamsd = np.full(K, np.nan, dtype=float)
    # For each lag n, average over all windows i where i+n < N
    for n in range(1, K + 1):
        dx = x[n:] - x[:-n]
        dy = y[n:] - y[:-n]
        vals = dx * dx + dy * dy
        # N - n windows; protected mean
        tamsd[n - 1] = float(np.mean(vals)) if vals.size else float("nan")

    tracks_per_lag = np.ones(K, dtype=int)

    return MSDResult(
        tau=tau,
        msd=tamsd,
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


def average_across_trajectories(per_track_means: List[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    """Average per-trajectory MSD values equally across trajectories.

    Args:
        per_track_means: List of arrays (one per trajectory), each shaped (K,)
            containing time-averaged MSD values per lag. Entries may be NaN for
            lags not supported by that trajectory's length.

    Returns:
        Tuple (msd_values, tracks_per_lag) where:
          - msd_values: (K,) ensemble-averaged MSD values per lag (nanmean).
          - tracks_per_lag: (K,) count of trajectories contributing a finite value per lag.
    """
    if not per_track_means:
        return np.asarray([], dtype=float), np.asarray([], dtype=int)

    stacked = np.vstack(per_track_means)
    msd_values = np.nanmean(stacked, axis=0)
    tracks_per_lag = np.sum(np.isfinite(stacked), axis=0).astype(int)
    return msd_values, tracks_per_lag


def MSDEmptyResult() -> MSDResult:
    """Return an empty MSDResult for edge cases (no data or too few points)."""
    return MSDResult(
        np.asarray([], dtype=float),
        np.asarray([], dtype=float),
        np.asarray([], dtype=int),
        float("nan"),
        0,
        0,
        0,
    )


def run_from_csv(csv_path: str, max_lag_fraction: Optional[float] = None) -> MSDResult:
    """Convenience function: read trajectories then compute ensemble MSD."""
    trajectories = read_trajectories_from_csv(csv_path)
    return calculate_ensemble_msd(trajectories, max_lag_fraction=max_lag_fraction)


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

    res = run_from_csv(args.csv, max_lag_fraction=args.max_lag_fraction)
    print(f"Δt (global) = {res.dt}")
    print(f"n_max (steps) = {res.n_max}, τ_max = {res.tau[-1] if res.tau.size else float('nan')}")
    print(f"Total trajectories (M) = {res.total_trajectories}")
    print(f"Longest trajectory length (points) = {res.longest_trajectory_points}")
    for t, m in zip(res.tau, res.msd):
        print(f"tau={t:.6g}, MSD={m:.6g}")

# Backward-compatible alias for earlier API name
def compute_ensemble_msd(
    trajectories: Mapping[Union[int, str], Trajectory],
    *,
    max_lag_fraction: float = 0.10,
    global_dt: Optional[float] = None,
) -> MSDResult:
    return calculate_ensemble_msd(trajectories, max_lag_fraction=max_lag_fraction, global_dt=global_dt)
