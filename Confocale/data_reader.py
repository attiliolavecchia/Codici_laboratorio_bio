"""
Data reader utilities for trajectory CSV files.

This module focuses on reading exactly four required columns from a CSV file:
        - 'Track ID'    (identifier used to group points into trajectories)
        - 'X position'  (x coordinate)
        - 'Y position'  (y coordinate)
        - 'Time'        (time point)

The output is a simple in-memory representation:
        Mapping[Track ID, Trajectory]

Each Trajectory holds NumPy arrays for time, x, and y, sorted by time. We also
estimate a per-trajectory time step (dt) using the median of consecutive time
differences; a robust choice for mostly uniform sampling.

Why this structure?
        - Keeping arrays per track makes later vectorized MSD calculations simple.
        - Median dt is robust to rare glitches in acquisition timing.

Dependencies: pandas, numpy

References (APIs verified):
    - pandas.read_csv: https://pandas.pydata.org/docs/reference/api/pandas.read_csv.html
    - pandas.DataFrame.groupby: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.groupby.html
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Tuple, Union

import numpy as np
import pandas as pd


REQUIRED_COLUMNS = [
    "Track ID",
    "X position",
    "Y position",
    "Time",
]

# Accept common synonyms as seen in TrackMate and similar exports.
# We normalize headers by lowering, replacing non-alphanumerics with spaces, and collapsing spaces.
_HEADER_SYNONYMS: Dict[str, List[str]] = {
    "Track ID": ["track id", "track_id", "trackid", "id track"],
    "X position": ["x position", "x", "position x", "pos x", "position_x"],
    "Y position": ["y position", "y", "position y", "pos y", "position_y"],
    "Time": ["time", "t", "position t", "pos t", "position_t"],
}

def _normalize_header(name: str) -> str:
    import re
    s = str(name).strip().lower()
    s = re.sub(r"[^0-9a-z]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


@dataclass(frozen=True)
class Trajectory:
    """Simple container for one particle trajectory.

    Attributes:
        track_id: Identifier of the trajectory (string as read from CSV).
        time: 1D numpy array of time points (sorted ascending).
        x: 1D numpy array of x coordinates aligned with ``time``.
        y: 1D numpy array of y coordinates aligned with ``time``.
        dt: Estimated time step for this trajectory based on the median
            difference of consecutive time values. Can be NaN if fewer than two
            points are present.
    """

    track_id: Union[int, str]
    time: np.ndarray
    x: np.ndarray
    y: np.ndarray
    dt: float

    @property
    def n_points(self) -> int:
        """Return the number of time points for this trajectory."""
        return int(self.time.shape[0])


def normalize_and_select_required_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Return a DataFrame with exactly the required four columns in order.

    We accept case/whitespace variations in the incoming CSV headers and
    normalize them to the exact required names. If any required column is
    missing, we raise a clear error message.

    Args:
        df: A pandas DataFrame as initially read from the CSV.

    Returns:
        A new DataFrame with columns labeled exactly as in ``REQUIRED_COLUMNS``.

    Raises:
        ValueError: If any required column is not present in the input.
    """
    # Build mapping from normalized key -> original exact column name
    normalized_to_original: Dict[str, str] = {}
    for column_name in df.columns:
        key = _normalize_header(column_name)
        normalized_to_original[key] = column_name

    missing: List[str] = []
    selected_original: List[str] = []
    for required_exact in REQUIRED_COLUMNS:
        candidates = _HEADER_SYNONYMS.get(required_exact, [required_exact])
        found_col = None
        for cand in candidates:
            cand_key = _normalize_header(cand)
            if cand_key in normalized_to_original:
                found_col = normalized_to_original[cand_key]
                break
        if found_col is None:
            missing.append(required_exact)
        else:
            selected_original.append(found_col)

    if missing:
        raise ValueError(
            "CSV file is missing required columns: " + ", ".join(missing)
        )

    # Select and rename to the exact required names and order
    subset = df[selected_original].copy()
    subset.columns = REQUIRED_COLUMNS
    return subset


def coerce_required_column_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce dtypes for the four required columns.

    - 'Track ID' is treated as string (trimmed). This ensures consistent keys
      even if the CSV mixes int-like and string-like identifiers.
    - 'X position', 'Y position', and 'Time' are coerced to numeric floats with
      non-parsable values set to NaN (which we later drop).

    Args:
        df: DataFrame guaranteed to have exactly ``REQUIRED_COLUMNS``.

    Returns:
        A copy of the DataFrame with coerced dtypes.
    """
    df = df.copy()
    df["Track ID"] = df["Track ID"].astype("string").str.strip()
    for numeric_column in ("X position", "Y position", "Time"):
        df[numeric_column] = pd.to_numeric(df[numeric_column], errors="coerce")
    return df


def read_trajectories_from_csv(
    csv_path: str,
    encoding: Optional[str] = None,
    on_bad_lines: str = "error",
) -> Dict[Union[int, str], Trajectory]:
    """Read the CSV file and build a dictionary of trajectories.

    Main steps:
        1) Read only the required columns for safety and efficiency.
        2) Normalize column names and coerce expected dtypes.
        3) Drop rows with missing values in required fields.
        4) Group rows by 'Track ID', sort by 'Time', remove duplicate times, and
           convert to NumPy arrays per trajectory.
        5) Estimate a per-trajectory time step as the median difference of time.

    Args:
        csv_path: Path to the input CSV file.
        encoding: Optional encoding to pass to pandas.read_csv.
        on_bad_lines: Behavior on malformed lines; forwarded to pandas.read_csv.

    Returns:
        Dict mapping track ID (string) to a Trajectory dataclass.
    """
    # Read full CSV, then select required columns via robust normalization/synonyms.
    # This accommodates TrackMate-style headers like TRACK_ID, POSITION_X, T, etc.
    df = pd.read_csv(
        csv_path,
        encoding=encoding,
        on_bad_lines=on_bad_lines,
    )

    # 1) Normalize headers, 2) Coerce dtypes
    df = normalize_and_select_required_columns(df)
    df = coerce_required_column_dtypes(df)

    # 3) Drop rows with missing essential values
    df = df.dropna(subset=["Track ID", "X position", "Y position", "Time"]).copy()

    # 4) Group by Track ID and build Trajectory objects
    trajectories: Dict[Union[int, str], Trajectory] = {}
    for track_id, group in df.groupby("Track ID", sort=True):
        trajectory = _build_trajectory_from_group(track_id, group)
        trajectories[track_id] = trajectory

    return trajectories


def _build_trajectory_from_group(track_id: Union[int, str], group: pd.DataFrame) -> Trajectory:
    """Convert a per-track DataFrame into a Trajectory dataclass.

    We sort by time to ensure monotonic order, drop duplicate time stamps to
    avoid zero-lag duplicates, convert columns to NumPy arrays, and compute the
    median dt as a robust estimate of the time step.

    Args:
        track_id: Identifier of this trajectory.
        group: Sub-DataFrame containing only rows from one 'Track ID'.

    Returns:
        Trajectory dataclass with arrays and estimated dt.
    """
    group = group.sort_values("Time", kind="mergesort").copy()
    group = group.loc[~group["Time"].duplicated(keep="first")]

    time_values = np.asarray(group["Time"].to_numpy(), dtype=float)
    x_values = np.asarray(group["X position"].to_numpy(), dtype=float)
    y_values = np.asarray(group["Y position"].to_numpy(), dtype=float)

    if time_values.size < 2:
        dt_value = float("nan")
    else:
        dt_value = _median_time_step(time_values)

    return Trajectory(track_id=track_id, time=time_values, x=x_values, y=y_values, dt=dt_value)


def _median_time_step(time_values: np.ndarray) -> float:
    """Return the median time step from a sorted array of time values."""
    dtime = np.diff(time_values)
    return float(np.median(dtime))


def estimate_global_time_step(trajectories: Mapping[Union[int, str], Trajectory]) -> float:
    """Estimate a global Δt across all trajectories using the median of per-track dt.

    We ignore non-finite or non-positive dt values. If no valid dt is available,
    we fall back to 1.0 as a neutral default.

    Args:
        trajectories: Mapping from track id to Trajectory.

    Returns:
        A single float Δt to be used to convert lag steps (n) into τ = n Δt.
    """
    dt_values = [t.dt for t in trajectories.values() if np.isfinite(t.dt) and t.dt > 0]
    if not dt_values:
        return 1.0
    return float(np.median(np.asarray(dt_values, dtype=float)))


def read_trajectories(
    csv_path: str,
    encoding: Optional[str] = None,
    on_bad_lines: str = "error",
) -> Dict[Union[int, str], Trajectory]:
    """Backward-compatible alias for ``read_trajectories_from_csv``."""
    return read_trajectories_from_csv(csv_path, encoding=encoding, on_bad_lines=on_bad_lines)


def estimate_global_dt(trajectories: Mapping[Union[int, str], Trajectory]) -> float:
    """Backward-compatible alias for ``estimate_global_time_step``."""
    return estimate_global_time_step(trajectories)


__all__ = [
    "Trajectory",
    "read_trajectories_from_csv",
    "read_trajectories",
    "estimate_global_time_step",
    "estimate_global_dt",
    "normalize_and_select_required_columns",
    "coerce_required_column_dtypes",
    "REQUIRED_COLUMNS",
]
