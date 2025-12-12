"""
FLIM Data Reader Module

This module provides functions to read FLIM (Fluorescence Lifetime Imaging Microscopy)
data from CSV files and convert time bin indices to actual time values.

In TCSPC (Time-Correlated Single-Photon Counting) FLIM, the data represents a histogram
of photon arrival times. Each "frame" or "bin" corresponds to a time channel in the
decay histogram.

For an 80 MHz pulsed laser:
- Laser period = 1 / 80 MHz = 12.5 ns
- With 128 time bins: Δt = 12.5 ns / 128 ≈ 0.0977 ns per bin

Author: Generated for FLIM analysis
Date: December 2025
"""

import numpy as np
import pandas as pd
from pathlib import Path
from dataclasses import dataclass
from typing import Tuple, Optional


# Default TCSPC parameters for 80 MHz laser with 128 time bins
DEFAULT_LASER_REP_RATE_MHZ = 80.0  # MHz
DEFAULT_NUM_TIME_BINS = 128
DEFAULT_TIME_RANGE_NS = 1000.0 / DEFAULT_LASER_REP_RATE_MHZ  # 12.5 ns


@dataclass(frozen=True)
class FLIMData:
    """
    Container for FLIM decay data.
    
    Attributes:
        time_ns: Array of time values in nanoseconds
        intensity: Array of fluorescence intensity values
        time_bin_width_ns: Width of each time bin in nanoseconds
        num_bins: Total number of time bins
        peak_index: Index of the maximum intensity (decay start)
        peak_time_ns: Time at peak intensity in nanoseconds
    """
    time_ns: np.ndarray
    intensity: np.ndarray
    time_bin_width_ns: float
    num_bins: int
    peak_index: int
    peak_time_ns: float


def read_flim_csv(
    csv_path: str,
    laser_rep_rate_mhz: float = DEFAULT_LASER_REP_RATE_MHZ,
    num_time_bins: Optional[int] = None,
    x_column: str = 'X',
    y_column: str = 'Y'
) -> FLIMData:
    """
    Read FLIM decay data from a CSV file and convert to time domain.
    
    The CSV file should contain two columns:
    - X: Time bin index (frame number)
    - Y: Mean fluorescence intensity
    
    Parameters:
        csv_path: Path to the CSV file containing FLIM data
        laser_rep_rate_mhz: Laser repetition rate in MHz (default: 80 MHz)
        num_time_bins: Number of time bins (default: auto-detected from data)
        x_column: Name of the column containing time bin indices
        y_column: Name of the column containing intensity values
    
    Returns:
        FLIMData object containing time and intensity arrays with metadata
    
    Raises:
        FileNotFoundError: If the CSV file does not exist
        ValueError: If required columns are not found in the CSV
    """
    csv_path = Path(csv_path)
    
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV file not found: {csv_path}")
    
    # Read the CSV file
    df = pd.read_csv(csv_path)
    
    # Validate columns
    if x_column not in df.columns:
        raise ValueError(f"Column '{x_column}' not found in CSV. Available columns: {list(df.columns)}")
    if y_column not in df.columns:
        raise ValueError(f"Column '{y_column}' not found in CSV. Available columns: {list(df.columns)}")
    
    # Extract data
    bin_indices = df[x_column].values.astype(float)
    intensity = df[y_column].values.astype(float)
    
    # Determine number of time bins
    if num_time_bins is None:
        num_time_bins = len(bin_indices)
    
    # Calculate time parameters
    # For 80 MHz laser: period = 12.5 ns
    time_range_ns = 1000.0 / laser_rep_rate_mhz  # Convert MHz to ns
    time_bin_width_ns = time_range_ns / num_time_bins
    
    # Convert bin indices to time in nanoseconds
    time_ns = bin_indices * time_bin_width_ns
    
    # Find peak (start of decay)
    peak_index = int(np.argmax(intensity))
    peak_time_ns = time_ns[peak_index]
    
    return FLIMData(
        time_ns=time_ns,
        intensity=intensity,
        time_bin_width_ns=time_bin_width_ns,
        num_bins=num_time_bins,
        peak_index=peak_index,
        peak_time_ns=peak_time_ns
    )


def extract_decay_region(
    data: FLIMData,
    start_index: Optional[int] = None,
    end_index: Optional[int] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract the decay region from FLIM data for fitting.
    
    By default, extracts data from the peak (maximum intensity) to the end,
    as this represents the fluorescence decay after excitation.
    
    Parameters:
        data: FLIMData object containing the full dataset
        start_index: Starting bin index (default: peak index)
        end_index: Ending bin index (default: last bin)
    
    Returns:
        Tuple of (time_ns, intensity) arrays for the decay region
    """
    if start_index is None:
        start_index = data.peak_index
    
    if end_index is None:
        end_index = len(data.intensity)
    
    # Ensure valid range
    start_index = max(0, start_index)
    end_index = min(len(data.intensity), end_index)
    
    time_decay = data.time_ns[start_index:end_index]
    intensity_decay = data.intensity[start_index:end_index]
    
    return time_decay, intensity_decay


def get_data_summary(data: FLIMData) -> str:
    """
    Generate a summary string of the FLIM data.
    
    Parameters:
        data: FLIMData object
    
    Returns:
        Formatted string with data summary
    """
    summary = f"""
FLIM Data Summary
=================
Number of time bins: {data.num_bins}
Time bin width: {data.time_bin_width_ns:.4f} ns
Total time range: {data.time_ns[-1]:.2f} ns
Peak intensity: {np.max(data.intensity):.4f}
Peak position: bin {data.peak_index} ({data.peak_time_ns:.2f} ns)
Intensity range: [{np.min(data.intensity):.4f}, {np.max(data.intensity):.4f}]
"""
    return summary.strip()


if __name__ == "__main__":
    # Example usage / test
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Read and display FLIM data from CSV file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python flim_data_reader.py "Plot Values.csv"
    python flim_data_reader.py "Plot Values.csv" --laser-rate 80
        """
    )
    
    parser.add_argument('csv', type=str, help='Path to the CSV file containing FLIM data')
    parser.add_argument('--laser-rate', type=float, default=DEFAULT_LASER_REP_RATE_MHZ,
                        help=f'Laser repetition rate in MHz (default: {DEFAULT_LASER_REP_RATE_MHZ})')
    
    args = parser.parse_args()
    
    # Read and display data
    data = read_flim_csv(args.csv, laser_rep_rate_mhz=args.laser_rate)
    print(get_data_summary(data))
    
    # Show decay region info
    time_decay, intensity_decay = extract_decay_region(data)
    print(f"\nDecay region: {len(time_decay)} points from {time_decay[0]:.2f} ns to {time_decay[-1]:.2f} ns")
