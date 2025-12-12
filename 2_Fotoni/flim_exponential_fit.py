"""
FLIM Exponential Decay Fitting Script

This script fits fluorescence decay data from FLIM (Fluorescence Lifetime Imaging 
Microscopy) experiments using exponential decay models:

Mono-exponential: I(t) = I₀ × exp(-t/τ)
Bi-exponential:   I(t) = A₁ × exp(-t/τ₁) + A₂ × exp(-t/τ₂)

The fitting uses the Levenberg-Marquardt algorithm via scipy.optimize.curve_fit.

For TCSPC data with an 80 MHz laser:
- Laser period = 12.5 ns
- Time bin width = 12.5 ns / 128 bins ≈ 0.0977 ns

Author: Generated for FLIM analysis
Date: December 2025
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Tuple
import argparse

# Import the data reader module
from flim_data_reader import (
    read_flim_csv, 
    extract_decay_region, 
    get_data_summary,
    FLIMData,
    DEFAULT_LASER_REP_RATE_MHZ
)


# Output directory for plots
OUTPUT_DIR = Path(__file__).parent / "img" / "Rodamina_FLIM"


# =============================================================================
# Exponential Decay Models
# =============================================================================

def mono_exponential(t: np.ndarray, I0: float, tau: float) -> np.ndarray:
    """
    Mono-exponential decay model.
    
    I(t) = I₀ × exp(-t/τ)
    
    Parameters:
        t: Time array (ns)
        I0: Initial intensity at t=0
        tau: Fluorescence lifetime (ns)
    
    Returns:
        Intensity array
    """
    return I0 * np.exp(-t / tau)


def bi_exponential(t: np.ndarray, A1: float, tau1: float, A2: float, tau2: float) -> np.ndarray:
    """
    Bi-exponential decay model.
    
    I(t) = A₁ × exp(-t/τ₁) + A₂ × exp(-t/τ₂)
    
    This model is used when multiple decay components are present,
    e.g., different fluorophore environments, FRET, or multiple species.
    
    Parameters:
        t: Time array (ns)
        A1: Amplitude of first decay component
        tau1: Lifetime of first component (ns)
        A2: Amplitude of second decay component
        tau2: Lifetime of second component (ns)
    
    Returns:
        Intensity array
    """
    return A1 * np.exp(-t / tau1) + A2 * np.exp(-t / tau2)


# =============================================================================
# Fit Result Dataclasses
# =============================================================================

@dataclass(frozen=True)
class MonoExpFitResult:
    """
    Results from mono-exponential decay fitting.
    
    Attributes:
        I0: Initial intensity
        I0_error: Standard error of I0
        tau: Fluorescence lifetime (ns)
        tau_error: Standard error of tau (ns)
        R_squared: Coefficient of determination (R²)
        fit_start_ns: Start time of fitted region (ns)
        fit_end_ns: End time of fitted region (ns)
    """
    I0: float
    I0_error: float
    tau: float
    tau_error: float
    R_squared: float
    fit_start_ns: float
    fit_end_ns: float


@dataclass(frozen=True)
class BiExpFitResult:
    """
    Results from bi-exponential decay fitting.
    
    Attributes:
        A1: Amplitude of first component
        A1_error: Standard error of A1
        tau1: Lifetime of first component (ns)
        tau1_error: Standard error of tau1 (ns)
        A2: Amplitude of second component
        A2_error: Standard error of A2
        tau2: Lifetime of second component (ns)
        tau2_error: Standard error of tau2 (ns)
        tau_avg: Amplitude-weighted average lifetime (ns)
        R_squared: Coefficient of determination (R²)
        fit_start_ns: Start time of fitted region (ns)
        fit_end_ns: End time of fitted region (ns)
    """
    A1: float
    A1_error: float
    tau1: float
    tau1_error: float
    A2: float
    A2_error: float
    tau2: float
    tau2_error: float
    tau_avg: float
    R_squared: float
    fit_start_ns: float
    fit_end_ns: float


# =============================================================================
# Fitting Functions
# =============================================================================

def calculate_r_squared(y_observed: np.ndarray, y_predicted: np.ndarray) -> float:
    """
    Calculate the coefficient of determination (R²).
    
    R² = 1 - SS_res / SS_tot
    
    Parameters:
        y_observed: Observed data values
        y_predicted: Predicted values from the model
    
    Returns:
        R² value (0 to 1, where 1 is perfect fit)
    """
    ss_res = np.sum((y_observed - y_predicted) ** 2)
    ss_tot = np.sum((y_observed - np.mean(y_observed)) ** 2)
    
    if ss_tot == 0:
        return 0.0
    
    return 1.0 - (ss_res / ss_tot)


def auto_fit_indices_from_peak_to_last_nonzero(
    data: FLIMData,
    epsilon: float = 0.0
) -> Tuple[int, int]:
    """
    Pick fit window starting at the peak and ending at the last non-zero point.

    Parameters:
        data: FLIMData containing the full trace
        epsilon: Values greater than this are treated as non-zero

    Returns:
        (start_idx, end_idx) suitable for slicing (end is exclusive)
    """
    start_idx = data.peak_index

    # Look for the last non-zero sample at or after the peak.
    after_peak_nonzero = np.where(data.intensity[start_idx:] > epsilon)[0]

    if len(after_peak_nonzero) == 0:
        # Fallback: use full array to avoid empty window
        end_idx = len(data.intensity)
    else:
        last_nonzero_idx = start_idx + after_peak_nonzero[-1]
        end_idx = min(len(data.intensity), last_nonzero_idx + 1)

    # Ensure at least two points for fitting
    if end_idx - start_idx < 2:
        end_idx = min(len(data.intensity), start_idx + 2)

    return start_idx, end_idx


def fit_mono_exponential(
    time_ns: np.ndarray,
    intensity: np.ndarray,
    p0: Optional[Tuple[float, float]] = None
) -> MonoExpFitResult:
    """
    Fit mono-exponential decay model to FLIM data using Levenberg-Marquardt.
    
    Parameters:
        time_ns: Time array in nanoseconds
        intensity: Intensity array
        p0: Initial guess for [I0, tau] (default: auto-estimated)
    
    Returns:
        MonoExpFitResult with fitted parameters and R²
    """
    # Shift time to start at 0 for fitting
    t_shifted = time_ns - time_ns[0]
    
    # Auto-estimate initial parameters if not provided
    if p0 is None:
        I0_init = intensity[0]  # Initial intensity
        # Estimate tau from 1/e decay point
        threshold = I0_init / np.e
        try:
            idx_1e = np.where(intensity <= threshold)[0][0]
            tau_init = t_shifted[idx_1e]
        except IndexError:
            tau_init = t_shifted[-1] / 3  # Fallback estimate
        
        p0 = [I0_init, tau_init]
    
    # Perform fitting with Levenberg-Marquardt
    try:
        popt, pcov = curve_fit(
            mono_exponential,
            t_shifted,
            intensity,
            p0=p0,
            method='lm',  # Levenberg-Marquardt
            maxfev=10000
        )
    except RuntimeError as e:
        raise RuntimeError(f"Mono-exponential fit failed to converge: {e}")
    
    # Extract parameters and errors
    I0_fit, tau_fit = popt
    errors = np.sqrt(np.diag(pcov))
    I0_error, tau_error = errors
    
    # Calculate R²
    y_predicted = mono_exponential(t_shifted, I0_fit, tau_fit)
    r_squared = calculate_r_squared(intensity, y_predicted)
    
    return MonoExpFitResult(
        I0=I0_fit,
        I0_error=I0_error,
        tau=tau_fit,
        tau_error=tau_error,
        R_squared=r_squared,
        fit_start_ns=time_ns[0],
        fit_end_ns=time_ns[-1]
    )


def fit_bi_exponential(
    time_ns: np.ndarray,
    intensity: np.ndarray,
    p0: Optional[Tuple[float, float, float, float]] = None,
    mono_result: Optional[MonoExpFitResult] = None
) -> BiExpFitResult:
    """
    Fit bi-exponential decay model to FLIM data using Levenberg-Marquardt.
    
    Parameters:
        time_ns: Time array in nanoseconds
        intensity: Intensity array
        p0: Initial guess for [A1, tau1, A2, tau2] (default: auto-estimated)
        mono_result: Mono-exponential fit result to help initialize bi-exp guess
    
    Returns:
        BiExpFitResult with fitted parameters and R²
    """
    # Shift time to start at 0 for fitting
    t_shifted = time_ns - time_ns[0]
    
    # Auto-estimate initial parameters if not provided
    if p0 is None:
        if mono_result is not None:
            # Use mono-exponential result to initialize
            tau_mono = mono_result.tau
            I0_mono = mono_result.I0
            
            A1_init = I0_mono * 0.5
            A2_init = I0_mono * 0.5
            tau1_init = tau_mono * 0.5  # Shorter component
            tau2_init = tau_mono * 1.5  # Longer component
        else:
            # Estimate from data
            A1_init = intensity[0] * 0.5
            A2_init = intensity[0] * 0.5
            tau1_init = t_shifted[-1] / 6
            tau2_init = t_shifted[-1] / 2
        
        p0 = [A1_init, tau1_init, A2_init, tau2_init]
    
    # Set bounds to ensure positive parameters
    bounds = (
        [0, 0.001, 0, 0.001],  # Lower bounds
        [np.inf, np.inf, np.inf, np.inf]  # Upper bounds
    )
    
    # Perform fitting - use 'trf' method with bounds for stability
    try:
        popt, pcov = curve_fit(
            bi_exponential,
            t_shifted,
            intensity,
            p0=p0,
            bounds=bounds,
            method='trf',  # Trust Region Reflective (supports bounds)
            maxfev=10000
        )
    except RuntimeError as e:
        raise RuntimeError(f"Bi-exponential fit failed to converge: {e}")
    
    # Extract parameters and errors
    A1_fit, tau1_fit, A2_fit, tau2_fit = popt
    errors = np.sqrt(np.diag(pcov))
    A1_error, tau1_error, A2_error, tau2_error = errors
    
    # Ensure tau1 < tau2 for consistency (short component first)
    if tau1_fit > tau2_fit:
        A1_fit, A2_fit = A2_fit, A1_fit
        tau1_fit, tau2_fit = tau2_fit, tau1_fit
        A1_error, A2_error = A2_error, A1_error
        tau1_error, tau2_error = tau2_error, tau1_error
    
    # Calculate amplitude-weighted average lifetime
    tau_avg = (A1_fit * tau1_fit + A2_fit * tau2_fit) / (A1_fit + A2_fit)
    
    # Calculate R²
    y_predicted = bi_exponential(t_shifted, A1_fit, tau1_fit, A2_fit, tau2_fit)
    r_squared = calculate_r_squared(intensity, y_predicted)
    
    return BiExpFitResult(
        A1=A1_fit,
        A1_error=A1_error,
        tau1=tau1_fit,
        tau1_error=tau1_error,
        A2=A2_fit,
        A2_error=A2_error,
        tau2=tau2_fit,
        tau2_error=tau2_error,
        tau_avg=tau_avg,
        R_squared=r_squared,
        fit_start_ns=time_ns[0],
        fit_end_ns=time_ns[-1]
    )


# =============================================================================
# Plotting Functions
# =============================================================================

def plot_mono_exponential_fit(
    data: FLIMData,
    fit_result: MonoExpFitResult,
    time_decay: np.ndarray,
    intensity_decay: np.ndarray,
    output_path: Optional[Path] = None,
    show_plot: bool = True
) -> None:
    """
    Create publication-quality plot of mono-exponential fit.
    
    Parameters:
        data: Full FLIM data object
        fit_result: MonoExpFitResult from fitting
        time_decay: Time array used for fitting
        intensity_decay: Intensity array used for fitting
        output_path: Path to save SVG file (optional)
        show_plot: Whether to display the plot
    """
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Separate fitted vs not-fitted points for clarity
    fitted_mask = (data.time_ns >= time_decay[0]) & (data.time_ns <= time_decay[-1])
    ax.scatter(data.time_ns[~fitted_mask], data.intensity[~fitted_mask], s=15, alpha=0.5,
               color='gray', label='Not fitted', zorder=1)
    ax.scatter(data.time_ns[fitted_mask], data.intensity[fitted_mask], s=20, color='blue', 
               label='Fitted points', zorder=2)
    
    # Generate fitted curve
    t_fit = np.linspace(time_decay[0], time_decay[-1], 500)
    t_shifted = t_fit - time_decay[0]
    I_fit = mono_exponential(t_shifted, fit_result.I0, fit_result.tau)
    
    ax.plot(t_fit, I_fit, 'r-', linewidth=2, 
            label=f'Mono-exp fit: τ = {fit_result.tau:.3f} ns', zorder=3)
    
    # Mark peak position
    ax.axvline(x=data.peak_time_ns, color='green', linestyle='--', alpha=0.5,
               label=f'Peak (t = {data.peak_time_ns:.2f} ns)')
    
    # Labels and formatting (NO TITLE)
    ax.set_xlabel('Time (ns)', fontsize=12)
    ax.set_ylabel('Fluorescence Intensity (a.u.)', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    # Set axis limits to the fit window (start at peak)
    if time_decay[-1] == time_decay[0]:
        ax.set_xlim(time_decay[0], time_decay[0] + data.time_bin_width_ns)
    else:
        ax.set_xlim(time_decay[0], time_decay[-1])
    ax.set_ylim(0, np.max(data.intensity) * 1.1)
    
    # Add text box with fit parameters (outside plot area on right)
    textstr = '\n'.join([
        r'$\mathbf{Mono\text{-}exponential\ Fit}$',
        r'$I(t) = I_0 \cdot e^{-t/\tau}$',
        '',
        fr'$I_0 = {fit_result.I0:.4f} \pm {fit_result.I0_error:.4f}$',
        fr'$\tau = {fit_result.tau:.3f} \pm {fit_result.tau_error:.3f}\,\mathrm{{ns}}$',
        '',
        fr'$R^2 = {fit_result.R_squared:.6f}$'
    ])
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.95, edgecolor='gray')
    # Aligned at x=0.99, textbox at top (y=0.97)
    ax.text(0.99, 0.97, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', horizontalalignment='right', bbox=props)
    
    plt.tight_layout()
    
    # Save figure
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, format='svg', dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_path}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()


def plot_bi_exponential_fit(
    data: FLIMData,
    fit_result: BiExpFitResult,
    time_decay: np.ndarray,
    intensity_decay: np.ndarray,
    output_path: Optional[Path] = None,
    show_plot: bool = True
) -> None:
    """
    Create publication-quality plot of bi-exponential fit.
    
    Parameters:
        data: Full FLIM data object
        fit_result: BiExpFitResult from fitting
        time_decay: Time array used for fitting
        intensity_decay: Intensity array used for fitting
        output_path: Path to save SVG file (optional)
        show_plot: Whether to display the plot
    """
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Separate fitted vs not-fitted points for clarity
    fitted_mask = (data.time_ns >= time_decay[0]) & (data.time_ns <= time_decay[-1])
    ax.scatter(data.time_ns[~fitted_mask], data.intensity[~fitted_mask], s=15, alpha=0.5,
               color='gray', label='Not fitted', zorder=1)
    ax.scatter(data.time_ns[fitted_mask], data.intensity[fitted_mask], s=20, color='blue', 
               label='Fitted points', zorder=2)
    
    # Generate fitted curve
    t_fit = np.linspace(time_decay[0], time_decay[-1], 500)
    t_shifted = t_fit - time_decay[0]
    I_fit = bi_exponential(t_shifted, fit_result.A1, fit_result.tau1, 
                           fit_result.A2, fit_result.tau2)
    
    ax.plot(t_fit, I_fit, 'r-', linewidth=2, 
            label=f'Bi-exp fit: ⟨τ⟩ = {fit_result.tau_avg:.3f} ns', zorder=3)
    
    # Mark peak position
    ax.axvline(x=data.peak_time_ns, color='green', linestyle='--', alpha=0.5,
               label=f'Peak (t = {data.peak_time_ns:.2f} ns)')
    
    # Labels and formatting (NO TITLE)
    ax.set_xlabel('Time (ns)', fontsize=12)
    ax.set_ylabel('Fluorescence Intensity (a.u.)', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    # Set axis limits to the fit window (start at peak)
    if time_decay[-1] == time_decay[0]:
        ax.set_xlim(time_decay[0], time_decay[0] + data.time_bin_width_ns)
    else:
        ax.set_xlim(time_decay[0], time_decay[-1])
    ax.set_ylim(0, np.max(data.intensity) * 1.1)
    
    # Add text box with fit parameters (outside plot area on right)
    textstr = '\n'.join([
        r'$\mathbf{Bi\text{-}exponential\ Fit}$',
        r'$I(t) = A_1 e^{-t/\tau_1} + A_2 e^{-t/\tau_2}$',
        '',
        fr'$A_1 = {fit_result.A1:.4f} \pm {fit_result.A1_error:.4f}$',
        fr'$\tau_1 = {fit_result.tau1:.3f} \pm {fit_result.tau1_error:.3f}\,\mathrm{{ns}}$',
        '',
        fr'$A_2 = {fit_result.A2:.4f} \pm {fit_result.A2_error:.4f}$',
        fr'$\tau_2 = {fit_result.tau2:.3f} \pm {fit_result.tau2_error:.3f}\,\mathrm{{ns}}$',
        '',
        fr'$\langle\tau\rangle = {fit_result.tau_avg:.3f}\,\mathrm{{ns}}$',
        fr'$R^2 = {fit_result.R_squared:.6f}$'
    ])
    
    props = dict(boxstyle='round', facecolor='white', alpha=0.95, edgecolor='gray')
    # Aligned at x=0.99, textbox at top (y=0.97)
    ax.text(0.99, 0.97, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', horizontalalignment='right', bbox=props)
    
    plt.tight_layout()
    
    # Save figure
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, format='svg', dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_path}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()


# =============================================================================
# Result Printing Functions
# =============================================================================

def print_mono_result(result: MonoExpFitResult) -> None:
    """Print mono-exponential fit results."""
    print("\n" + "=" * 60)
    print("MONO-EXPONENTIAL FIT RESULTS")
    print("=" * 60)
    print(f"Model: I(t) = I₀ × exp(-t/τ)")
    print("-" * 60)
    print(f"I₀     = {result.I0:.6f} ± {result.I0_error:.6f}")
    print(f"τ      = {result.tau:.4f} ± {result.tau_error:.4f} ns")
    print("-" * 60)
    print(f"R²     = {result.R_squared:.8f}")
    print(f"Fit range: {result.fit_start_ns:.2f} - {result.fit_end_ns:.2f} ns")
    print("=" * 60)


def print_bi_result(result: BiExpFitResult) -> None:
    """Print bi-exponential fit results."""
    print("\n" + "=" * 60)
    print("BI-EXPONENTIAL FIT RESULTS")
    print("=" * 60)
    print(f"Model: I(t) = A₁×exp(-t/τ₁) + A₂×exp(-t/τ₂)")
    print("-" * 60)
    print(f"Component 1 (short lifetime):")
    print(f"  A₁   = {result.A1:.6f} ± {result.A1_error:.6f}")
    print(f"  τ₁   = {result.tau1:.4f} ± {result.tau1_error:.4f} ns")
    print(f"\nComponent 2 (long lifetime):")
    print(f"  A₂   = {result.A2:.6f} ± {result.A2_error:.6f}")
    print(f"  τ₂   = {result.tau2:.4f} ± {result.tau2_error:.4f} ns")
    print("-" * 60)
    print(f"⟨τ⟩    = {result.tau_avg:.4f} ns (amplitude-weighted average)")
    print(f"R²     = {result.R_squared:.8f}")
    print(f"Fit range: {result.fit_start_ns:.2f} - {result.fit_end_ns:.2f} ns")
    print("=" * 60)


# =============================================================================
# Main Function
# =============================================================================

def main():
    """Main function to run FLIM decay fitting."""
    parser = argparse.ArgumentParser(
        description='Fit FLIM fluorescence decay data with exponential models',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python flim_exponential_fit.py "Plot Values.csv" --model mono
    python flim_exponential_fit.py "Plot Values.csv" --model bi
    python flim_exponential_fit.py "Plot Values.csv" --model mono --fit-start 2.0
    python flim_exponential_fit.py "Plot Values.csv" --model bi --laser-rate 80

FLIM Theory:
    Fluorescence decay follows I(t) = I₀ × exp(-t/τ) for a single-component system.
    The fluorescence lifetime τ is independent of excitation intensity and provides
    information about the local molecular environment.
    
    For an 80 MHz laser: period = 12.5 ns, time bin width ≈ 0.0977 ns (128 bins)
        """
    )
    
    parser.add_argument('csv', type=str, 
                        help='Path to the CSV file containing FLIM data')
    
    parser.add_argument('--model', type=str, choices=['mono', 'bi'], default='mono',
                        help='Decay model: mono (mono-exponential) or bi (bi-exponential). Default: mono')
    
    parser.add_argument('--laser-rate', type=float, default=DEFAULT_LASER_REP_RATE_MHZ,
                        help=f'Laser repetition rate in MHz (default: {DEFAULT_LASER_REP_RATE_MHZ})')
    
    parser.add_argument('--fit-start', type=float, default=None,
                        help='Manual start time for fit (ns). If omitted, auto-detect peak.')
    
    parser.add_argument('--fit-end', type=float, default=None,
                        help='Manual end time for fit (ns). If omitted, use all data after start.')
    
    parser.add_argument('--output-dir', type=str, default=None,
                        help=f'Output directory for plots (default: {OUTPUT_DIR})')
    
    parser.add_argument('--no-show', action='store_true',
                        help='Do not display plot (just save to file)')
    
    parser.add_argument('--output-name', type=str, default=None,
                        help='Custom output filename prefix (without extension)')
    
    parser.add_argument('--x-column', type=str, default='X',
                        help='Name of the X column in CSV (default: X)')
    
    parser.add_argument('--y-column', type=str, default='Y',
                        help='Name of the Y column in CSV (default: Y)')
    
    args = parser.parse_args()
    
    # Set output directory
    output_dir = Path(args.output_dir) if args.output_dir else OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read FLIM data
    print(f"\nReading FLIM data from: {args.csv}")
    print(f"Laser repetition rate: {args.laser_rate} MHz")
    
    data = read_flim_csv(args.csv, laser_rep_rate_mhz=args.laser_rate, 
                         x_column=args.x_column, y_column=args.y_column)
    print(get_data_summary(data))
    
    # Generate base name from input file for output naming
    input_stem = Path(args.csv).stem  # e.g., "Plot Values" -> "Plot Values"
    # Clean the stem for filename (replace spaces with underscores)
    input_stem_clean = input_stem.replace(' ', '_')
    
    # Save time-intensity CSV for reference (helps choosing fit ranges)
    csv_output_path = output_dir / f"{input_stem_clean}_time_intensity.csv"
    time_intensity_df = pd.DataFrame({
        'time_ns': data.time_ns,
        'intensity': data.intensity
    })
    time_intensity_df.to_csv(csv_output_path, index=False)
    print(f"\nTime-intensity data saved to: {csv_output_path}")
    
    # Determine fit range
    if args.fit_start is not None:
        # Find index closest to specified start time
        start_idx = np.argmin(np.abs(data.time_ns - args.fit_start))
        auto_used = False
    else:
        start_idx, _ = auto_fit_indices_from_peak_to_last_nonzero(data)
        auto_used = True
        print(f"\nAuto-detected peak at bin {start_idx} (t = {data.peak_time_ns:.2f} ns)")

    if args.fit_end is not None:
        end_idx = np.argmin(np.abs(data.time_ns - args.fit_end))
    else:
        # Auto end: last non-zero after peak (inclusive -> exclusive slice)
        _, end_idx = auto_fit_indices_from_peak_to_last_nonzero(data)

    # Extract decay region
    time_decay, intensity_decay = extract_decay_region(data, start_idx, end_idx)
    if auto_used:
        print(
            f"Auto-selected decay window: peak at {time_decay[0]:.2f} ns "
            f"to last non-zero at {time_decay[-1]:.2f} ns ({len(time_decay)} points)"
        )
    else:
        print(f"Fitting decay region: {time_decay[0]:.2f} - {time_decay[-1]:.2f} ns ({len(time_decay)} points)")
    
    # Generate descriptive output filename
    # Format: {input_name}_{model}_fit_t{start}-{end}ns.svg
    fit_range_str = f"t{time_decay[0]:.1f}-{time_decay[-1]:.1f}ns"
    
    # Perform fitting based on model choice
    if args.model == 'mono':
        print("\n>>> Fitting MONO-EXPONENTIAL model...")
        result = fit_mono_exponential(time_decay, intensity_decay)
        print_mono_result(result)
        
        # Generate output filename
        if args.output_name:
            output_file = output_dir / f"{args.output_name}_mono_{fit_range_str}.svg"
        else:
            output_file = output_dir / f"{input_stem_clean}_mono_{fit_range_str}.svg"
        
        # Plot
        plot_mono_exponential_fit(
            data, result, time_decay, intensity_decay,
            output_path=output_file,
            show_plot=not args.no_show
        )
        
    else:  # bi-exponential
        # First fit mono to get initial estimates
        print("\n>>> Fitting mono-exponential first for initial estimates...")
        mono_result = fit_mono_exponential(time_decay, intensity_decay)
        print(f"    Mono-exp τ = {mono_result.tau:.3f} ns (used for bi-exp initialization)")
        
        print("\n>>> Fitting BI-EXPONENTIAL model...")
        result = fit_bi_exponential(time_decay, intensity_decay, mono_result=mono_result)
        print_bi_result(result)
        
        # Generate output filename
        if args.output_name:
            output_file = output_dir / f"{args.output_name}_biexp_{fit_range_str}.svg"
        else:
            output_file = output_dir / f"{input_stem_clean}_biexp_{fit_range_str}.svg"
        
        # Plot
        plot_bi_exponential_fit(
            data, result, time_decay, intensity_decay,
            output_path=output_file,
            show_plot=not args.no_show
        )
    
    print("\nDone!")


if __name__ == "__main__":
    main()
