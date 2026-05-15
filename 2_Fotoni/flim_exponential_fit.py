"""
FLIM Exponential Decay Fitting Library
=======================================

Core library module for FLIM (Fluorescence Lifetime Imaging Microscopy) analysis.
Provides data structures, fitting functions, and plotting utilities used by the
TIFF-based pipeline in ``flim_tiff_fit.py``.

Data structures:
    FLIMData               — container for a 1-D decay trace with time metadata

Fitting functions:
    fit_mono_exponential   — mono-exp: I(t) = I₀ × exp(-t/τ)
    fit_bi_exponential     — bi-exp: I(t) = A₁×exp(-t/τ₁) + A₂×exp(-t/τ₂)

Fit uses the Levenberg-Marquardt algorithm (scipy.optimize.curve_fit).
Goodness-of-fit: χ²_red (Poisson weights) when parameter errors are available;
falls back to R² when the covariance matrix is singular.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Tuple


# Default TCSPC laser repetition rate
DEFAULT_LASER_REP_RATE_MHZ = 80.0  # MHz


@dataclass(frozen=True)
class FLIMData:
    """
    Container for a 1-D FLIM decay trace with time-axis metadata.

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


def extract_decay_region(
    data: FLIMData,
    start_index: Optional[int] = None,
    end_index: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Slice the decay region from FLIM data.

    Parameters:
        data: FLIMData object
        start_index: First bin (default: peak index)
        end_index: Last bin, exclusive (default: end of array)

    Returns:
        (time_ns, intensity) arrays for the selected region
    """
    if start_index is None:
        start_index = data.peak_index
    if end_index is None:
        end_index = len(data.intensity)
    start_index = max(0, start_index)
    end_index = min(len(data.intensity), end_index)
    return data.time_ns[start_index:end_index], data.intensity[start_index:end_index]


def get_data_summary(data: FLIMData) -> str:
    """Return a formatted one-line-per-field summary of the FLIM data."""
    return (
        f"FLIM Data Summary\n"
        f"=================\n"
        f"Number of time bins: {data.num_bins}\n"
        f"Time bin width: {data.time_bin_width_ns:.4f} ns\n"
        f"Total time range: {data.time_ns[-1]:.2f} ns\n"
        f"Peak intensity: {np.max(data.intensity):.4f}\n"
        f"Peak position: bin {data.peak_index} ({data.peak_time_ns:.2f} ns)\n"
        f"Intensity range: [{np.min(data.intensity):.4f}, {np.max(data.intensity):.4f}]"
    )


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
        chi2_red: Reduced chi-squared (Poisson weights; ≈1 for good fit).
                  NaN when parameter errors are not available.
        r2: Coefficient of determination R².  Computed only when chi2_red
            is NaN (i.e. errors unavailable); None otherwise.
        fit_start_ns: Start time of fitted region (ns)
        fit_end_ns: End time of fitted region (ns)
    """
    I0: float
    I0_error: float
    tau: float
    tau_error: float
    chi2_red: float
    fit_start_ns: float
    fit_end_ns: float
    r2: Optional[float] = None


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
        chi2_red: Reduced chi-squared (Poisson weights; ≈1 for good fit).
                  NaN when parameter errors are not available.
        r2: Coefficient of determination R².  Computed only when chi2_red
            is NaN (i.e. errors unavailable); None otherwise.
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
    chi2_red: float
    fit_start_ns: float
    fit_end_ns: float
    r2: Optional[float] = None


# =============================================================================
# Fitting Functions
# =============================================================================

def calculate_chi2_reduced(
    counts: np.ndarray,
    y_predicted: np.ndarray,
    n_params: int,
) -> float:
    """
    Compute reduced chi-squared using Poisson statistics.

    σ_i = sqrt(max(counts_i, 1))  (floor at 1 avoids div-by-zero)
    χ²_red = (1/ν) × Σ [(C_i - f_i)² / σ_i²]
    where ν = n_points - n_params is the number of degrees of freedom.

    A value ≈ 1 indicates a good fit; >> 1 means the model is inadequate
    or errors are underestimated; << 1 suggests over-fitting or
    over-estimated errors.
    """
    sigma = np.sqrt(np.maximum(counts, 1.0))
    chi2 = float(np.sum(((counts - y_predicted) / sigma) ** 2))
    dof = len(counts) - n_params
    return chi2 / dof if dof > 0 else np.inf


def calculate_r2(y_observed: np.ndarray, y_predicted: np.ndarray) -> float:
    """
    Compute the coefficient of determination R².

    R² = 1 - SS_res / SS_tot

    A value of 1 indicates a perfect fit; values close to 1 are good.
    Used as a fallback goodness-of-fit metric when Poisson-weighted
    χ²_red is not meaningful (e.g. covariance matrix is singular).
    """
    ss_res = float(np.sum((y_observed - y_predicted) ** 2))
    ss_tot = float(np.sum((y_observed - np.mean(y_observed)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 0.0


def auto_fit_indices_from_peak_to_last_nonzero(
    data: FLIMData,
    min_counts: int = 1,
) -> Tuple[int, int]:
    """
    Pick fit window starting at the peak and ending at the last bin with
    sufficient counts (>= min_counts).  Use min_counts=5 (or 10) when
    intensity is in raw photon counts to exclude the low-statistics tail
    where the Gaussian approximation to Poisson breaks down.

    Parameters:
        data: FLIMData containing the full trace
        min_counts: Minimum bin value to include in the fit window

    Returns:
        (start_idx, end_idx) suitable for slicing (end is exclusive)
    """
    start_idx = data.peak_index

    # Look for the last sample at or after the peak with enough counts.
    after_peak_nonzero = np.where(data.intensity[start_idx:] >= min_counts)[0]

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


def _closest_time_index(time_ns: np.ndarray, target_ns: float) -> int:
    """Return index of the sample closest to target time (ns)."""
    return int(np.argmin(np.abs(time_ns - target_ns)))


def resolve_fit_indices(
    data: FLIMData,
    fit_start_ns: Optional[float],
    fit_end_ns: Optional[float],
    min_counts: int = 1,
) -> Tuple[int, int, bool]:
    """
    Resolve fit slice indices from optional manual times.

    Parameters:
        min_counts: Passed to auto_fit_indices_from_peak_to_last_nonzero when
                    fit_end_ns is None.  Set to 5 when intensity is raw photon
                    counts (reduce='sum') to exclude the low-statistics tail.

    Returns:
        (start_idx, end_idx, auto_used), where end_idx is exclusive.
    """
    auto_start_idx, auto_end_idx = auto_fit_indices_from_peak_to_last_nonzero(
        data, min_counts=min_counts
    )

    if fit_start_ns is None:
        start_idx = auto_start_idx
        auto_used = True
    else:
        start_idx = _closest_time_index(data.time_ns, fit_start_ns)
        auto_used = False

    if fit_end_ns is None:
        end_idx = auto_end_idx
    else:
        end_idx = _closest_time_index(data.time_ns, fit_end_ns)

    # Keep behavior stable and avoid empty/invalid slices.
    if end_idx <= start_idx:
        end_idx = min(len(data.intensity), start_idx + 2)

    return start_idx, end_idx, auto_used


def fit_mono_exponential(
    time_ns: np.ndarray,
    intensity: np.ndarray,
    p0: Optional[Tuple[float, float]] = None,
    counts: Optional[np.ndarray] = None,
) -> MonoExpFitResult:
    """
    Fit mono-exponential decay model to FLIM data using Levenberg-Marquardt.

    Parameters:
        time_ns: Time array in nanoseconds
        intensity: Intensity array (photon counts per bin when counts is None)
        p0: Initial guess for [I0, tau] (default: auto-estimated)
        counts: Raw photon counts per bin (same shape as intensity).  When
                provided, Poisson weights σ_i = sqrt(max(N_i, 1)) are passed
                to curve_fit and χ²_red is computed from them.  If None,
                intensity is used as a proxy for counts.

    Returns:
        MonoExpFitResult with fitted parameters and χ²_red
    """
    # Shift time to start at 0 for fitting
    t_shifted = time_ns - time_ns[0]
    
    # Auto-estimate initial parameters if not provided
    if p0 is None:
        I0_init = intensity[0]  # Initial intensity
        # Estimate tau from 1/e decay point (I(t) = I0/e) or fallback to 1/3 of the time range
        threshold = I0_init / np.e
        try:
            idx_1e = np.where(intensity <= threshold)[0][0]
            tau_init = t_shifted[idx_1e]
        except IndexError:
            tau_init = t_shifted[-1] / 3  # Fallback estimate
        
        p0 = [I0_init, tau_init]
    
    # Poisson weights: σ_i = sqrt(max(N_i, 1))
    _counts = counts if counts is not None else intensity
    sigma = np.sqrt(np.maximum(_counts, 1.0))

    # Perform fitting with Levenberg-Marquardt
    try:
        popt, pcov = curve_fit(
            mono_exponential,
            t_shifted,
            intensity,
            p0=p0,
            sigma=sigma,
            absolute_sigma=True,
            method='lm',  # Levenberg-Marquardt
            maxfev=10000
        )
    except RuntimeError as e:
        raise RuntimeError(f"Mono-exponential fit failed to converge: {e}")

    # Extract parameters and errors
    I0_fit, tau_fit = popt
    diag = np.diag(pcov)
    errors_available = not np.any(np.isinf(diag))
    errors = np.sqrt(np.abs(diag))  # abs guards against tiny negatives from numerics
    I0_error, tau_error = errors

    # Goodness-of-fit: χ²_red when errors are reliable, R² otherwise
    y_predicted = mono_exponential(t_shifted, I0_fit, tau_fit)
    if errors_available:
        chi2_red = calculate_chi2_reduced(_counts, y_predicted, n_params=2)
        r2 = None
    else:
        chi2_red = np.nan
        r2 = calculate_r2(intensity, y_predicted)

    return MonoExpFitResult(
        I0=I0_fit,
        I0_error=I0_error,
        tau=tau_fit,
        tau_error=tau_error,
        chi2_red=chi2_red,
        fit_start_ns=time_ns[0],
        fit_end_ns=time_ns[-1],
        r2=r2
    )


def fit_bi_exponential(
    time_ns: np.ndarray,
    intensity: np.ndarray,
    p0: Optional[Tuple[float, float, float, float]] = None,
    mono_result: Optional[MonoExpFitResult] = None,
    counts: Optional[np.ndarray] = None,
) -> BiExpFitResult:
    """
    Fit bi-exponential decay model to FLIM data using Levenberg-Marquardt.

    Parameters:
        time_ns: Time array in nanoseconds
        intensity: Intensity array (photon counts per bin when counts is None)
        p0: Initial guess for [A1, tau1, A2, tau2] (default: auto-estimated)
        mono_result: Mono-exponential fit result to help initialize bi-exp guess
        counts: Raw photon counts per bin (same shape as intensity).  When
                provided, Poisson weights σ_i = sqrt(max(N_i, 1)) are passed
                to curve_fit and χ²_red is computed from them.

    Returns:
        BiExpFitResult with fitted parameters and χ²_red
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
    
    # Poisson weights: σ_i = sqrt(max(N_i, 1))
    _counts = counts if counts is not None else intensity
    sigma = np.sqrt(np.maximum(_counts, 1.0))

    # Set bounds to keep parameters within physically meaningful ranges.
    # τ ∈ [0.05, 12] ns: lower bound avoids sub-resolution artefacts;
    # upper bound is the laser period (80 MHz → 12.5 ns).
    # Without tight bounds the optimizer often finds degenerate solutions
    # (e.g. τ₁ ≈ 0 or τ ≫ 12 ns) that are then rejected as non-physical.
    bounds = (
        [0, 0.05, 0, 0.05],           # A ≥ 0, τ ≥ 0.05 ns
        [np.inf, 12.0, np.inf, 12.0]  # A unlimited, τ ≤ 12 ns
    )

    # Perform fitting - use 'trf' method with bounds for stability
    try:
        popt, pcov = curve_fit(
            bi_exponential,
            t_shifted,
            intensity,
            p0=p0,
            sigma=sigma,
            absolute_sigma=True,
            bounds=bounds,
            method='trf',  # Trust Region Reflective (supports bounds)
            maxfev=10000
        )
    except RuntimeError as e:
        raise RuntimeError(f"Bi-exponential fit failed to converge: {e}")

    # Extract parameters and errors
    A1_fit, tau1_fit, A2_fit, tau2_fit = popt
    diag = np.diag(pcov)
    errors_available = not np.any(np.isinf(diag))
    errors = np.sqrt(np.abs(diag))
    A1_error, tau1_error, A2_error, tau2_error = errors

    # Ensure tau1 < tau2 for consistency (short component first)
    if tau1_fit > tau2_fit:
        A1_fit, A2_fit = A2_fit, A1_fit
        tau1_fit, tau2_fit = tau2_fit, tau1_fit
        A1_error, A2_error = A2_error, A1_error
        tau1_error, tau2_error = tau2_error, tau1_error

    # Calculate amplitude-weighted average lifetime
    tau_avg = (A1_fit * tau1_fit + A2_fit * tau2_fit) / (A1_fit + A2_fit)

    # Goodness-of-fit: χ²_red when errors are reliable, R² otherwise
    y_predicted = bi_exponential(t_shifted, A1_fit, tau1_fit, A2_fit, tau2_fit)
    if errors_available:
        chi2_red = calculate_chi2_reduced(_counts, y_predicted, n_params=4)
        r2 = None
    else:
        chi2_red = np.nan
        r2 = calculate_r2(intensity, y_predicted)

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
        chi2_red=chi2_red,
        fit_start_ns=time_ns[0],
        fit_end_ns=time_ns[-1],
        r2=r2
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

    # Poisson error bars on fitted points: σ_i = √max(N_i, 1)
    sigma_fitted = np.sqrt(np.maximum(data.intensity[fitted_mask], 1.0))
    ax.errorbar(data.time_ns[fitted_mask], data.intensity[fitted_mask],
                yerr=sigma_fitted, fmt='o', ms=3, color='blue', ecolor='royalblue',
                elinewidth=0.8, capsize=2, alpha=0.7, label='Fitted points ± √N', zorder=2)

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
    if fit_result.r2 is not None:
        gof_line = f'$R^2$ = {fit_result.r2:.6f}'
    else:
        gof_line = f'$\\chi^2_{{\\mathrm{{red}}}}$ = {fit_result.chi2_red:.4f}'
    textstr = '\n'.join([
        r'$\mathbf{Mono-exponential\ Fit}$',
        r'$I(t) = I_0 \cdot e^{-t/\tau}$',
        '',
        f'$I_0$ = {fit_result.I0:.4f} ± {fit_result.I0_error:.4f}',
        f'$\\tau$ = {fit_result.tau:.3f} ± {fit_result.tau_error:.3f} ns',
        '',
        gof_line
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

    # Poisson error bars on fitted points: σ_i = √max(N_i, 1)
    sigma_fitted = np.sqrt(np.maximum(data.intensity[fitted_mask], 1.0))
    ax.errorbar(data.time_ns[fitted_mask], data.intensity[fitted_mask],
                yerr=sigma_fitted, fmt='o', ms=3, color='blue', ecolor='royalblue',
                elinewidth=0.8, capsize=2, alpha=0.7, label='Fitted points ± √N', zorder=2)

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
    if fit_result.r2 is not None:
        gof_line = f'$R^2$ = {fit_result.r2:.6f}'
    else:
        gof_line = f'$\\chi^2_{{\\mathrm{{red}}}}$ = {fit_result.chi2_red:.4f}'
    textstr = '\n'.join([
        r'$\mathbf{Bi-exponential\ Fit}$',
        r'$I(t) = A_1 e^{-t/\tau_1} + A_2 e^{-t/\tau_2}$',
        '',
        f'$A_1$ = {fit_result.A1:.4f} ± {fit_result.A1_error:.4f}',
        f'$\\tau_1$ = {fit_result.tau1:.3f} ± {fit_result.tau1_error:.3f} ns',
        '',
        f'$A_2$ = {fit_result.A2:.4f} ± {fit_result.A2_error:.4f}',
        f'$\\tau_2$ = {fit_result.tau2:.3f} ± {fit_result.tau2_error:.3f} ns',
        '',
        f'$\\langle\\tau\\rangle$ = {fit_result.tau_avg:.3f} ns',
        gof_line
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
    if result.r2 is not None:
        print("(Errors unavailable — R² used instead of χ²_red)")
        print(f"R²     = {result.r2:.6f}")
    else:
        print(f"χ²_red = {result.chi2_red:.4f}")
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
    if result.r2 is not None:
        print("(Errors unavailable — R² used instead of χ²_red)")
        print(f"R²     = {result.r2:.6f}")
    else:
        print(f"χ²_red = {result.chi2_red:.4f}")
    print(f"Fit range: {result.fit_start_ns:.2f} - {result.fit_end_ns:.2f} ns")
    print("=" * 60)


