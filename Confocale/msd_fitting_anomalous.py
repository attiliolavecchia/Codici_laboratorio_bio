"""
Anomalous diffusion MSD fitting module (non-drifting model).

This module implements fitting to the anomalous diffusion model without drift:
    
    MSD(τ) = 4D_α·τ^α
    
where:
    - MSD is the Mean Squared Displacement
    - τ (tau) is the time lag
    - D_α is the generalized diffusion coefficient
    - α (alpha) is the anomalous diffusion exponent
    
The exponent α characterizes the type of diffusion:
    - α < 1: Subdiffusion (confined/hindered motion, crowded environments)
    - α = 1: Normal Brownian diffusion
    - α > 1: Superdiffusion (active transport, ballistic components)
    - α = 2: Pure ballistic motion

The fitting is performed over multiple time intervals (10%-90% of data) to find
the optimal range based on R² and residual analysis.

Physical Context:
    For particles in viscoelastic media or crowded environments, subdiffusion
    (α < 1) is commonly observed. The generalized diffusion coefficient D_α
    depends on α and has units [length²/time^α].

Dependencies: numpy, scipy, matplotlib (optional)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Mapping, Tuple, Union

import numpy as np
from scipy.optimize import curve_fit

from data_reader import Trajectory, read_trajectories_from_csv
from msd_analyzer import calculate_ensemble_msd
from msd_fitting import calculate_r_squared


@dataclass(frozen=True)
class AnomalousFitResult:
    """Container for anomalous MSD fit results (non-drifting model).
    
    Attributes:
        D_alpha: Generalized diffusion coefficient [μm²/s^α]
        D_alpha_error: Standard error on D_α
        alpha: Anomalous diffusion exponent (dimensionless)
        alpha_error: Standard error on α
        pcov: Full covariance matrix from the fit
        tau_fit: Time lag values used in the fit
        msd_fit: MSD values used in the fit
        msd_predicted: Predicted MSD values from the fitted model
        R_squared: Coefficient of determination (R²)
        RSS: Residual sum of squares
        optimal_fraction: Optimal fraction of data used (0.1-0.9)
        n_fit_steps: Number of lag steps used for fitting
        interval_results: Dictionary mapping fraction to (R², RSS) for all tested intervals
    """
    
    D_alpha: float
    D_alpha_error: float
    alpha: float
    alpha_error: float
    pcov: np.ndarray
    tau_fit: np.ndarray
    msd_fit: np.ndarray
    msd_predicted: np.ndarray
    R_squared: float
    RSS: float
    optimal_fraction: float
    n_fit_steps: int
    interval_results: Dict[float, Tuple[float, float]]


def anomalous_msd_model(tau: np.ndarray, D_alpha: float, alpha: float) -> np.ndarray:
    """Anomalous MSD model: MSD(τ) = 4D_α·τ^α
    
    Captures anomalous diffusion behavior through the exponent α.
    
    Args:
        tau: Time lag values [s]
        D_alpha: Generalized diffusion coefficient [μm²/s^α]
        alpha: Anomalous diffusion exponent (dimensionless)
    
    Returns:
        MSD values predicted by the model [μm²]
    """
    return 4.0 * D_alpha * np.power(tau, alpha)


def calculate_rss(y_observed: np.ndarray, y_predicted: np.ndarray) -> float:
    """Calculate residual sum of squares (RSS).
    
    RSS = Σ(y_observed - y_predicted)²
    
    Args:
        y_observed: Observed data values
        y_predicted: Predicted values from the model
    
    Returns:
        RSS value (lower is better)
    """
    return float(np.sum((y_observed - y_predicted) ** 2))


def fit_msd_anomalous_at_fraction(
    tau: np.ndarray,
    msd: np.ndarray,
    n_max: int,
    dt: float,
    fit_fraction: float,
    D_alpha_initial: float,
    D_alpha_bounds: Tuple[float, float],
    alpha_initial: float,
    alpha_bounds: Tuple[float, float],
) -> Tuple[float, float, float, float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]:
    """Fit MSD to anomalous diffusion model at a specific fraction of data.
    
    Returns:
        Tuple of (D_alpha, D_alpha_error, alpha, alpha_error, pcov, 
                  tau_fit, msd_fit, msd_predicted, R_squared, RSS)
    """
    # Calculate number of steps for this fraction
    n_fit_steps = max(2, int(fit_fraction * n_max))
    tau_max = n_fit_steps * dt
    
    # Select data subset
    mask = tau <= tau_max
    tau_fit = tau[mask]
    msd_fit = msd[mask]
    
    # Remove NaN/inf
    valid = np.isfinite(tau_fit) & np.isfinite(msd_fit) & (tau_fit > 0)
    tau_fit = tau_fit[valid]
    msd_fit = msd_fit[valid]
    
    if tau_fit.size < 3:  # Need at least 3 points for 2-parameter fit
        raise ValueError(f"Insufficient points ({tau_fit.size}) for anomalous fit at fraction {fit_fraction}")
    
    # Perform the fit
    try:
        popt, pcov = curve_fit(
            anomalous_msd_model,
            tau_fit,
            msd_fit,
            p0=[D_alpha_initial, alpha_initial],
            bounds=([D_alpha_bounds[0], alpha_bounds[0]], [D_alpha_bounds[1], alpha_bounds[1]]),
            method='trf',
            maxfev=5000,
        )
    except RuntimeError as e:
        raise RuntimeError(f"Fitting failed at fraction {fit_fraction}: {e}")
    
    D_alpha_opt = float(popt[0])
    alpha_opt = float(popt[1])
    
    # Calculate errors
    if pcov is not None and np.all(np.isfinite(pcov)):
        D_alpha_error = float(np.sqrt(pcov[0, 0]))
        alpha_error = float(np.sqrt(pcov[1, 1]))
    else:
        D_alpha_error = float('nan')
        alpha_error = float('nan')
    
    # Calculate metrics
    msd_predicted = anomalous_msd_model(tau_fit, D_alpha_opt, alpha_opt)
    R_squared = calculate_r_squared(msd_fit, msd_predicted)
    RSS = calculate_rss(msd_fit, msd_predicted)
    
    return D_alpha_opt, D_alpha_error, alpha_opt, alpha_error, pcov, tau_fit, msd_fit, msd_predicted, R_squared, RSS


def fit_msd_anomalous(
    tau: np.ndarray,
    msd: np.ndarray,
    n_max: int,
    dt: float,
    D_alpha_initial: float = 1e-2,
    D_alpha_bounds: Tuple[float, float] = (1e-6, 1e2),
    alpha_initial: float = 1.0,
    alpha_bounds: Tuple[float, float] = (0.01, 2.0),
    interval_step: float = 0.10,
) -> AnomalousFitResult:
    """Fit MSD to anomalous diffusion model testing multiple intervals.
    
    Tests intervals from 10% to 90% in steps of interval_step (default 10%).
    Selects optimal interval based on:
    1. Maximum R²
    2. If multiple intervals have R² within 0.01, select minimum RSS
    
    Args:
        tau: Array of time lag values [s]
        msd: Array of MSD values [μm²]
        n_max: Maximum number of lag steps available
        dt: Time step [s]
        D_alpha_initial: Initial guess for D_α [μm²/s^α]
        D_alpha_bounds: Tuple (lower, upper) bounds for D_α [μm²/s^α]
        alpha_initial: Initial guess for α (default: 1.0 = normal diffusion)
        alpha_bounds: Tuple (lower, upper) bounds for α (default: 0.01-2.0)
        interval_step: Step size for interval search (default: 0.10 = 10%)
    
    Returns:
        AnomalousFitResult with optimal fit parameters and all interval results
    """
    # Generate fractions to test: 0.1, 0.2, ..., 0.9
    fractions = np.arange(interval_step, 1.0, interval_step)
    
    # Store results for each fraction
    results: Dict[float, Tuple[float, float, float, float, float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]] = {}
    interval_results: Dict[float, Tuple[float, float]] = {}
    
    print(f"\nTesting intervals from {fractions[0]:.0%} to {fractions[-1]:.0%} in {interval_step:.0%} steps...")
    
    for frac in fractions:
        try:
            D_alpha, D_alpha_err, alpha, alpha_err, pcov, tau_f, msd_f, msd_p, R2, rss = fit_msd_anomalous_at_fraction(
                tau, msd, n_max, dt, frac,
                D_alpha_initial, D_alpha_bounds,
                alpha_initial, alpha_bounds
            )
            results[frac] = (D_alpha, D_alpha_err, alpha, alpha_err, R2, pcov, tau_f, msd_f, msd_p, rss, int(frac * n_max))
            interval_results[frac] = (R2, rss)
            print(f"  {frac:4.0%}: R^2 = {R2:.6f}, RSS = {rss:.6e}, D_alpha = {D_alpha:.6e}, alpha = {alpha:.4f}")
        except (ValueError, RuntimeError) as e:
            print(f"  {frac:4.0%}: Failed - {e}")
            continue
    
    if len(results) == 0:
        raise RuntimeError("No valid fits obtained across any interval")
    
    # Find optimal fraction
    # Step 1: Find maximum R²
    max_R2 = max(r[4] for r in results.values())
    
    # Step 2: Find all fractions within 0.01 of max R²
    candidates = [(frac, r) for frac, r in results.items() if r[4] >= max_R2 - 0.01]
    
    # Step 3: Among candidates, select minimum RSS
    optimal_fraction = min(candidates, key=lambda x: x[1][9])[0]
    
    # Extract optimal result
    D_alpha, D_alpha_err, alpha, alpha_err, R2, pcov, tau_f, msd_f, msd_p, rss, n_steps = results[optimal_fraction]
    
    print(f"\nOptimal interval: {optimal_fraction:.0%} (R^2 = {R2:.6f}, RSS = {rss:.6e})")
    
    return AnomalousFitResult(
        D_alpha=D_alpha,
        D_alpha_error=D_alpha_err,
        alpha=alpha,
        alpha_error=alpha_err,
        pcov=pcov,
        tau_fit=tau_f,
        msd_fit=msd_f,
        msd_predicted=msd_p,
        R_squared=R2,
        RSS=rss,
        optimal_fraction=optimal_fraction,
        n_fit_steps=n_steps,
        interval_results=interval_results,
    )


def plot_anomalous_fit(
    tau_fit: np.ndarray,
    msd_fit: np.ndarray,
    msd_predicted: np.ndarray,
    D_alpha: float,
    D_alpha_error: float,
    alpha: float,
    alpha_error: float,
    R_squared: float,
    output_path: Path,
    optimal_fraction: float = 0.10,
    n_fit_steps: int = 0,
) -> None:
    """Create and save a plot of MSD data with anomalous diffusion fit.
    
    Plots data points and fit line for anomalous model: MSD = 4D_α·τ^α
    
    Args:
        tau_fit: Time lag values used in the fit
        msd_fit: MSD data values used in the fit
        msd_predicted: Predicted MSD values from the fitted model
        D_alpha: Fitted generalized diffusion coefficient [μm²/s^α]
        D_alpha_error: Standard error on D_α
        alpha: Fitted anomalous exponent
        alpha_error: Standard error on α
        R_squared: Coefficient of determination
        output_path: Path to save the figure
        optimal_fraction: Optimal fraction of data used for fitting
        n_fit_steps: Number of lag steps used for fitting
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Warning: matplotlib not available, skipping plot generation")
        return
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot data points used in the fit
    ax.plot(tau_fit, msd_fit, 'o', color='C0', markersize=8, 
            alpha=0.7, label='MSD Data', zorder=2)
    
    # Plot the fitted line
    ax.plot(tau_fit, msd_predicted, '-', color='C3', linewidth=2.5,
            label=r'Anomalous Fit: MSD = 4$D_\alpha \tau^\alpha$', zorder=3)
    
    # Labels
    ax.set_xlabel(r'Time Lag $\tau$ [s]', fontsize=12)
    ax.set_ylabel(r'MSD [μm$^2$]', fontsize=12)
    
    # Grid
    ax.grid(True, linestyle=':', alpha=0.4)
    
    # Legend
    ax.legend(loc='upper left', fontsize=10, framealpha=0.9)
    
    # Text box with fit results
    textstr = '\n'.join([
        r'$D_\alpha = (%.2e \pm %.1e)\ \mu m^2/s^\alpha$' % (D_alpha, D_alpha_error),
        r'$\alpha = (%.4f \pm %.4f)$' % (alpha, alpha_error),
        r'$R^2 = %.6f$' % R_squared,
    ])
    
    # Place text box in upper right with white background
    props = dict(boxstyle='round', facecolor='white', alpha=0.95, edgecolor='black', linewidth=1.2)
    ax.text(0.98, 0.97, textstr, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', horizontalalignment='right',
            bbox=props)
    
    # Tight layout and save
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plot saved to: {output_path.resolve()}")


def main() -> None:
    """Main function for anomalous diffusion MSD fitting analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Fit MSD data to anomalous diffusion model: MSD = 4D_α·τ^α',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    # Automatic optimization (tests 10%%-90%%)
    python msd_fitting_anomalous.py Data/Experiment8_spots_40minstep.csv
    
    # Manual interval selection (fit at 30%%)
    python msd_fitting_anomalous.py Data/Experiment8_spots_40minstep.csv --manual-fraction 0.3
    
    # Custom alpha bounds for subdiffusion
    python msd_fitting_anomalous.py Data/Experiment8_spots_40minstep.csv --alpha-bounds 0.1 0.99
    
    # Compare with normal diffusion model (α=1)
    python msd_fitting_anomalous.py Data/Experiment8_spots_40minstep.csv --compare-normal
        """
    )
    
    parser.add_argument('csv', type=str, help='Path to the trajectories CSV file')
    
    parser.add_argument(
        '--max-lag-fraction',
        type=float,
        default=None,
        help='Fraction (0<f<=1] of longest track to cap max lag for MSD calculation'
    )
    
    parser.add_argument(
        '--D-alpha-initial',
        type=float,
        default=1e-2,
        help='Initial guess for D_α [μm²/s^α] (default: 1e-2)'
    )
    
    parser.add_argument(
        '--D-alpha-bounds',
        type=float,
        nargs=2,
        default=[1e-6, 1e2],
        metavar=('LOWER', 'UPPER'),
        help='Bounds for D_α [μm²/s^α] (default: 1e-6 1e2)'
    )
    
    parser.add_argument(
        '--alpha-initial',
        type=float,
        default=1.0,
        help='Initial guess for α (default: 1.0 = normal diffusion)'
    )
    
    parser.add_argument(
        '--alpha-bounds',
        type=float,
        nargs=2,
        default=[0.01, 2.0],
        metavar=('LOWER', 'UPPER'),
        help='Bounds for α (default: 0.01 2.0, covering subdiffusion to ballistic)'
    )
    
    parser.add_argument(
        '--interval-step',
        type=float,
        default=0.10,
        help='Step size for interval search (default: 0.10 = 10%%)'
    )
    
    parser.add_argument(
        '--manual-fraction',
        type=float,
        default=None,
        metavar='FRACTION',
        help='Manually specify fit interval fraction (0.1-0.9). Skips automatic optimization.'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        default='anomalous_fits',
        help='Output directory for results (default: anomalous_fits/)'
    )
    
    parser.add_argument(
        '--compare-normal',
        action='store_true',
        help='Also run fit with α fixed at 1.0 (normal diffusion) for comparison'
    )
    
    args = parser.parse_args()
    
    # Validate manual-fraction if provided
    if args.manual_fraction is not None:
        if not (0.1 <= args.manual_fraction <= 0.9):
            print(f"Error: --manual-fraction must be between 0.1 and 0.9 (got {args.manual_fraction})")
            return
    
    # Validate alpha bounds
    if args.alpha_bounds[0] >= args.alpha_bounds[1]:
        print(f"Error: alpha lower bound must be less than upper bound")
        return
    if args.alpha_bounds[0] <= 0:
        print(f"Error: alpha lower bound must be positive")
        return
    
    # Create output directory relative to script location
    script_dir = Path(__file__).parent
    output_dir = script_dir / args.output_dir
    output_dir.mkdir(exist_ok=True)
    
    # Step 1: Read trajectories
    print(f"Reading trajectories from: {args.csv}")
    trajectories = read_trajectories_from_csv(args.csv)
    print(f"Total trajectories: {len(trajectories)}")
    
    # Step 2: Compute ensemble MSD
    print(f"\nComputing ensemble-averaged MSD...")
    msd_result = calculate_ensemble_msd(trajectories, max_lag_fraction=args.max_lag_fraction)
    
    if msd_result.tau.size == 0:
        print("Error: No MSD data to fit")
        return
    
    print(f"  dt (global):               {msd_result.dt:.6f} s")
    print(f"  n_max (steps):             {msd_result.n_max}")
    print(f"  tau_max (seconds):         {msd_result.tau[-1]:.2f} s")
    print(f"  Total trajectories:        {msd_result.total_trajectories}")
    
    # Step 3: Anomalous diffusion fit
    print(f"\n{'='*60}")
    print(f"Anomalous Diffusion MSD Fitting: MSD(tau) = 4*D_alpha*tau^alpha")
    print(f"{'='*60}")
    print(f"  D_alpha initial guess: {args.D_alpha_initial:.2e} um^2/s^alpha")
    print(f"  D_alpha bounds:        [{args.D_alpha_bounds[0]:.2e}, {args.D_alpha_bounds[1]:.2e}] um^2/s^alpha")
    print(f"  alpha initial guess:   {args.alpha_initial:.4f}")
    print(f"  alpha bounds:          [{args.alpha_bounds[0]:.4f}, {args.alpha_bounds[1]:.4f}]")
    
    # Check if manual fraction is specified
    if args.manual_fraction is not None:
        print(f"\nUsing manual fraction: {args.manual_fraction:.0%}")
        try:
            D_alpha, D_alpha_err, alpha, alpha_err, pcov, tau_f, msd_f, msd_p, R2, rss = fit_msd_anomalous_at_fraction(
                tau=msd_result.tau,
                msd=msd_result.msd,
                n_max=msd_result.n_max,
                dt=msd_result.dt,
                fit_fraction=args.manual_fraction,
                D_alpha_initial=args.D_alpha_initial,
                D_alpha_bounds=tuple(args.D_alpha_bounds),
                alpha_initial=args.alpha_initial,
                alpha_bounds=tuple(args.alpha_bounds),
            )
            
            n_steps = int(args.manual_fraction * msd_result.n_max)
            interval_results = {args.manual_fraction: (R2, rss)}
            
            # Create result object
            anomalous_result = AnomalousFitResult(
                D_alpha=D_alpha,
                D_alpha_error=D_alpha_err,
                alpha=alpha,
                alpha_error=alpha_err,
                pcov=pcov,
                tau_fit=tau_f,
                msd_fit=msd_f,
                msd_predicted=msd_p,
                R_squared=R2,
                RSS=rss,
                optimal_fraction=args.manual_fraction,
                n_fit_steps=n_steps,
                interval_results=interval_results,
            )
            
        except (ValueError, RuntimeError) as e:
            print(f"\nError during manual fraction fit: {e}")
            return
    else:
        # Run automatic optimization
        try:
            anomalous_result = fit_msd_anomalous(
                tau=msd_result.tau,
                msd=msd_result.msd,
                n_max=msd_result.n_max,
                dt=msd_result.dt,
                D_alpha_initial=args.D_alpha_initial,
                D_alpha_bounds=tuple(args.D_alpha_bounds),
                alpha_initial=args.alpha_initial,
                alpha_bounds=tuple(args.alpha_bounds),
                interval_step=args.interval_step,
            )
        except RuntimeError as e:
            print(f"\nError during anomalous fitting: {e}")
            return
    
    # Print anomalous fit results
    print(f"\n{'='*60}")
    print(f"Anomalous Diffusion Fit Results:")
    print(f"{'='*60}")
    interval_label = "Manual interval" if args.manual_fraction is not None else "Optimal interval"
    print(f"  {interval_label}:                {anomalous_result.optimal_fraction:.0%} ({anomalous_result.n_fit_steps} steps)")
    print(f"  Diffusion Coeff (D_alpha):     ({anomalous_result.D_alpha:.6e} +/- {anomalous_result.D_alpha_error:.2e}) um^2/s^alpha")
    print(f"  Anomalous exponent (alpha):    ({anomalous_result.alpha:.4f} +/- {anomalous_result.alpha_error:.4f})")
    print(f"  R^2:                           {anomalous_result.R_squared:.6f}")
    print(f"  RSS:                           {anomalous_result.RSS:.6e}")
    print(f"  Fit tau range:                 [{anomalous_result.tau_fit.min():.2f}, {anomalous_result.tau_fit.max():.2f}] s")
    
    # Interpret alpha value
    if anomalous_result.alpha < 0.95:
        diffusion_type = "Subdiffusion (confined/hindered motion)"
    elif anomalous_result.alpha < 1.05:
        diffusion_type = "Normal Brownian diffusion"
    elif anomalous_result.alpha < 1.95:
        diffusion_type = "Superdiffusion (active transport)"
    else:
        diffusion_type = "Ballistic motion"
    print(f"  Interpretation:                {diffusion_type}")
    print(f"{'='*60}")
    
    # Step 4: Generate plot
    csv_name = Path(args.csv).stem
    plot_filename = f"{csv_name}_anomalous_fit_{int(anomalous_result.optimal_fraction*100)}pct.svg"
    plot_path = output_dir / plot_filename
    
    print(f"\nGenerating plot...")
    plot_anomalous_fit(
        tau_fit=anomalous_result.tau_fit,
        msd_fit=anomalous_result.msd_fit,
        msd_predicted=anomalous_result.msd_predicted,
        D_alpha=anomalous_result.D_alpha,
        D_alpha_error=anomalous_result.D_alpha_error,
        alpha=anomalous_result.alpha,
        alpha_error=anomalous_result.alpha_error,
        R_squared=anomalous_result.R_squared,
        output_path=plot_path,
        optimal_fraction=anomalous_result.optimal_fraction,
        n_fit_steps=anomalous_result.n_fit_steps,
    )
    
    # Step 5: Optional comparison with normal diffusion (α=1)
    if args.compare_normal:
        print(f"\n{'='*60}")
        print(f"Normal Diffusion Comparison (alpha fixed at 1.0)")
        print(f"{'='*60}")
        
        try:
            # Fit with alpha fixed at 1.0
            from msd_fitting import fit_msd_linear
            
            linear_result = fit_msd_linear(
                tau=msd_result.tau,
                msd=msd_result.msd,
                n_max=msd_result.n_max,
                dt=msd_result.dt,
                fit_fraction=anomalous_result.optimal_fraction,
                D_initial=args.D_alpha_initial,
                D_bounds=tuple(args.D_alpha_bounds),
            )
            
            print(f"\nNormal diffusion (alpha=1, MSD = 4D*tau):")
            print(f"  D = ({linear_result.D:.6e} +/- {linear_result.D_error:.2e}) um^2/s")
            print(f"  R^2 = {linear_result.R_squared:.6f}")
            
            print(f"\nAnomalous diffusion (alpha free, MSD = 4*D_alpha*tau^alpha):")
            print(f"  D_alpha = ({anomalous_result.D_alpha:.6e} +/- {anomalous_result.D_alpha_error:.2e}) um^2/s^alpha")
            print(f"  alpha = ({anomalous_result.alpha:.4f} +/- {anomalous_result.alpha_error:.4f})")
            print(f"  R^2 = {anomalous_result.R_squared:.6f}")
            
            # Compare R² values
            R2_improvement = anomalous_result.R_squared - linear_result.R_squared
            if R2_improvement > 0.01:
                print(f"\nAnomalous model provides better fit (Delta_R^2 = {R2_improvement:+.4f}")
            else:
                print(f"\nBoth models fit similarly well (Delta_R^2 = {R2_improvement:+.4f}")
            
            print(f"{'='*60}")
            
        except Exception as e:
            print(f"Could not perform normal diffusion comparison: {e}")
    
    print(f"\nAnomalous diffusion fitting analysis complete!")
    print(f"Results saved to: {output_dir.resolve()}")


if __name__ == "__main__":
    main()
