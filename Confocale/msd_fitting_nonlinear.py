"""
Nonlinear MSD fitting module for extracting Diffusion Coefficient (D) and Velocity (v).

This module implements fitting to the nonlinear model:
    
    MSD(τ) = 4D·τ + v²·τ²
    
where:
    - MSD is the Mean Squared Displacement
    - τ (tau) is the time lag
    - D is the Diffusion Coefficient
    - v is the drift/active transport velocity

The model captures both diffusive (4D·τ) and ballistic (v²·τ²) motion components.
The fitting is performed over multiple time intervals (10%-90% of data) to find
the optimal range based on R² and residual analysis.

Physical Context:
    For particles in 84% Glycerol/Water solution, we expect:
    - D: low values, typically 10^-6 to 10^-1 μm²/s
    - v: depends on presence of active transport/flow, typically 0-1 μm/s

Dependencies: numpy, scipy, matplotlib (optional for velocity histogram)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Tuple, Union

import numpy as np
from scipy.optimize import curve_fit

from data_reader import Trajectory, read_trajectories_from_csv, estimate_global_time_step
from msd_analyzer import calculate_ensemble_msd
from msd_fitting import fit_msd_linear, calculate_r_squared
from plot_msd_fit import plot_msd_with_nonlinear_fit


@dataclass(frozen=True)
class VelocityStats:
    """Container for velocity analysis results.
    
    Attributes:
        velocities: Array of computed velocities for each trajectory [μm/s]
        mean: Mean velocity [μm/s]
        median: Median velocity [μm/s]
        std: Standard deviation [μm/s]
        v_initial: Initial guess for fitting (median) [μm/s]
        v_bounds: Tuple (lower, upper) bounds for v [μm/s]
        n_trajectories_total: Total number of trajectories before filtering
        n_trajectories_used: Number of trajectories after filtering (≥30 points)
    """
    
    velocities: np.ndarray
    mean: float
    median: float
    std: float
    v_initial: float
    v_bounds: Tuple[float, float]
    n_trajectories_total: int
    n_trajectories_used: int


@dataclass(frozen=True)
class NonlinearFitResult:
    """Container for nonlinear MSD fit results.
    
    Attributes:
        D: Diffusion coefficient [μm²/s]
        D_error: Standard error on D
        v: Drift velocity [μm/s]
        v_error: Standard error on v
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
    
    D: float
    D_error: float
    v: float
    v_error: float
    pcov: np.ndarray
    tau_fit: np.ndarray
    msd_fit: np.ndarray
    msd_predicted: np.ndarray
    R_squared: float
    RSS: float
    optimal_fraction: float
    n_fit_steps: int
    interval_results: Dict[float, Tuple[float, float]]


def compute_trajectory_velocity(traj: Trajectory) -> float:
    """Compute velocity as total path length / total time duration.
    
    Path length L = Σ sqrt((x[i+1]-x[i])² + (y[i+1]-y[i])²)
    Time duration T = time[-1] - time[0]
    Velocity v = L / T
    
    Args:
        traj: Single trajectory with x, y, time arrays
    
    Returns:
        Velocity in [μm/s] (or position units / time units)
    """
    if traj.n_points < 2:
        return float('nan')
    
    # Compute path length
    dx = np.diff(traj.x)
    dy = np.diff(traj.y)
    segment_lengths = np.sqrt(dx**2 + dy**2)
    total_path_length = float(np.sum(segment_lengths))
    
    # Compute time duration
    time_duration = float(traj.time[-1] - traj.time[0])
    
    if time_duration <= 0:
        return float('nan')
    
    return total_path_length / time_duration


def analyze_velocities(
    trajectories: Mapping[Union[int, str], Trajectory],
    min_points: int = 30,
) -> VelocityStats:
    """Analyze velocities from trajectories to establish fitting parameters.
    
    Filters trajectories with n_points ≥ min_points, computes velocity for each
    as v = (total path length) / (total time), and calculates statistics.
    
    Args:
        trajectories: Mapping of track ID to Trajectory
        min_points: Minimum number of points required (default: 30)
    
    Returns:
        VelocityStats containing statistics and recommended bounds/initial guess
    """
    n_total = len(trajectories)
    
    # Filter trajectories and compute velocities
    velocities_list = []
    for traj in trajectories.values():
        if traj.n_points >= min_points:
            v = compute_trajectory_velocity(traj)
            if np.isfinite(v) and v > 0:
                velocities_list.append(v)
    
    if len(velocities_list) == 0:
        raise ValueError(
            f"No valid trajectories found with ≥{min_points} points. "
            f"Total trajectories: {n_total}"
        )
    
    velocities = np.array(velocities_list)
    n_used = len(velocities)
    
    # Compute statistics
    mean_v = float(np.mean(velocities))
    median_v = float(np.median(velocities))
    std_v = float(np.std(velocities))
    
    # Set bounds: (0, mean + 3*std) to cover most of distribution
    v_lower = 0.0
    v_upper = mean_v + 3.0 * std_v
    
    # Initial guess: median (robust to outliers)
    v_initial = median_v
    
    return VelocityStats(
        velocities=velocities,
        mean=mean_v,
        median=median_v,
        std=std_v,
        v_initial=v_initial,
        v_bounds=(v_lower, v_upper),
        n_trajectories_total=n_total,
        n_trajectories_used=n_used,
    )


def nonlinear_msd_model(tau: np.ndarray, D: float, v: float) -> np.ndarray:
    """Nonlinear MSD model: MSD(τ) = 4D·τ + v²·τ²
    
    Captures both diffusive (linear in τ) and ballistic (quadratic in τ) behavior.
    
    Args:
        tau: Time lag values [s]
        D: Diffusion coefficient [μm²/s]
        v: Drift velocity [μm/s]
    
    Returns:
        MSD values predicted by the model [μm²]
    """
    return 4.0 * D * tau + v**2 * tau**2


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


def fit_msd_nonlinear_at_fraction(
    tau: np.ndarray,
    msd: np.ndarray,
    n_max: int,
    dt: float,
    fit_fraction: float,
    D_initial: float,
    D_bounds: Tuple[float, float],
    v_initial: float,
    v_bounds: Tuple[float, float],
) -> Tuple[float, float, float, float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]:
    """Fit MSD to nonlinear model at a specific fraction of data.
    
    Returns:
        Tuple of (D, D_error, v, v_error, pcov, tau_fit, msd_fit, msd_predicted, R_squared, RSS)
    """
    # Calculate number of steps for this fraction
    n_fit_steps = max(2, int(fit_fraction * n_max))
    tau_max = n_fit_steps * dt
    
    # Select data subset
    mask = tau <= tau_max
    tau_fit = tau[mask]
    msd_fit = msd[mask]
    
    # Remove NaN/inf
    valid = np.isfinite(tau_fit) & np.isfinite(msd_fit)
    tau_fit = tau_fit[valid]
    msd_fit = msd_fit[valid]
    
    if tau_fit.size < 3:  # Need at least 3 points for 2-parameter fit
        raise ValueError(f"Insufficient points ({tau_fit.size}) for nonlinear fit at fraction {fit_fraction}")
    
    # Perform the fit
    try:
        popt, pcov = curve_fit(
            nonlinear_msd_model,
            tau_fit,
            msd_fit,
            p0=[D_initial, v_initial],
            bounds=([D_bounds[0], v_bounds[0]], [D_bounds[1], v_bounds[1]]),
            method='trf',
        )
    except RuntimeError as e:
        raise RuntimeError(f"Fitting failed at fraction {fit_fraction}: {e}")
    
    D_opt = float(popt[0])
    v_opt = float(popt[1])
    
    # Calculate errors
    if pcov is not None and np.all(np.isfinite(pcov)):
        D_error = float(np.sqrt(pcov[0, 0]))
        v_error = float(np.sqrt(pcov[1, 1]))
    else:
        D_error = float('nan')
        v_error = float('nan')
    
    # Calculate metrics
    msd_predicted = nonlinear_msd_model(tau_fit, D_opt, v_opt)
    R_squared = calculate_r_squared(msd_fit, msd_predicted)
    RSS = calculate_rss(msd_fit, msd_predicted)
    
    return D_opt, D_error, v_opt, v_error, pcov, tau_fit, msd_fit, msd_predicted, R_squared, RSS


def fit_msd_nonlinear(
    tau: np.ndarray,
    msd: np.ndarray,
    n_max: int,
    dt: float,
    velocity_stats: VelocityStats,
    D_initial: float = 1e-2,
    D_bounds: Tuple[float, float] = (9e-4, 1.5e-1),
    interval_step: float = 0.10,
) -> NonlinearFitResult:
    """Fit MSD to nonlinear model testing multiple intervals.
    
    Tests intervals from 10% to 90% in steps of interval_step (default 10%).
    Selects optimal interval based on:
    1. Maximum R²
    2. If multiple intervals have R² within 0.01, select minimum RSS
    
    Args:
        tau: Array of time lag values [s]
        msd: Array of MSD values [μm²]
        n_max: Maximum number of lag steps available
        dt: Time step [s]
        velocity_stats: VelocityStats from velocity analysis
        D_initial: Initial guess for D [μm²/s]
        D_bounds: Tuple (lower, upper) bounds for D [μm²/s]
        interval_step: Step size for interval search (default: 0.10 = 10%)
    
    Returns:
        NonlinearFitResult with optimal fit parameters and all interval results
    """
    # Generate fractions to test: 0.1, 0.2, ..., 0.9
    fractions = np.arange(interval_step, 1.0, interval_step)
    
    # Store results for each fraction
    results: Dict[float, Tuple[float, float, float, float, float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]] = {}
    interval_results: Dict[float, Tuple[float, float]] = {}
    
    print(f"\nTesting intervals from {fractions[0]:.0%} to {fractions[-1]:.0%} in {interval_step:.0%} steps...")
    
    for frac in fractions:
        try:
            D, D_err, v, v_err, pcov, tau_f, msd_f, msd_p, R2, rss = fit_msd_nonlinear_at_fraction(
                tau, msd, n_max, dt, frac,
                D_initial, D_bounds,
                velocity_stats.v_initial, velocity_stats.v_bounds
            )
            results[frac] = (D, D_err, v, v_err, R2, pcov, tau_f, msd_f, msd_p, rss, int(frac * n_max))
            interval_results[frac] = (R2, rss)
            print(f"  {frac:4.0%}: R^2 = {R2:.6f}, RSS = {rss:.6e}, D = {D:.6e}, v = {v:.6e}")
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
    D, D_err, v, v_err, R2, pcov, tau_f, msd_f, msd_p, rss, n_steps = results[optimal_fraction]
    
    print(f"\nOptimal interval: {optimal_fraction:.0%} (R^2 = {R2:.6f}, RSS = {rss:.6e})")
    
    return NonlinearFitResult(
        D=D,
        D_error=D_err,
        v=v,
        v_error=v_err,
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


def plot_velocity_histogram(
    velocity_stats: VelocityStats,
    output_path: Path,
) -> None:
    """Plot histogram of velocity distribution with statistics.
    
    Args:
        velocity_stats: VelocityStats containing velocity data
        output_path: Path to save the figure
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Warning: matplotlib not available, skipping velocity histogram")
        return
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Histogram
    ax.hist(velocity_stats.velocities, bins=30, color='C0', alpha=0.7, edgecolor='black')
    
    # Mark mean and median
    ax.axvline(velocity_stats.mean, color='C3', linestyle='--', linewidth=2, label=f'Mean: {velocity_stats.mean:.4f} μm/s')
    ax.axvline(velocity_stats.median, color='C2', linestyle='--', linewidth=2, label=f'Median: {velocity_stats.median:.4f} μm/s')
    
    # Labels
    ax.set_xlabel(r'Velocity [μm/s]', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.grid(True, linestyle=':', alpha=0.4)
    ax.legend(loc='upper right', fontsize=10)
    
    # Minimal text annotation
    textstr = f'n = {velocity_stats.n_trajectories_used}'
    props = dict(boxstyle='round', facecolor='white', alpha=0.9)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Velocity histogram saved to: {output_path.resolve()}")


def main() -> None:
    """Main function for nonlinear MSD fitting analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Fit MSD data to nonlinear model: MSD = 4D·τ + v²·τ²',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    # Automatic optimization (tests 10%%-90%%)
    python msd_fitting_nonlinear.py Data/Experiment8_spots_40minstep.csv
    
    # Manual interval selection (fit at 30%%)
    python msd_fitting_nonlinear.py Data/Experiment8_spots_40minstep.csv --manual-fraction 0.3
    
    # Compare with linear model
    python msd_fitting_nonlinear.py Data/Experiment8_spots_40minstep.csv --compare-linear
    
    # Generate velocity histogram
    python msd_fitting_nonlinear.py Data/Experiment8_spots_40minstep.csv --plot-velocity
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
        '--min-points',
        type=int,
        default=30,
        help='Minimum number of points required per trajectory for velocity analysis (default: 30)'
    )
    
    parser.add_argument(
        '--D-initial',
        type=float,
        default=1e-2,
        help='Initial guess for D [μm²/s] (default: 1e-2)'
    )
    
    parser.add_argument(
        '--D-bounds',
        type=float,
        nargs=2,
        default=[9e-4, 1.5e-1],
        metavar=('LOWER', 'UPPER'),
        help='Bounds for D [μm²/s] (default: 9e-4 1.5e-1)'
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
        '--plot-velocity',
        action='store_true',
        help='Generate velocity distribution histogram'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        default='nonlinear_fits',
        help='Output directory for plots (default: nonlinear_fits/)'
    )
    
    parser.add_argument(
        '--compare-linear',
        action='store_true',
        help='Also run linear fit and compare D values'
    )
    
    args = parser.parse_args()
    
    # Validate manual-fraction if provided
    if args.manual_fraction is not None:
        if not (0.1 <= args.manual_fraction <= 0.9):
            print(f"Error: --manual-fraction must be between 0.1 and 0.9 (got {args.manual_fraction})")
            return
    
    # Create output directory relative to script location
    script_dir = Path(__file__).parent
    output_dir = script_dir / args.output_dir
    output_dir.mkdir(exist_ok=True)
    
    # Step 1: Read trajectories
    print(f"Reading trajectories from: {args.csv}")
    trajectories = read_trajectories_from_csv(args.csv)
    print(f"Total trajectories: {len(trajectories)}")
    
    # Step 2: Velocity analysis
    print(f"\n{'='*60}")
    print(f"Velocity Analysis (filtering trajectories with >={args.min_points} points)")
    print(f"{'='*60}")
    
    velocity_stats = analyze_velocities(trajectories, min_points=args.min_points)
    
    print(f"  Trajectories used: {velocity_stats.n_trajectories_used}/{velocity_stats.n_trajectories_total}")
    print(f"  Mean velocity:     {velocity_stats.mean:.6f} um/s")
    print(f"  Median velocity:   {velocity_stats.median:.6f} um/s")
    print(f"  Std deviation:     {velocity_stats.std:.6f} um/s")
    print(f"  Initial guess (v): {velocity_stats.v_initial:.6f} um/s")
    print(f"  Bounds for v:      [{velocity_stats.v_bounds[0]:.6f}, {velocity_stats.v_bounds[1]:.6f}] um/s")
    print(f"{'='*60}")
    
    # Optional: Plot velocity histogram
    if args.plot_velocity:
        csv_name = Path(args.csv).stem
        velocity_plot_path = output_dir / f"{csv_name}_velocity_distribution.svg"
        plot_velocity_histogram(velocity_stats, velocity_plot_path)
    
    # Step 3: Compute ensemble MSD
    print(f"\nComputing ensemble-averaged MSD...")
    msd_result = calculate_ensemble_msd(trajectories, max_lag_fraction=args.max_lag_fraction)
    
    if msd_result.tau.size == 0:
        print("Error: No MSD data to fit")
        return
    
    print(f"  dt (global):               {msd_result.dt:.6f} s")
    print(f"  n_max (steps):             {msd_result.n_max}")
    print(f"  tau_max (seconds):           {msd_result.tau[-1]:.2f} s")
    print(f"  Total trajectories:        {msd_result.total_trajectories}")
    
    # Step 4: Nonlinear fit
    print(f"\n{'='*60}")
    print(f"Nonlinear MSD Fitting: MSD(tau) = 4D*tau + v^2*tau^2")
    print(f"{'='*60}")
    print(f"  D initial guess: {args.D_initial:.2e} um^2/s")
    print(f"  D bounds:        [{args.D_bounds[0]:.2e}, {args.D_bounds[1]:.2e}] um^2/s")
    print(f"  v initial guess: {velocity_stats.v_initial:.6f} um/s")
    print(f"  v bounds:        [{velocity_stats.v_bounds[0]:.6f}, {velocity_stats.v_bounds[1]:.6f}] um/s")
    
    # Check if manual fraction is specified
    if args.manual_fraction is not None:
        print(f"\nUsing manual fraction: {args.manual_fraction:.0%}")
        try:
            D, D_err, v, v_err, pcov, tau_f, msd_f, msd_p, R2, rss = fit_msd_nonlinear_at_fraction(
                tau=msd_result.tau,
                msd=msd_result.msd,
                n_max=msd_result.n_max,
                dt=msd_result.dt,
                fit_fraction=args.manual_fraction,
                D_initial=args.D_initial,
                D_bounds=tuple(args.D_bounds),
                v_initial=velocity_stats.v_initial,
                v_bounds=velocity_stats.v_bounds,
            )
            
            n_steps = int(args.manual_fraction * msd_result.n_max)
            interval_results = {args.manual_fraction: (R2, rss)}
            
            # Create result object
            nonlinear_result = NonlinearFitResult(
                D=D,
                D_error=D_err,
                v=v,
                v_error=v_err,
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
            nonlinear_result = fit_msd_nonlinear(
                tau=msd_result.tau,
                msd=msd_result.msd,
                n_max=msd_result.n_max,
                dt=msd_result.dt,
                velocity_stats=velocity_stats,
                D_initial=args.D_initial,
                D_bounds=tuple(args.D_bounds),
                interval_step=args.interval_step,
            )
        except RuntimeError as e:
            print(f"\nError during nonlinear fitting: {e}")
            return
    
    # Print nonlinear fit results
    print(f"\n{'='*60}")
    print(f"Nonlinear Fit Results:")
    print(f"{'='*60}")
    interval_label = "Manual interval" if args.manual_fraction is not None else "Optimal interval"
    print(f"  {interval_label}:            {nonlinear_result.optimal_fraction:.0%} ({nonlinear_result.n_fit_steps} steps)")
    print(f"  Diffusion Coeff (D):       ({nonlinear_result.D:.6e} +/- {nonlinear_result.D_error:.2e}) um^2/s")
    print(f"  Drift Velocity (v):        ({nonlinear_result.v:.6e} +/- {nonlinear_result.v_error:.2e}) um/s")
    print(f"  R^2:                        {nonlinear_result.R_squared:.6f}")
    print(f"  RSS:                       {nonlinear_result.RSS:.6e}")
    print(f"  Fit tau range:               [{nonlinear_result.tau_fit.min():.2f}, {nonlinear_result.tau_fit.max():.2f}] s")
    print(f"{'='*60}")
    
    # Step 5: Optional comparison with linear fit
    if args.compare_linear:
        print(f"\n{'='*60}")
        print(f"Linear Fit Comparison (first 10% of data)")
        print(f"{'='*60}")
        
        try:
            linear_result = fit_msd_linear(
                tau=msd_result.tau,
                msd=msd_result.msd,
                n_max=msd_result.n_max,
                dt=msd_result.dt,
                fit_fraction=0.10,
                D_initial=args.D_initial,
                D_bounds=tuple(args.D_bounds),
            )
            
            print(f"\nD with linear model (MSD = 4D*tau) is:        ({linear_result.D:.6e} +/- {linear_result.D_error:.2e}) um^2/s")
            print(f"D with non-linear model (MSD = 4D*tau + v^2*tau^2) is: ({nonlinear_result.D:.6e} +/- {nonlinear_result.D_error:.2e}) um^2/s")
            print(f"\nLinear fit R^2:     {linear_result.R_squared:.6f}")
            print(f"Nonlinear fit R^2:  {nonlinear_result.R_squared:.6f}")
            print(f"{'='*60}")
            
        except Exception as e:
            print(f"Could not perform linear fit comparison: {e}")
    
    # Plot the nonlinear fit
    plot_filename = output_dir / f"{Path(args.csv).stem}_nonlinear_fit.svg"
    print(f"\nSaving nonlinear fit plot to: {plot_filename}")
    plot_msd_with_nonlinear_fit(
        tau_fit=nonlinear_result.tau_fit,
        msd_fit=nonlinear_result.msd_fit,
        msd_fit_line=nonlinear_result.msd_predicted,
        D=nonlinear_result.D,
        D_error=nonlinear_result.D_error,
        v=nonlinear_result.v,
        v_error=nonlinear_result.v_error,
        R_squared=nonlinear_result.R_squared,
        output_path=plot_filename,
        optimal_fraction=nonlinear_result.optimal_fraction,
        n_fit_steps=nonlinear_result.n_fit_steps
    )

    print(f"\nNonlinear fitting analysis complete!")
    print(f"Results saved to: {output_dir.resolve()}")


if __name__ == "__main__":
    main()
