"""
Anomalous diffusion MSD fitting module with drift.

This module implements fitting to the anomalous diffusion model with drift velocity:
    
    MSD(τ) = 4D_α·τ^α + v²·τ²
    
where:
    - MSD is the Mean Squared Displacement
    - τ (tau) is the time lag
    - D_α is the generalized diffusion coefficient
    - α (alpha) is the anomalous diffusion exponent
    - v is the drift/active transport velocity
    
The model combines:
    - Anomalous diffusion component: 4D_α·τ^α (with exponent α)
    - Ballistic drift component: v²·τ² (quadratic in time)

The exponent α characterizes the type of diffusion:
    - α < 1: Subdiffusion (confined/hindered motion)
    - α = 1: Normal Brownian diffusion
    - α > 1: Superdiffusion (active transport)
    
The fitting is performed over multiple time intervals (10%-90% of data) to find
the optimal range based on R² and residual analysis. Velocity bounds are estimated
from trajectory path lengths.

Physical Context:
    This model is appropriate for systems exhibiting both anomalous diffusion and
    directed motion, such as particles in viscoelastic media under flow or active
    transport in crowded biological environments.

Dependencies: numpy, scipy, matplotlib (optional for velocity histogram)
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
        n_trajectories_used: Number of trajectories after filtering (>=30 points)
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
class DriftAnomalousFitResult:
    """Container for anomalous MSD fit results with drift.
    
    Attributes:
        D_alpha: Generalized diffusion coefficient [μm²/s^α]
        D_alpha_error: Standard error on D_α
        alpha: Anomalous diffusion exponent (dimensionless)
        alpha_error: Standard error on α
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
    
    D_alpha: float
    D_alpha_error: float
    alpha: float
    alpha_error: float
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
    
    Filters trajectories with n_points >= min_points, computes velocity for each
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
            f"No valid trajectories found with >={min_points} points. "
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


def anomalous_drift_msd_model(tau: np.ndarray, D_alpha: float, alpha: float, v: float) -> np.ndarray:
    """Anomalous diffusion with drift model: MSD(τ) = 4D_α·τ^α + v²·τ²
    
    Combines anomalous diffusion (with exponent α) and ballistic drift.
    
    Args:
        tau: Time lag values [s]
        D_alpha: Generalized diffusion coefficient [μm²/s^α]
        alpha: Anomalous diffusion exponent (dimensionless)
        v: Drift velocity [μm/s]
    
    Returns:
        MSD values predicted by the model [μm²]
    """
    return 4.0 * D_alpha * np.power(tau, alpha) + v**2 * tau**2


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


def fit_msd_anomalous_drift_at_fraction(
    tau: np.ndarray,
    msd: np.ndarray,
    n_max: int,
    dt: float,
    fit_fraction: float,
    D_alpha_initial: float,
    D_alpha_bounds: Tuple[float, float],
    alpha_initial: float,
    alpha_bounds: Tuple[float, float],
    v_initial: float,
    v_bounds: Tuple[float, float],
) -> Tuple[float, float, float, float, float, float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]:
    """Fit MSD to anomalous drift model at a specific fraction of data.
    
    Returns:
        Tuple of (D_alpha, D_alpha_error, alpha, alpha_error, v, v_error, pcov,
                  tau_fit, msd_fit, msd_predicted, R_squared, RSS)
    """
    # Calculate number of steps for this fraction
    n_fit_steps = max(2, int(fit_fraction * n_max))
    tau_max = n_fit_steps * dt
    
    # Select data subset
    mask = tau <= tau_max
    tau_fit = tau[mask]
    msd_fit = msd[mask]
    
    # Remove NaN/inf and ensure tau > 0 for power operation
    valid = np.isfinite(tau_fit) & np.isfinite(msd_fit) & (tau_fit > 0)
    tau_fit = tau_fit[valid]
    msd_fit = msd_fit[valid]
    
    if tau_fit.size < 4:  # Need at least 4 points for 3-parameter fit
        raise ValueError(f"Insufficient points ({tau_fit.size}) for anomalous drift fit at fraction {fit_fraction}")
    
    # Perform the fit
    try:
        popt, pcov = curve_fit(
            anomalous_drift_msd_model,
            tau_fit,
            msd_fit,
            p0=[D_alpha_initial, alpha_initial, v_initial],
            bounds=(
                [D_alpha_bounds[0], alpha_bounds[0], v_bounds[0]],
                [D_alpha_bounds[1], alpha_bounds[1], v_bounds[1]]
            ),
            method='trf',
            maxfev=10000,
        )
    except RuntimeError as e:
        raise RuntimeError(f"Fitting failed at fraction {fit_fraction}: {e}")
    
    D_alpha_opt = float(popt[0])
    alpha_opt = float(popt[1])
    v_opt = float(popt[2])
    
    # Calculate errors
    if pcov is not None and np.all(np.isfinite(pcov)):
        D_alpha_error = float(np.sqrt(pcov[0, 0]))
        alpha_error = float(np.sqrt(pcov[1, 1]))
        v_error = float(np.sqrt(pcov[2, 2]))
    else:
        D_alpha_error = float('nan')
        alpha_error = float('nan')
        v_error = float('nan')
    
    # Calculate metrics
    msd_predicted = anomalous_drift_msd_model(tau_fit, D_alpha_opt, alpha_opt, v_opt)
    R_squared = calculate_r_squared(msd_fit, msd_predicted)
    RSS = calculate_rss(msd_fit, msd_predicted)
    
    return D_alpha_opt, D_alpha_error, alpha_opt, alpha_error, v_opt, v_error, pcov, tau_fit, msd_fit, msd_predicted, R_squared, RSS


def fit_msd_anomalous_drift(
    tau: np.ndarray,
    msd: np.ndarray,
    n_max: int,
    dt: float,
    velocity_stats: VelocityStats,
    D_alpha_initial: float = 1e-2,
    D_alpha_bounds: Tuple[float, float] = (1e-6, 1e2),
    alpha_initial: float = 1.0,
    alpha_bounds: Tuple[float, float] = (0.01, 2.0),
    interval_step: float = 0.10,
) -> DriftAnomalousFitResult:
    """Fit MSD to anomalous drift model testing multiple intervals.
    
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
        D_alpha_initial: Initial guess for D_α [μm²/s^α]
        D_alpha_bounds: Tuple (lower, upper) bounds for D_α [μm²/s^α]
        alpha_initial: Initial guess for α (default: 1.0 = normal diffusion)
        alpha_bounds: Tuple (lower, upper) bounds for α (default: 0.01-2.0)
        interval_step: Step size for interval search (default: 0.10 = 10%)
    
    Returns:
        DriftAnomalousFitResult with optimal fit parameters and all interval results
    """
    # Generate fractions to test: 0.1, 0.2, ..., 0.9
    fractions = np.arange(interval_step, 1.0, interval_step)
    
    # Store results for each fraction
    results: Dict[float, Tuple[float, float, float, float, float, float, float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, float]] = {}
    interval_results: Dict[float, Tuple[float, float]] = {}
    
    print(f"\nTesting intervals from {fractions[0]:.0%} to {fractions[-1]:.0%} in {interval_step:.0%} steps...")
    
    for frac in fractions:
        try:
            D_alpha, D_alpha_err, alpha, alpha_err, v, v_err, pcov, tau_f, msd_f, msd_p, R2, rss = fit_msd_anomalous_drift_at_fraction(
                tau, msd, n_max, dt, frac,
                D_alpha_initial, D_alpha_bounds,
                alpha_initial, alpha_bounds,
                velocity_stats.v_initial, velocity_stats.v_bounds
            )
            results[frac] = (D_alpha, D_alpha_err, alpha, alpha_err, v, v_err, R2, pcov, tau_f, msd_f, msd_p, rss, int(frac * n_max))
            interval_results[frac] = (R2, rss)
            print(f"  {frac:4.0%}: R^2 = {R2:.6f}, RSS = {rss:.6e}, D_alpha = {D_alpha:.6e}, alpha = {alpha:.4f}, v = {v:.6e}")
        except (ValueError, RuntimeError) as e:
            print(f"  {frac:4.0%}: Failed - {e}")
            continue
    
    if len(results) == 0:
        raise RuntimeError("No valid fits obtained across any interval")
    
    # Find optimal fraction
    # Step 1: Find maximum R²
    max_R2 = max(r[6] for r in results.values())
    
    # Step 2: Find all fractions within 0.01 of max R²
    candidates = [(frac, r) for frac, r in results.items() if r[6] >= max_R2 - 0.01]
    
    # Step 3: Among candidates, select minimum RSS
    optimal_fraction = min(candidates, key=lambda x: x[1][11])[0]
    
    # Extract optimal result
    D_alpha, D_alpha_err, alpha, alpha_err, v, v_err, R2, pcov, tau_f, msd_f, msd_p, rss, n_steps = results[optimal_fraction]
    
    print(f"\nOptimal interval: {optimal_fraction:.0%} (R^2 = {R2:.6f}, RSS = {rss:.6e})")
    
    return DriftAnomalousFitResult(
        D_alpha=D_alpha,
        D_alpha_error=D_alpha_err,
        alpha=alpha,
        alpha_error=alpha_err,
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


def plot_anomalous_drift_fit(
    tau_fit: np.ndarray,
    msd_fit: np.ndarray,
    msd_predicted: np.ndarray,
    D_alpha: float,
    D_alpha_error: float,
    alpha: float,
    alpha_error: float,
    v: float,
    v_error: float,
    R_squared: float,
    output_path: Path,
    optimal_fraction: float = 0.10,
    n_fit_steps: int = 0,
) -> None:
    """Create and save a plot of MSD data with anomalous drift fit.
    
    Plots data points and fit line for anomalous drift model: MSD = 4D_α·τ^α + v²·τ²
    
    Args:
        tau_fit: Time lag values used in the fit
        msd_fit: MSD data values used in the fit
        msd_predicted: Predicted MSD values from the fitted model
        D_alpha: Fitted generalized diffusion coefficient [μm²/s^α]
        D_alpha_error: Standard error on D_α
        alpha: Fitted anomalous exponent
        alpha_error: Standard error on α
        v: Fitted drift velocity [μm/s]
        v_error: Standard error on v
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
            label=r'Anomalous Drift Fit: MSD = 4$D_\alpha \tau^\alpha$ + $v^2\tau^2$', zorder=3)
    
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
        r'$v = (%.2e \pm %.1e)\ \mu m/s$' % (v, v_error),
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
    """Main function for anomalous diffusion with drift MSD fitting analysis."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Fit MSD data to anomalous diffusion with drift model: MSD = 4D_α·τ^α + v²·τ²',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    # Automatic optimization (tests 10%%-90%%)
    python msd_fitting_anomalous_drift.py Data/Experiment8_spots_40minstep.csv
    
    # Manual interval selection (fit at 30%%)
    python msd_fitting_anomalous_drift.py Data/Experiment8_spots_40minstep.csv --manual-fraction 0.3
    
    # Custom alpha bounds for subdiffusion
    python msd_fitting_anomalous_drift.py Data/Experiment8_spots_40minstep.csv --alpha-bounds 0.1 0.99
    
    # Generate velocity histogram
    python msd_fitting_anomalous_drift.py Data/Experiment8_spots_40minstep.csv --plot-velocity
    
    # Compare with non-drifting anomalous model
    python msd_fitting_anomalous_drift.py Data/Experiment8_spots_40minstep.csv --compare-no-drift
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
        '--plot-velocity',
        action='store_true',
        help='Generate velocity distribution histogram'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        default='anomalous_drift_fits',
        help='Output directory for results (default: anomalous_drift_fits/)'
    )
    
    parser.add_argument(
        '--compare-no-drift',
        action='store_true',
        help='Also run fit without drift term (MSD = 4D_α·τ^α) for comparison'
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
    print(f"  tau_max (seconds):         {msd_result.tau[-1]:.2f} s")
    print(f"  Total trajectories:        {msd_result.total_trajectories}")
    
    # Step 4: Anomalous diffusion with drift fit
    print(f"\n{'='*60}")
    print(f"Anomalous Diffusion with Drift Fitting: MSD(tau) = 4*D_alpha*tau^alpha + v^2*tau^2")
    print(f"{'='*60}")
    print(f"  D_alpha initial guess: {args.D_alpha_initial:.2e} um^2/s^alpha")
    print(f"  D_alpha bounds:        [{args.D_alpha_bounds[0]:.2e}, {args.D_alpha_bounds[1]:.2e}] um^2/s^alpha")
    print(f"  alpha initial guess:   {args.alpha_initial:.4f}")
    print(f"  alpha bounds:          [{args.alpha_bounds[0]:.4f}, {args.alpha_bounds[1]:.4f}]")
    print(f"  v initial guess:   {velocity_stats.v_initial:.6f} um/s")
    print(f"  v bounds:          [{velocity_stats.v_bounds[0]:.6f}, {velocity_stats.v_bounds[1]:.6f}] um/s")
    
    # Check if manual fraction is specified
    if args.manual_fraction is not None:
        print(f"\nUsing manual fraction: {args.manual_fraction:.0%}")
        try:
            D_alpha, D_alpha_err, alpha, alpha_err, v, v_err, pcov, tau_f, msd_f, msd_p, R2, rss = fit_msd_anomalous_drift_at_fraction(
                tau=msd_result.tau,
                msd=msd_result.msd,
                n_max=msd_result.n_max,
                dt=msd_result.dt,
                fit_fraction=args.manual_fraction,
                D_alpha_initial=args.D_alpha_initial,
                D_alpha_bounds=tuple(args.D_alpha_bounds),
                alpha_initial=args.alpha_initial,
                alpha_bounds=tuple(args.alpha_bounds),
                v_initial=velocity_stats.v_initial,
                v_bounds=velocity_stats.v_bounds,
            )
            
            n_steps = int(args.manual_fraction * msd_result.n_max)
            interval_results = {args.manual_fraction: (R2, rss)}
            
            # Create result object
            drift_anomalous_result = DriftAnomalousFitResult(
                D_alpha=D_alpha,
                D_alpha_error=D_alpha_err,
                alpha=alpha,
                alpha_error=alpha_err,
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
            drift_anomalous_result = fit_msd_anomalous_drift(
                tau=msd_result.tau,
                msd=msd_result.msd,
                n_max=msd_result.n_max,
                dt=msd_result.dt,
                velocity_stats=velocity_stats,
                D_alpha_initial=args.D_alpha_initial,
                D_alpha_bounds=tuple(args.D_alpha_bounds),
                alpha_initial=args.alpha_initial,
                alpha_bounds=tuple(args.alpha_bounds),
                interval_step=args.interval_step,
            )
        except RuntimeError as e:
            print(f"\nError during anomalous drift fitting: {e}")
            return
    
    # Print drift anomalous fit results
    print(f"\n{'='*60}")
    print(f"Anomalous Diffusion with Drift Fit Results:")
    print(f"{'='*60}")
    interval_label = "Manual interval" if args.manual_fraction is not None else "Optimal interval"
    print(f"  {interval_label}:                {drift_anomalous_result.optimal_fraction:.0%} ({drift_anomalous_result.n_fit_steps} steps)")
    print(f"  Diffusion Coeff (D_alpha):     ({drift_anomalous_result.D_alpha:.6e} +/- {drift_anomalous_result.D_alpha_error:.2e}) um^2/s^alpha")
    print(f"  Anomalous exponent (alpha):    ({drift_anomalous_result.alpha:.4f} +/- {drift_anomalous_result.alpha_error:.4f})")
    print(f"  Drift Velocity (v):            ({drift_anomalous_result.v:.6e} +/- {drift_anomalous_result.v_error:.2e}) um/s")
    print(f"  R^2:                           {drift_anomalous_result.R_squared:.6f}")
    print(f"  RSS:                           {drift_anomalous_result.RSS:.6e}")
    print(f"  Fit tau range:                 [{drift_anomalous_result.tau_fit.min():.2f}, {drift_anomalous_result.tau_fit.max():.2f}] s")
    
    # Interpret alpha value
    if drift_anomalous_result.alpha < 0.95:
        diffusion_type = "Subdiffusion (confined/hindered motion)"
    elif drift_anomalous_result.alpha < 1.05:
        diffusion_type = "Normal Brownian diffusion"
    elif drift_anomalous_result.alpha < 1.95:
        diffusion_type = "Superdiffusion (active transport)"
    else:
        diffusion_type = "Ballistic motion"
    print(f"  Interpretation:                {diffusion_type} with drift")
    print(f"{'='*60}")
    
    # Step 5: Generate plot
    csv_name = Path(args.csv).stem
    plot_filename = f"{csv_name}_anomalous_drift_fit_{int(drift_anomalous_result.optimal_fraction*100)}pct.svg"
    plot_path = output_dir / plot_filename
    
    print(f"\nGenerating plot...")
    plot_anomalous_drift_fit(
        tau_fit=drift_anomalous_result.tau_fit,
        msd_fit=drift_anomalous_result.msd_fit,
        msd_predicted=drift_anomalous_result.msd_predicted,
        D_alpha=drift_anomalous_result.D_alpha,
        D_alpha_error=drift_anomalous_result.D_alpha_error,
        alpha=drift_anomalous_result.alpha,
        alpha_error=drift_anomalous_result.alpha_error,
        v=drift_anomalous_result.v,
        v_error=drift_anomalous_result.v_error,
        R_squared=drift_anomalous_result.R_squared,
        output_path=plot_path,
        optimal_fraction=drift_anomalous_result.optimal_fraction,
        n_fit_steps=drift_anomalous_result.n_fit_steps,
    )
    
    # Step 6: Optional comparison with non-drifting anomalous model
    if args.compare_no_drift:
        print(f"\n{'='*60}")
        print(f"Non-Drifting Anomalous Model Comparison")
        print(f"{'='*60}")
        
        try:
            from msd_fitting_anomalous import fit_msd_anomalous
            
            no_drift_result = fit_msd_anomalous(
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
            
            print(f"\nNon-drifting model (MSD = 4*D_alpha*tau^alpha):")
            print(f"  D_alpha = ({no_drift_result.D_alpha:.6e} +/- {no_drift_result.D_alpha_error:.2e}) um^2/s^alpha")
            print(f"  alpha = ({no_drift_result.alpha:.4f} +/- {no_drift_result.alpha_error:.4f})")
            print(f"  R^2 = {no_drift_result.R_squared:.6f}")
            
            print(f"\nDrifting model (MSD = 4*D_alpha*tau^alpha + v^2*tau^2):")
            print(f"  D_alpha = ({drift_anomalous_result.D_alpha:.6e} +/- {drift_anomalous_result.D_alpha_error:.2e}) um^2/s^alpha")
            print(f"  alpha = ({drift_anomalous_result.alpha:.4f} +/- {drift_anomalous_result.alpha_error:.4f})")
            print(f"  v = ({drift_anomalous_result.v:.6e} +/- {drift_anomalous_result.v_error:.2e}) um/s")
            print(f"  R^2 = {drift_anomalous_result.R_squared:.6f}")
            
            # Compare R² values
            R2_improvement = drift_anomalous_result.R_squared - no_drift_result.R_squared
            if R2_improvement > 0.01:
                print(f"\nDrifting model provides better fit (Delta_R^2 = {R2_improvement:+.4f}")
            else:
                print(f"\nBoth models fit similarly well (Delta_R^2 = {R2_improvement:+.4f}")
            
            print(f"{'='*60}")
            
        except Exception as e:
            print(f"Could not perform non-drifting comparison: {e}")
    
    print(f"\nAnomalous diffusion with drift fitting analysis complete!")
    print(f"Results saved to: {output_dir.resolve()}")


if __name__ == "__main__":
    main()
