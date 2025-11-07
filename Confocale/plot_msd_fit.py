"""
Plot MSD with linear fit and extract Diffusion Coefficient.

This script performs the following steps:
1. Reads trajectory data from CSV and computes ensemble-averaged MSD
2. Fits the MSD data (for τ ≤ 15s) to the linear model: MSD(τ) = 2D·τ
3. Generates a publication-quality plot showing:
   - Original MSD data points
   - Fitted linear model (over the fitting range)
   - Text box with D and R² values
4. Saves the plot to a file

Physical Context:
    This analysis is designed for particles diffusing in highly viscous media
    (e.g., 80% Glycerol/Water solution), where the Diffusion Coefficient is
    expected to be in the range 10^-6 to 10^-2 μm²/s.

Dependencies: matplotlib, numpy, scipy (via msd_fitting)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from data_reader import read_trajectories_from_csv
from msd_analyzer import calculate_ensemble_msd
from msd_fitting import fit_msd_linear


def plot_msd_with_fit(
    tau_fit: np.ndarray,
    msd_fit: np.ndarray,
    msd_fit_line: np.ndarray,
    D: float,
    D_error: float,
    R_squared: float,
    output_path: Path,
    fit_fraction: float = 0.10,
    n_fit_steps: int = 0,
) -> None:
    """Create and save a plot of MSD data with linear fit.
    
    Only plots the data points and fit line in the relevant range (first fit_fraction of data).
    
    Args:
        tau_fit: Time lag values used in the fit
        msd_fit: MSD data values used in the fit
        msd_fit_line: Fitted MSD values (model predictions)
        D: Fitted Diffusion Coefficient
        D_error: Standard error on D
        R_squared: Coefficient of determination
        output_path: Path to save the figure
        fit_fraction: Fraction of data used for fitting
        n_fit_steps: Number of lag steps used for fitting
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot ONLY the data points used in the fit
    ax.plot(tau_fit, msd_fit, 'o', color='C0', markersize=8, 
            alpha=0.7, label='MSD Data', zorder=2)
    
    # Plot the fitted line (over the same range)
    ax.plot(tau_fit, msd_fit_line, '-', color='C3', linewidth=2.5,
            label=f'Linear Fit: MSD = 2Dτ', zorder=3)
    
    # Labels (no title for publication standards)
    ax.set_xlabel(r'Time Lag $\tau$ [s]', fontsize=12)
    ax.set_ylabel(r'MSD [μm$^2$]', fontsize=12)
    
    # Grid
    ax.grid(True, linestyle=':', alpha=0.4)
    
    # Legend
    ax.legend(loc='upper left', fontsize=10, framealpha=0.9)
    
    # Text box with fit results - minimal, white background
    textstr = '\n'.join([
        r'$D = (%.2e \pm %.1e)\ \mu m^2/s$' % (D, D_error),
        r'$R^2 = %.4f$' % R_squared,
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
    
    print(f"\nPlot saved to: {output_path.resolve()}")


def main() -> None:
    """Main function to run the MSD fitting analysis."""
    
    parser = argparse.ArgumentParser(
        description='Fit MSD data to linear model and extract Diffusion Coefficient',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python plot_msd_fit.py Data/Experiment8_spots_40minstep.csv
    python plot_msd_fit.py Data/Experiment8_spots_40minstep.csv --output msd_fit.png --tau-max 20
    python plot_msd_fit.py Data/Experiment8_spots_40minstep.csv --max-lag-fraction 0.5
        """
    )
    
    parser.add_argument(
        'csv',
        type=str,
        help='Path to the trajectories CSV file'
    )
    
    parser.add_argument(
        '--max-lag-fraction',
        type=float,
        default=None,
        help='Fraction (0<f<=1] of the longest track to cap maximum lag for MSD calculation. '
             'If omitted, uses full range (N_max-1).'
    )
    
    parser.add_argument(
        '--fit-fraction',
        type=float,
        default=0.10,
        help='Fraction of lag steps to use for linear fit (default: 0.10 = 10%%). '
             'Theoretically, the first ~10%% shows pure linear diffusion.'
    )
    
    parser.add_argument(
        '--D-initial',
        type=float,
        default=1e-2,
        help='Initial guess for Diffusion Coefficient in μm²/s (default: 1e-2)'
    )
    
    parser.add_argument(
        '--D-bounds',
        type=float,
        nargs=2,
        default=[9e-4, 1.5e-1],
        metavar=('LOWER', 'UPPER'),
        help='Bounds for D: lower and upper (default: 9e-4 1.5e-1 μm²/s). '
             'Based on Stokes-Einstein for 240nm particles with safety margins.'
    )
    
    parser.add_argument(
        '--output',
        type=str,
        default='msd_fit_plot.png',
        help='Output filename for the plot (default: msd_fit_plot.png)'
    )
    
    args = parser.parse_args()
    
    # Step 1: Read trajectories and compute ensemble-averaged MSD
    print(f"Reading trajectories from: {args.csv}")
    trajectories = read_trajectories_from_csv(args.csv)
    
    print(f"Computing ensemble-averaged MSD...")
    result = calculate_ensemble_msd(trajectories, max_lag_fraction=args.max_lag_fraction)
    
    if result.tau.size == 0:
        print("Error: No MSD data to fit (empty or insufficient points).")
        return
    
    # Print MSD calculation info
    print(f"\n{'='*60}")
    print(f"MSD Calculation Results:")
    print(f"{'='*60}")
    print(f"  Δt (global):               {result.dt:.6f} s")
    print(f"  n_max (steps):             {result.n_max}")
    print(f"  τ_max (seconds):           {result.tau[-1]:.2f} s")
    print(f"  Total trajectories (M):    {result.total_trajectories}")
    print(f"  Longest trajectory length: {result.longest_trajectory_points} points")
    print(f"{'='*60}\n")
    
    # Step 2: Fit MSD to linear model using first fit_fraction of steps
    n_fit_steps = int(args.fit_fraction * result.n_max)
    tau_fit_max = n_fit_steps * result.dt
    
    print(f"Fitting MSD to linear model: MSD(τ) = 2D·τ")
    print(f"  Using first {args.fit_fraction:.0%} of lag steps ({n_fit_steps}/{result.n_max} steps)")
    print(f"  Fit range: τ ≤ {tau_fit_max:.2f} s")
    print(f"  Initial guess: D = {args.D_initial:.2e} μm²/s")
    print(f"  Bounds: D ∈ [{args.D_bounds[0]:.2e}, {args.D_bounds[1]:.2e}] μm²/s")
    
    try:
        fit_result = fit_msd_linear(
            tau=result.tau,
            msd=result.msd,
            n_max=result.n_max,
            dt=result.dt,
            fit_fraction=args.fit_fraction,
            D_initial=args.D_initial,
            D_bounds=tuple(args.D_bounds),
        )
    except (ValueError, RuntimeError) as e:
        print(f"\nError during fitting: {e}")
        return
    
    # Print fit results
    print(f"\n{'='*60}")
    print(f"Fit Results:")
    print(f"{'='*60}")
    print(f"  Diffusion Coefficient (D): ({fit_result.D:.6e} ± {fit_result.D_error:.2e}) μm²/s")
    print(f"  R² (coefficient of determination): {fit_result.R_squared:.6f}")
    print(f"  Points used in fit: {fit_result.tau_fit.size}")
    print(f"  Fit τ range: [{fit_result.tau_fit.min():.2f}, {fit_result.tau_fit.max():.2f}] s")
    print(f"{'='*60}\n")
    
    # Step 3: Create and save plot
    output_path = Path(args.output)
    print(f"Creating plot...")
    
    plot_msd_with_fit(
        tau_fit=fit_result.tau_fit,
        msd_fit=fit_result.msd_fit,
        msd_fit_line=fit_result.msd_predicted,
        D=fit_result.D,
        D_error=fit_result.D_error,
        R_squared=fit_result.R_squared,
        output_path=output_path,
        fit_fraction=args.fit_fraction,
        n_fit_steps=n_fit_steps,
    )
    
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
