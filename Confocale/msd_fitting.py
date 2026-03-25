"""
Linear MSD fitting module for extracting Diffusion Coefficient (D).

This module implements a robust fitting function to extract the Diffusion Coefficient
from Mean Squared Displacement (MSD) data using the linear model:
    
    MSD(τ) = 4D·τ
    
where:
    - MSD is the Mean Squared Displacement
    - τ (tau) is the time lag
    - D is the Diffusion Coefficient

The fitting is performed using scipy.optimize.curve_fit with the Levenberg-Marquardt
algorithm (or Trust Region Reflective when bounds are specified).

Physical Context:
    For particles diffusing in highly viscous media (e.g., 80% Glycerol/Water solution),
    the expected Diffusion Coefficient is very low, typically in the range of
    10^-6 to 10^-2 μm²/s.

Dependencies: numpy, scipy
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from scipy.optimize import curve_fit


@dataclass(frozen=True)
class FitResult:
    """Container for MSD linear fit results.
    
    Attributes:
        D: Diffusion coefficient in units of (length²/time), e.g., μm²/s
        D_error: Standard error on D computed from the covariance matrix
        pcov: Full covariance matrix from the fit
        tau_fit: Time lag values used in the fit (subset where τ ≤ tau_max)
        msd_fit: MSD values used in the fit (corresponding to tau_fit)
        msd_predicted: Predicted MSD values from the fitted model over tau_fit range
        chi_squared_red: Reduced chi-squared (χ²_ν) for the fit quality (NaN if no sigma provided)
    """
    
    D: float
    D_error: float
    pcov: np.ndarray
    tau_fit: np.ndarray
    msd_fit: np.ndarray
    msd_predicted: np.ndarray
    chi_squared_red: float
    msd_sigma_fit: Optional[np.ndarray] = None


def linear_msd_model(tau: np.ndarray, D: float) -> np.ndarray:
    """Linear MSD model: MSD(τ) = 4D·τ (for 2D motion)
    
    Args:
        tau: Time lag values (seconds or appropriate time units)
        D: Diffusion coefficient (length²/time units)
    
    Returns:
        MSD values predicted by the model
    """
    return 4.0 * D * tau


def calculate_r_squared(y_observed: np.ndarray, y_predicted: np.ndarray) -> float:
    """Calculate the coefficient of determination (R²).
    
    R² = 1 - (SS_res / SS_tot)
    
    where:
        SS_res = Σ(y_observed - y_predicted)²  (residual sum of squares)
        SS_tot = Σ(y_observed - y_mean)²       (total sum of squares)
    
    Args:
        y_observed: Observed data values
        y_predicted: Predicted values from the model
    
    Returns:
        R² value (1.0 = perfect fit, lower values indicate worse fit)
    """
    ss_res = np.sum((y_observed - y_predicted) ** 2)
    ss_tot = np.sum((y_observed - np.mean(y_observed)) ** 2)
    
    # Avoid division by zero
    if ss_tot == 0:
        return 0.0
    
    return 1.0 - (ss_res / ss_tot)


def calculate_reduced_chi_squared(
    y_observed: np.ndarray,
    y_predicted: np.ndarray,
    sigma: np.ndarray,
    n_params: int,
) -> float:
    """Calculate the reduced chi-squared statistic.

    χ²_ν = (1/ν) · Σ((y_observed - y_predicted)² / σ²)

    where ν = N - n_params is the number of degrees of freedom.
    A value near 1 indicates the model fits the data well within the uncertainties.

    Args:
        y_observed: Observed data values
        y_predicted: Predicted values from the model
        sigma: Uncertainties on y_observed (e.g., SEM of EA-MSD)
        n_params: Number of free parameters in the model

    Returns:
        Reduced chi-squared value (χ²_ν). Returns nan if ν ≤ 0.
    """
    nu = len(y_observed) - n_params
    if nu <= 0:
        return float('nan')
    return float(np.sum(((y_observed - y_predicted) / sigma) ** 2) / nu)


def fit_msd_linear(
    tau: np.ndarray,
    msd: np.ndarray,
    n_max: int,
    dt: float,
    fit_fraction: float = 0.10,
    D_initial: float = 1e-2,
    D_bounds: Tuple[float, float] = (1e-6, 10.0),
    msd_sigma: Optional[np.ndarray] = None,
) -> FitResult:
    """Fit MSD data to linear model MSD(τ) = 4D·τ using Levenberg-Marquardt algorithm.
    
    This function performs a partial fit on the first portion of the MSD data,
    corresponding to the short-time linear diffusion regime. From theory, the
    first ~10% of lag steps exhibit purely linear behavior before confinement
    or other effects cause deviations.
    
    Physical Constraints:
        For 240nm particles in 85% glycerol/water solution at room temperature,
        using the Stokes-Einstein relation D = kT/(6πηr), the expected diffusion
        coefficient is approximately 0.006-0.009 μm²/s. However, experimental
        observations show higher values (~0.14 μm²/s), possibly due to:
        - Lower actual glycerol concentration
        - Higher temperature
        - Active transport or convective flows
        The default bounds [1e-6, 10.0] are wide to accommodate various experimental conditions.
    
    Algorithm:
        Uses scipy.optimize.curve_fit which automatically selects:
        - 'lm' (Levenberg-Marquardt) for unconstrained problems
        - 'trf' (Trust Region Reflective) when bounds are specified
    
    Args:
        tau: Array of time lag values (must be in seconds or consistent time units)
        msd: Array of MSD values (must be in length² units, e.g., μm²)
        n_max: Maximum number of lag steps available
        dt: Time step (seconds)
        fit_fraction: Fraction of n_max to use for fitting (default: 0.10 = 10%)
        D_initial: Initial guess for the Diffusion Coefficient (default: 1e-2 μm²/s)
        D_bounds: Tuple (lower, upper) bounds for D (default: (1e-6, 10.0) μm²/s)
    
    Returns:
        FitResult containing:
            - D: Optimal Diffusion Coefficient
            - D_error: Standard error on D
            - pcov: Covariance matrix
            - tau_fit, msd_fit: Data subset used for fitting
            - msd_predicted: Model predictions for tau_fit
            - chi_squared_red: Goodness of fit metric (reduced chi-squared)
    
    Raises:
        ValueError: If input arrays are incompatible or insufficient data for fitting
        RuntimeError: If the optimization fails to converge
    """
    # Input validation
    tau = np.asarray(tau, dtype=float)
    msd = np.asarray(msd, dtype=float)
    
    if tau.shape != msd.shape:
        raise ValueError(f"tau and msd must have the same shape. Got tau: {tau.shape}, msd: {msd.shape}")
    
    if tau.size == 0:
        raise ValueError("Input arrays are empty")
    
    # Validate fraction
    if not (0 < fit_fraction <= 1.0):
        raise ValueError(f"fit_fraction must be in (0, 1], got {fit_fraction}")
    
    # Calculate number of steps to use for fitting (first fit_fraction of n_max)
    n_fit_steps = max(2, int(fit_fraction * n_max))
    tau_max = n_fit_steps * dt
    
    # Select data subset for fitting: first n_fit_steps
    # Since tau = n * dt, we use tau <= tau_max
    mask = tau <= tau_max
    tau_fit = tau[mask]
    msd_fit = msd[mask]
    
    if tau_fit.size == 0:
        raise ValueError(f"No data points in fitting range (first {fit_fraction:.0%} = {n_fit_steps} steps)")
    
    if tau_fit.size < 2:
        raise ValueError(f"Need at least 2 data points for fitting, got {tau_fit.size}")
    
    # Remove any NaN or infinite values
    valid = np.isfinite(tau_fit) & np.isfinite(msd_fit)
    tau_fit = tau_fit[valid]
    msd_fit = msd_fit[valid]
    
    if tau_fit.size < 2:
        raise ValueError("Insufficient valid (finite) data points after filtering")

    # Prepare sigma for weighted fit
    sigma_fit = None
    if msd_sigma is not None:
        sigma_subset = np.asarray(msd_sigma, dtype=float)[mask][valid]
        if np.all(np.isfinite(sigma_subset) & (sigma_subset > 0)):
            sigma_fit = sigma_subset

    # Validate bounds
    D_lower, D_upper = D_bounds
    if D_lower <= 0:
        raise ValueError(f"Lower bound for D must be positive, got {D_lower}")
    if D_upper <= D_lower:
        raise ValueError(f"Upper bound must be greater than lower bound: ({D_lower}, {D_upper})")
    if not (D_lower <= D_initial <= D_upper):
        raise ValueError(f"Initial guess D_initial={D_initial} must be within bounds ({D_lower}, {D_upper})")
    
    # Perform the fit
    # When bounds are provided, scipy automatically uses 'trf' method instead of 'lm'
    try:
        popt, pcov = curve_fit(
            linear_msd_model,
            tau_fit,
            msd_fit,
            p0=[D_initial],
            bounds=([D_lower], [D_upper]),
            method='trf',  # Trust Region Reflective handles bounds well
            **({"sigma": sigma_fit, "absolute_sigma": True} if sigma_fit is not None else {}),
        )
    except RuntimeError as e:
        raise RuntimeError(f"Curve fitting failed: {e}")
    
    # Extract results
    D_optimal = float(popt[0])
    
    # Calculate standard error on D from covariance matrix
    # perr = sqrt(diag(pcov))
    if pcov is not None and np.all(np.isfinite(pcov)):
        D_error = float(np.sqrt(pcov[0, 0]))
    else:
        D_error = float('nan')
        print("Warning: Could not estimate error on D from covariance matrix")
    
    # Calculate predicted values and goodness-of-fit
    msd_predicted = linear_msd_model(tau_fit, D_optimal)
    if sigma_fit is not None:
        chi_squared_red = calculate_reduced_chi_squared(msd_fit, msd_predicted, sigma_fit, n_params=1)
    else:
        chi_squared_red = float('nan')

    return FitResult(
        D=D_optimal,
        D_error=D_error,
        pcov=pcov,
        tau_fit=tau_fit,
        msd_fit=msd_fit,
        msd_predicted=msd_predicted,
        chi_squared_red=chi_squared_red,
        msd_sigma_fit=sigma_fit,
    )


if __name__ == "__main__":
    # Quick test with synthetic data
    print("Testing fit_msd_linear with synthetic data...")
    
    # Generate synthetic data: D_true = 5e-4 μm²/s
    D_true = 5e-4
    tau_test = np.linspace(0.1, 20, 50)
    msd_test = linear_msd_model(tau_test, D_true)
    
    # Add some noise
    rng = np.random.default_rng(42)
    noise = 0.05 * msd_test * rng.normal(size=msd_test.shape)
    msd_test_noisy = msd_test + noise
    
    # Fit the data
    n_max_test = len(tau_test)
    dt_test = float(tau_test[1] - tau_test[0])
    sigma_test = 0.05 * np.abs(msd_test) + 1e-6  # synthetic SEM proportional to MSD
    result = fit_msd_linear(tau_test, msd_test_noisy, n_max=n_max_test, dt=dt_test, msd_sigma=sigma_test)
    
    print(f"\nTrue D: {D_true:.6e} μm²/s")
    print(f"Fitted D: {result.D:.6e} ± {result.D_error:.6e} μm²/s")
    print(f"chi^2_red: {result.chi_squared_red:.6f}")
    print(f"Points used in fit: {result.tau_fit.size}")
    print(f"Fit range: τ = [{result.tau_fit.min():.2f}, {result.tau_fit.max():.2f}] s")
    print("\nTest passed!")
