"""
MSD fitting module — all four diffusion models.

Models:
    1.  Linear:           MSD(τ) = 4D·τ
    1b. Linear+offset:   MSD(τ) = 4D·τ + c   (c = 4σ², localization error)
    2.  Nonlinear (drift): MSD(τ) = 4D·τ + v²·τ²
    3.  Anomalous:        MSD(τ) = 4D_α·τ^α
    4.  Anomalous+drift:  MSD(τ) = 4D_α·τ^α + v²·τ²

Each model provides:
    - A model function
    - A FitResult dataclass
    - A fit function (with optional interval optimisation for models 2-4)

Shared utilities: calculate_r_squared, calculate_reduced_chi_squared, calculate_rss,
    velocity analysis (VelocityStats, compute_trajectory_velocity, analyze_velocities).

Dependencies: numpy, scipy
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Mapping, Optional, Tuple, Union

import numpy as np
from scipy.optimize import curve_fit

from data_reader import Trajectory


# ---------------------------------------------------------------------------
# Shared utilities
# ---------------------------------------------------------------------------

def calculate_r_squared(y_observed: np.ndarray, y_predicted: np.ndarray) -> float:
    """R² = 1 − SS_res / SS_tot."""
    ss_res = np.sum((y_observed - y_predicted) ** 2)
    ss_tot = np.sum((y_observed - np.mean(y_observed)) ** 2)
    if ss_tot == 0:
        return 0.0
    return 1.0 - (ss_res / ss_tot)


def calculate_reduced_chi_squared(
    y_observed: np.ndarray,
    y_predicted: np.ndarray,
    sigma: np.ndarray,
    n_params: int,
) -> float:
    """χ²_ν = (1/ν) Σ ((y_obs − y_pred) / σ)²  with ν = N − n_params."""
    nu = len(y_observed) - n_params
    if nu <= 0:
        return float("nan")
    return float(np.sum(((y_observed - y_predicted) / sigma) ** 2) / nu)


def calculate_rss(y_observed: np.ndarray, y_predicted: np.ndarray) -> float:
    """Residual sum of squares: Σ(y_obs − y_pred)²."""
    return float(np.sum((y_observed - y_predicted) ** 2))


# ---------------------------------------------------------------------------
# Velocity analysis (used by nonlinear and anomalous+drift models)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class VelocityStats:
    """Statistics of trajectory velocities (v = path_length / duration)."""
    velocities: np.ndarray
    mean: float
    median: float
    std: float
    v_initial: float
    v_bounds: Tuple[float, float]
    n_trajectories_total: int
    n_trajectories_used: int


def compute_trajectory_velocity(traj: Trajectory) -> float:
    """v = total_path_length / total_time for a single trajectory."""
    if traj.n_points < 2:
        return float("nan")
    dx = np.diff(traj.x)
    dy = np.diff(traj.y)
    total_path = float(np.sum(np.sqrt(dx**2 + dy**2)))
    duration = float(traj.time[-1] - traj.time[0])
    if duration <= 0:
        return float("nan")
    return total_path / duration


def analyze_velocities(
    trajectories: Mapping[Union[int, str], Trajectory],
    min_points: int = 30,
) -> VelocityStats:
    """Compute velocity statistics from trajectories with ≥ min_points."""
    n_total = len(trajectories)
    vels = [
        v
        for traj in trajectories.values()
        if traj.n_points >= min_points
        for v in [compute_trajectory_velocity(traj)]
        if np.isfinite(v) and v > 0
    ]
    if not vels:
        raise ValueError(
            f"No valid trajectories with ≥{min_points} points. Total: {n_total}"
        )
    velocities = np.array(vels)
    mean_v = float(np.mean(velocities))
    median_v = float(np.median(velocities))
    std_v = float(np.std(velocities))
    return VelocityStats(
        velocities=velocities,
        mean=mean_v,
        median=median_v,
        std=std_v,
        v_initial=median_v,
        v_bounds=(0.0, mean_v + 3.0 * std_v),
        n_trajectories_total=n_total,
        n_trajectories_used=len(velocities),
    )


# ---------------------------------------------------------------------------
# Helper: select data subset for fitting
# ---------------------------------------------------------------------------

def _select_fit_data(
    tau: np.ndarray,
    msd: np.ndarray,
    n_max: int,
    dt: float,
    fit_fraction: float,
    min_points: int,
    msd_sigma: Optional[np.ndarray] = None,
    require_positive_tau: bool = False,
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """Return (tau_fit, msd_fit, sigma_fit) after fraction-based subsetting."""
    n_fit = max(2, int(fit_fraction * n_max))
    tau_max = n_fit * dt
    mask = tau <= tau_max
    tau_f = tau[mask]
    msd_f = msd[mask]

    valid = np.isfinite(tau_f) & np.isfinite(msd_f)
    if require_positive_tau:
        valid &= tau_f > 0
    tau_f = tau_f[valid]
    msd_f = msd_f[valid]

    if tau_f.size < min_points:
        raise ValueError(
            f"Insufficient points ({tau_f.size}) for fit at fraction {fit_fraction}"
        )

    sigma_f = None
    if msd_sigma is not None:
        s = np.asarray(msd_sigma, dtype=float)[mask][valid]
        good_sigma = np.isfinite(s) & (s > 0)
        if np.all(good_sigma):
            sigma_f = s
        elif np.sum(good_sigma) >= min_points:
            # Keep only the points that have valid sigma
            tau_f = tau_f[good_sigma]
            msd_f = msd_f[good_sigma]
            sigma_f = s[good_sigma]

    return tau_f, msd_f, sigma_f


# ===================================================================
# Model 1 — Linear: MSD = 4D·τ
# ===================================================================

@dataclass(frozen=True)
class FitResult:
    """Result of a linear MSD fit."""
    D: float
    D_error: float
    pcov: np.ndarray
    tau_fit: np.ndarray
    msd_fit: np.ndarray
    msd_predicted: np.ndarray
    chi_squared_red: float
    msd_sigma_fit: Optional[np.ndarray] = None


def linear_msd_model(tau: np.ndarray, D: float) -> np.ndarray:
    """MSD(τ) = 4D·τ"""
    return 4.0 * D * tau


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
    """Fit MSD to linear model using the first *fit_fraction* of lag steps."""
    tau = np.asarray(tau, dtype=float)
    msd = np.asarray(msd, dtype=float)
    if tau.shape != msd.shape or tau.size == 0:
        raise ValueError("tau and msd must be non-empty arrays of the same shape")
    if not (0 < fit_fraction <= 1.0):
        raise ValueError(f"fit_fraction must be in (0, 1], got {fit_fraction}")

    tau_f, msd_f, sigma_f = _select_fit_data(
        tau, msd, n_max, dt, fit_fraction, min_points=2, msd_sigma=msd_sigma,
    )

    popt, pcov = curve_fit(
        linear_msd_model, tau_f, msd_f,
        p0=[D_initial],
        bounds=([D_bounds[0]], [D_bounds[1]]),
        method="trf",
        **({"sigma": sigma_f, "absolute_sigma": True} if sigma_f is not None else {}),
    )

    D_opt = float(popt[0])
    D_err = float(np.sqrt(pcov[0, 0])) if pcov is not None and np.all(np.isfinite(pcov)) else float("nan")
    msd_pred = linear_msd_model(tau_f, D_opt)
    chi2 = calculate_reduced_chi_squared(msd_f, msd_pred, sigma_f, 1) if sigma_f is not None else float("nan")

    return FitResult(
        D=D_opt, D_error=D_err, pcov=pcov,
        tau_fit=tau_f, msd_fit=msd_f, msd_predicted=msd_pred,
        chi_squared_red=chi2, msd_sigma_fit=sigma_f,
    )


# ===================================================================
# Model 1b — Linear with offset: MSD = 4D·τ + c   (c = 4σ²)
# ===================================================================

@dataclass(frozen=True)
class LinearOffsetFitResult:
    """Result of a linear+offset MSD fit (localization error)."""
    D: float
    D_error: float
    offset: float            # c in MSD = 4Dτ + c
    offset_error: float
    sigma_loc: float         # σ = sqrt(c/4) if c > 0, else NaN
    sigma_loc_error: float   # propagated error
    pcov: np.ndarray
    tau_fit: np.ndarray
    msd_fit: np.ndarray
    msd_predicted: np.ndarray
    chi_squared_red: float
    msd_sigma_fit: Optional[np.ndarray] = None


def linear_offset_msd_model(tau: np.ndarray, D: float, c: float) -> np.ndarray:
    """MSD(τ) = 4D·τ + c"""
    return 4.0 * D * tau + c


def fit_msd_linear_offset(
    tau: np.ndarray,
    msd: np.ndarray,
    n_max: int,
    dt: float,
    fit_fraction: float = 0.10,
    D_initial: float = 1e-2,
    D_bounds: Tuple[float, float] = (1e-6, 10.0),
    c_initial: float = 0.0,
    c_bounds: Tuple[float, float] = (-1.0, 10.0),
    msd_sigma: Optional[np.ndarray] = None,
) -> LinearOffsetFitResult:
    """Fit MSD to linear+offset model: MSD = 4D·τ + c  (c = 4σ²)."""
    tau = np.asarray(tau, dtype=float)
    msd = np.asarray(msd, dtype=float)
    if tau.shape != msd.shape or tau.size == 0:
        raise ValueError("tau and msd must be non-empty arrays of the same shape")
    if not (0 < fit_fraction <= 1.0):
        raise ValueError(f"fit_fraction must be in (0, 1], got {fit_fraction}")

    tau_f, msd_f, sigma_f = _select_fit_data(
        tau, msd, n_max, dt, fit_fraction, min_points=3, msd_sigma=msd_sigma,
    )

    popt, pcov = curve_fit(
        linear_offset_msd_model, tau_f, msd_f,
        p0=[D_initial, c_initial],
        bounds=([D_bounds[0], c_bounds[0]], [D_bounds[1], c_bounds[1]]),
        method="trf",
        **({"sigma": sigma_f, "absolute_sigma": True} if sigma_f is not None else {}),
    )

    D_opt = float(popt[0])
    c_opt = float(popt[1])
    fin = pcov is not None and np.all(np.isfinite(pcov))
    D_err = float(np.sqrt(pcov[0, 0])) if fin else float("nan")
    c_err = float(np.sqrt(pcov[1, 1])) if fin else float("nan")

    # Derive localization error: σ = sqrt(c/4)
    if c_opt > 0:
        sigma_loc = float(np.sqrt(c_opt / 4.0))
        # dσ/dc = 1/(4σ)  →  σ_σ = c_err / (4σ)
        sigma_loc_err = c_err / (4.0 * sigma_loc) if sigma_loc > 0 else float("nan")
    else:
        sigma_loc = float("nan")
        sigma_loc_err = float("nan")

    msd_pred = linear_offset_msd_model(tau_f, D_opt, c_opt)
    chi2 = (calculate_reduced_chi_squared(msd_f, msd_pred, sigma_f, 2)
            if sigma_f is not None else float("nan"))

    return LinearOffsetFitResult(
        D=D_opt, D_error=D_err,
        offset=c_opt, offset_error=c_err,
        sigma_loc=sigma_loc, sigma_loc_error=sigma_loc_err,
        pcov=pcov,
        tau_fit=tau_f, msd_fit=msd_f, msd_predicted=msd_pred,
        chi_squared_red=chi2, msd_sigma_fit=sigma_f,
    )


# ===================================================================
# Model 2 — Nonlinear (drift): MSD = 4D·τ + v²·τ²
# ===================================================================

@dataclass(frozen=True)
class NonlinearFitResult:
    """Result of a nonlinear (drift) MSD fit."""
    D: float
    D_error: float
    v: float
    v_error: float
    pcov: np.ndarray
    tau_fit: np.ndarray
    msd_fit: np.ndarray
    msd_predicted: np.ndarray
    chi_squared_red: float
    RSS: float
    optimal_fraction: float
    n_fit_steps: int
    interval_results: Dict[float, Tuple[float, float]]
    msd_sigma_fit: Optional[np.ndarray] = None


def nonlinear_msd_model(tau: np.ndarray, D: float, v: float) -> np.ndarray:
    """MSD(τ) = 4D·τ + v²·τ²"""
    return 4.0 * D * tau + v**2 * tau**2


def _fit_nonlinear_at_fraction(
    tau, msd, n_max, dt, frac,
    D_initial, D_bounds, v_initial, v_bounds,
    msd_sigma,
):
    tau_f, msd_f, sigma_f = _select_fit_data(
        tau, msd, n_max, dt, frac, min_points=3, msd_sigma=msd_sigma,
    )
    popt, pcov = curve_fit(
        nonlinear_msd_model, tau_f, msd_f,
        p0=[D_initial, v_initial],
        bounds=([D_bounds[0], v_bounds[0]], [D_bounds[1], v_bounds[1]]),
        method="trf",
        **({"sigma": sigma_f, "absolute_sigma": True} if sigma_f is not None else {}),
    )
    D, v = float(popt[0]), float(popt[1])
    D_err = float(np.sqrt(pcov[0, 0])) if np.all(np.isfinite(pcov)) else float("nan")
    v_err = float(np.sqrt(pcov[1, 1])) if np.all(np.isfinite(pcov)) else float("nan")
    msd_pred = nonlinear_msd_model(tau_f, D, v)
    rss = calculate_rss(msd_f, msd_pred)
    chi2 = calculate_reduced_chi_squared(msd_f, msd_pred, sigma_f, 2) if sigma_f is not None else float("nan")
    return D, D_err, v, v_err, pcov, tau_f, msd_f, msd_pred, chi2, rss, sigma_f


def fit_msd_nonlinear(
    tau: np.ndarray,
    msd: np.ndarray,
    n_max: int,
    dt: float,
    velocity_stats: VelocityStats,
    D_initial: float = 1e-2,
    D_bounds: Tuple[float, float] = (9e-4, 1.5e-1),
    interval_step: float = 0.10,
    msd_sigma: Optional[np.ndarray] = None,
) -> NonlinearFitResult:
    """Fit nonlinear model testing intervals from 10 % to 90 %, pick best χ²_ν."""
    fractions = np.arange(interval_step, 1.0, interval_step)
    results: Dict[float, tuple] = {}
    sigma_fits: Dict[float, Optional[np.ndarray]] = {}
    interval_results: Dict[float, Tuple[float, float]] = {}

    print(f"\nTesting intervals {fractions[0]:.0%}–{fractions[-1]:.0%} …")
    for frac in fractions:
        try:
            D, D_err, v, v_err, pcov, tau_f, msd_f, msd_p, chi2, rss, sig = _fit_nonlinear_at_fraction(
                tau, msd, n_max, dt, frac, D_initial, D_bounds,
                velocity_stats.v_initial, velocity_stats.v_bounds, msd_sigma,
            )
            results[frac] = (D, D_err, v, v_err, chi2, pcov, tau_f, msd_f, msd_p, rss, int(frac * n_max))
            sigma_fits[frac] = sig
            interval_results[frac] = (chi2, rss)
            print(f"  {frac:4.0%}: χ²_ν={chi2:.4f}, RSS={rss:.3e}, D={D:.3e}, v={v:.3e}")
        except (ValueError, RuntimeError) as e:
            print(f"  {frac:4.0%}: Failed – {e}")

    if not results:
        raise RuntimeError("No valid fits across any interval")

    optimal = _pick_best(results, chi_idx=4, rss_idx=9)
    D, D_err, v, v_err, chi2, pcov, tau_f, msd_f, msd_p, rss, n_steps = results[optimal]
    print(f"\nOptimal: {optimal:.0%} (χ²_ν={chi2:.4f})")

    return NonlinearFitResult(
        D=D, D_error=D_err, v=v, v_error=v_err, pcov=pcov,
        tau_fit=tau_f, msd_fit=msd_f, msd_predicted=msd_p,
        chi_squared_red=chi2, RSS=rss,
        optimal_fraction=optimal, n_fit_steps=n_steps,
        interval_results=interval_results, msd_sigma_fit=sigma_fits.get(optimal),
    )


# ===================================================================
# Model 3 — Anomalous: MSD = 4D_α·τ^α
# ===================================================================

@dataclass(frozen=True)
class AnomalousFitResult:
    """Result of an anomalous MSD fit (no drift)."""
    D_alpha: float
    D_alpha_error: float
    alpha: float
    alpha_error: float
    pcov: np.ndarray
    tau_fit: np.ndarray
    msd_fit: np.ndarray
    msd_predicted: np.ndarray
    chi_squared_red: float
    RSS: float
    optimal_fraction: float
    n_fit_steps: int
    interval_results: Dict[float, Tuple[float, float]]


def anomalous_msd_model(tau: np.ndarray, D_alpha: float, alpha: float) -> np.ndarray:
    """MSD(τ) = 4D_α·τ^α"""
    return 4.0 * D_alpha * np.power(tau, alpha)


def _fit_anomalous_at_fraction(
    tau, msd, n_max, dt, frac,
    D_alpha_initial, D_alpha_bounds, alpha_initial, alpha_bounds,
    msd_sigma,
):
    tau_f, msd_f, sigma_f = _select_fit_data(
        tau, msd, n_max, dt, frac, min_points=3, msd_sigma=msd_sigma,
        require_positive_tau=True,
    )
    popt, pcov = curve_fit(
        anomalous_msd_model, tau_f, msd_f,
        p0=[D_alpha_initial, alpha_initial],
        bounds=([D_alpha_bounds[0], alpha_bounds[0]], [D_alpha_bounds[1], alpha_bounds[1]]),
        method="trf", maxfev=5000,
        **({"sigma": sigma_f, "absolute_sigma": True} if sigma_f is not None else {}),
    )
    Da, a = float(popt[0]), float(popt[1])
    Da_err = float(np.sqrt(pcov[0, 0])) if np.all(np.isfinite(pcov)) else float("nan")
    a_err = float(np.sqrt(pcov[1, 1])) if np.all(np.isfinite(pcov)) else float("nan")
    msd_pred = anomalous_msd_model(tau_f, Da, a)
    rss = calculate_rss(msd_f, msd_pred)
    chi2 = calculate_reduced_chi_squared(msd_f, msd_pred, sigma_f, 2) if sigma_f is not None else float("nan")
    return Da, Da_err, a, a_err, pcov, tau_f, msd_f, msd_pred, chi2, rss


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
    msd_sigma: Optional[np.ndarray] = None,
) -> AnomalousFitResult:
    """Fit anomalous model testing intervals 10 %–90 %, pick best χ²_ν."""
    fractions = np.arange(interval_step, 1.0, interval_step)
    results: Dict[float, tuple] = {}
    interval_results: Dict[float, Tuple[float, float]] = {}

    print(f"\nTesting intervals {fractions[0]:.0%}–{fractions[-1]:.0%} …")
    for frac in fractions:
        try:
            Da, Da_err, a, a_err, pcov, tau_f, msd_f, msd_p, chi2, rss = _fit_anomalous_at_fraction(
                tau, msd, n_max, dt, frac,
                D_alpha_initial, D_alpha_bounds, alpha_initial, alpha_bounds, msd_sigma,
            )
            results[frac] = (Da, Da_err, a, a_err, chi2, pcov, tau_f, msd_f, msd_p, rss, int(frac * n_max))
            interval_results[frac] = (chi2, rss)
            print(f"  {frac:4.0%}: χ²_ν={chi2:.4f}, RSS={rss:.3e}, D_α={Da:.3e}, α={a:.4f}")
        except (ValueError, RuntimeError) as e:
            print(f"  {frac:4.0%}: Failed – {e}")

    if not results:
        raise RuntimeError("No valid fits across any interval")

    optimal = _pick_best(results, chi_idx=4, rss_idx=9)
    Da, Da_err, a, a_err, chi2, pcov, tau_f, msd_f, msd_p, rss, n_steps = results[optimal]
    print(f"\nOptimal: {optimal:.0%} (χ²_ν={chi2:.4f})")

    return AnomalousFitResult(
        D_alpha=Da, D_alpha_error=Da_err, alpha=a, alpha_error=a_err, pcov=pcov,
        tau_fit=tau_f, msd_fit=msd_f, msd_predicted=msd_p,
        chi_squared_red=chi2, RSS=rss,
        optimal_fraction=optimal, n_fit_steps=n_steps,
        interval_results=interval_results,
    )


# ===================================================================
# Model 4 — Anomalous + drift: MSD = 4D_α·τ^α + v²·τ²
# ===================================================================

@dataclass(frozen=True)
class DriftAnomalousFitResult:
    """Result of an anomalous+drift MSD fit."""
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
    chi_squared_red: float
    RSS: float
    optimal_fraction: float
    n_fit_steps: int
    interval_results: Dict[float, Tuple[float, float]]


def anomalous_drift_msd_model(tau: np.ndarray, D_alpha: float, alpha: float, v: float) -> np.ndarray:
    """MSD(τ) = 4D_α·τ^α + v²·τ²"""
    return 4.0 * D_alpha * np.power(tau, alpha) + v**2 * tau**2


def _fit_anomalous_drift_at_fraction(
    tau, msd, n_max, dt, frac,
    D_alpha_initial, D_alpha_bounds, alpha_initial, alpha_bounds,
    v_initial, v_bounds, msd_sigma,
):
    tau_f, msd_f, sigma_f = _select_fit_data(
        tau, msd, n_max, dt, frac, min_points=4, msd_sigma=msd_sigma,
        require_positive_tau=True,
    )
    popt, pcov = curve_fit(
        anomalous_drift_msd_model, tau_f, msd_f,
        p0=[D_alpha_initial, alpha_initial, v_initial],
        bounds=(
            [D_alpha_bounds[0], alpha_bounds[0], v_bounds[0]],
            [D_alpha_bounds[1], alpha_bounds[1], v_bounds[1]],
        ),
        method="trf", maxfev=10000,
        **({"sigma": sigma_f, "absolute_sigma": True} if sigma_f is not None else {}),
    )
    Da, a, v = float(popt[0]), float(popt[1]), float(popt[2])
    Da_err = float(np.sqrt(pcov[0, 0])) if np.all(np.isfinite(pcov)) else float("nan")
    a_err = float(np.sqrt(pcov[1, 1])) if np.all(np.isfinite(pcov)) else float("nan")
    v_err = float(np.sqrt(pcov[2, 2])) if np.all(np.isfinite(pcov)) else float("nan")
    msd_pred = anomalous_drift_msd_model(tau_f, Da, a, v)
    rss = calculate_rss(msd_f, msd_pred)
    chi2 = calculate_reduced_chi_squared(msd_f, msd_pred, sigma_f, 3) if sigma_f is not None else float("nan")
    return Da, Da_err, a, a_err, v, v_err, pcov, tau_f, msd_f, msd_pred, chi2, rss


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
    msd_sigma: Optional[np.ndarray] = None,
) -> DriftAnomalousFitResult:
    """Fit anomalous+drift model testing intervals 10 %–90 %, pick best χ²_ν."""
    fractions = np.arange(interval_step, 1.0, interval_step)
    results: Dict[float, tuple] = {}
    interval_results: Dict[float, Tuple[float, float]] = {}

    print(f"\nTesting intervals {fractions[0]:.0%}–{fractions[-1]:.0%} …")
    for frac in fractions:
        try:
            Da, Da_err, a, a_err, v, v_err, pcov, tau_f, msd_f, msd_p, chi2, rss = _fit_anomalous_drift_at_fraction(
                tau, msd, n_max, dt, frac,
                D_alpha_initial, D_alpha_bounds, alpha_initial, alpha_bounds,
                velocity_stats.v_initial, velocity_stats.v_bounds, msd_sigma,
            )
            results[frac] = (Da, Da_err, a, a_err, v, v_err, chi2, pcov, tau_f, msd_f, msd_p, rss, int(frac * n_max))
            interval_results[frac] = (chi2, rss)
            print(f"  {frac:4.0%}: χ²_ν={chi2:.4f}, RSS={rss:.3e}, D_α={Da:.3e}, α={a:.4f}, v={v:.3e}")
        except (ValueError, RuntimeError) as e:
            print(f"  {frac:4.0%}: Failed – {e}")

    if not results:
        raise RuntimeError("No valid fits across any interval")

    optimal = _pick_best_drift(results)
    Da, Da_err, a, a_err, v, v_err, chi2, pcov, tau_f, msd_f, msd_p, rss, n_steps = results[optimal]
    print(f"\nOptimal: {optimal:.0%} (χ²_ν={chi2:.4f})")

    return DriftAnomalousFitResult(
        D_alpha=Da, D_alpha_error=Da_err, alpha=a, alpha_error=a_err,
        v=v, v_error=v_err, pcov=pcov,
        tau_fit=tau_f, msd_fit=msd_f, msd_predicted=msd_p,
        chi_squared_red=chi2, RSS=rss,
        optimal_fraction=optimal, n_fit_steps=n_steps,
        interval_results=interval_results,
    )


# ---------------------------------------------------------------------------
# Internal: pick the best interval
# ---------------------------------------------------------------------------

def _pick_best(results: dict, chi_idx: int, rss_idx: int) -> float:
    """Select fraction with minimum χ²_ν (fallback: minimum RSS)."""
    chi_vals = {f: r[chi_idx] for f, r in results.items()}
    finite = {f: c for f, c in chi_vals.items() if np.isfinite(c)}
    if finite:
        return min(finite, key=finite.get)
    return min(results, key=lambda f: results[f][rss_idx])


def _pick_best_drift(results: dict) -> float:
    """Same logic for the anomalous+drift tuple layout (chi at idx 6, rss at idx 11)."""
    chi_vals = {f: r[6] for f, r in results.items()}
    finite = {f: c for f, c in chi_vals.items() if np.isfinite(c)}
    if finite:
        return min(finite, key=finite.get)
    return min(results, key=lambda f: results[f][11])
