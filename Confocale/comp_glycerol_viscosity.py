# Python Code to compute the viscosity of a glycerol-water mixture and the
# diffusion coefficient of a spherical particle via the Einstein-Stokes relation.
# Valid for mass concentrations 0-100% and temperatures 0-100 °C.
# Viscosity model: Cheng (2008), Ind. Eng. Chem. Res. 47, 3285-3288
# Source: https://pubs.acs.org/doi/10.1021/ie071349z

import numpy as np


def glycerol_water_viscosity(T_C, c_m):
    """
    Compute the dynamic viscosity of a glycerol-water mixture.

    Parameters
    ----------
    T_C : float
        Temperature in degrees Celsius.
    c_m : float
        Mass fraction of glycerol (0 to 1).

    Returns
    -------
    eta_mix : float
        Dynamic viscosity of the mixture in Pa·s.
    """
    # Pure glycerol viscosity [Pa·s]
    eta_glycerol = 12.100 * np.exp((-1233 + T_C) * T_C / (9900 + 70 * T_C))
    # Pure water viscosity [Pa·s]
    eta_water = 1.790e-3 * np.exp((-1230 + T_C) * T_C / (36100 + 360 * T_C))

    # Mixing rule (Cheng 2008)
    a = 0.705 - 0.0017 * T_C
    b = (4.9 + 0.036 * T_C) * a ** 2.5
    # Weight factor for the solution
    alpha = 1 - c_m + (a * b * c_m * (1 - c_m)) / (a * c_m + b * (1 - c_m))
    A = np.log(eta_water / eta_glycerol)

    eta_mix = eta_glycerol * np.exp(A * alpha)
    return eta_mix


def diffusion_coefficient(T_C, c_m, r):
    """
    Compute the diffusion coefficient via the Einstein-Stokes relation:
        D = k_B T / (6 pi eta r)

    Parameters
    ----------
    T_C : float
        Temperature in degrees Celsius.
    c_m : float
        Mass fraction of glycerol (0 to 1).
    r : float
        Particle radius in metres.

    Returns
    -------
    D : float
        Diffusion coefficient in m^2/s.
    """
    k_B = 1.38e-23  # Boltzmann constant [J/K]
    T_K = T_C + 273.15

    eta = glycerol_water_viscosity(T_C, c_m)
    D = k_B * T_K / (6 * np.pi * eta * r)
    return D


if __name__ == "__main__":
    # --- Experimental parameters ---
    T_C = 23.3       # Temperature [°C]
    c_m = 0.59      # Glycerol mass fraction
    r = 120e-9       # Particle radius [m]
    k_B = 1.38e-23   # Boltzmann constant [J/K]
    T_K = T_C + 273.15

    eta = glycerol_water_viscosity(T_C, c_m)
    D = diffusion_coefficient(T_C, c_m, r)

    print("=" * 55)
    print("VISCOSITY (Cheng 2008) & DIFFUSION (Einstein-Stokes)")
    print("=" * 55)
    print(f"Temperature:          T = {T_C} °C = {T_K:.2f} K")
    print(f"Glycerol mass frac.:  c_m = {c_m}")
    print(f"Particle radius:      r = {r*1e9:.0f} nm")
    print(f"k_B:                  {k_B} J/K")
    print()
    print(f"eta_mix ({c_m*100:.0f}%, {T_C}°C) = {eta:.6e} Pa·s    ({eta*1e3:.2f} mPa·s)")
    print()
    print(f"D = k_B T / (6 pi eta r)")
    print(f"  = {D:.6e} m²/s")
    print(f"  = {D*1e12:.6f} µm²/s")

