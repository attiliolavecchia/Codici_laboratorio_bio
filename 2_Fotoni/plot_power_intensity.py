import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


DEFAULT_OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "img")

def read_data(filename):
    """
    Read data from Excel file (.xlsx).
    
    Parameters:
    -----------
    filename : str
        Path to the data file (.xlsx)
        
    Returns:
    --------
    power : array
        Power values (mW)
    intensity : array
        Mean Intensity values (a.u.)
    error : array
        Error values (currently not used in plotting)
    """
    # Read Excel file, skipping the header row
    df = pd.read_excel(filename)
    
    # Extract columns (columns 2 and 3 are power and intensity)
    # Column index 1 is X(potenza, mW), column index 2 is Y(intensity), column index 3 is error
    power = df.iloc[:, 1].values
    intensity = df.iloc[:, 2].values
    error = df.iloc[:, 3].values
    
    return power, intensity, error

def power_law(P, alpha, n):
    return alpha * P**n


def print_data_table(power, intensity, error):
    """Print a formatted table of all data points with significance flags."""
    mask = intensity > error
    print(f"{'#':>4}  {'P (mW)':>10}  {'I (a.u.)':>12}  {'err (a.u.)':>12}  {'included':>9}")
    print("-" * 58)
    for i, (p, I, e, ok) in enumerate(zip(power, intensity, error, mask)):
        flag = "YES" if ok else "NO  (I <= err)"
        print(f"{i+1:>4}  {p:>10.2f}  {I:>12.4f}  {e:>12.4f}  {flag}")
    print("-" * 58)
    n_excl = int((~mask).sum())
    print(f"Excluded: {n_excl} point(s) with I ≤ err  |  Kept: {int(mask.sum())} point(s)\n")

def plot_data(power, intensity, error, output_filename=None):
    """
    Plot Power vs Mean Intensity.
    
    Parameters:
    -----------
    power : array
        Power values (mW)
    intensity : array
        Mean Intensity values (a.u.)
    error : array
        Error values
    output_filename : str, optional
        If provided, save the plot to this filename
    """
    # --- skip first 3 points (not physically meaningful) ---
    power = power[3:]
    intensity = intensity[3:]
    error = error[3:]

    # --- filter: keep only statistically significant points ---
    mask = intensity > error

    plt.figure(figsize=(10, 6))
    plt.errorbar(power[mask], intensity[mask], yerr=error[mask],
                 fmt='o-', markersize=6, linewidth=1.5, capsize=4, label='Data')
    plt.xlabel('Power (mW)', fontsize=12)
    plt.ylabel('Mean Intensity (a.u.)', fontsize=12)
    #plt.title('Mean Intensity vs Power', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    if output_filename:
        save_dir = DEFAULT_OUTPUT_DIR
        os.makedirs(save_dir, exist_ok=True)
        output_path = os.path.join(save_dir, output_filename)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_path}")
    
    #plt.show()

def plot_data_with_fit(power, intensity, error, output_filename=None):
    """
    Plot Power vs Mean Intensity with power-law fit I = alpha * P^n,
    error bars, and reduced chi-square.
    Only data points with intensity > error are used for the fit.
    """
    # --- skip first 3 points (not physically meaningful) ---
    power = power[3:]
    intensity = intensity[3:]
    error = error[3:]

    # --- filter: keep only statistically significant points ---
    mask = intensity > error
    P_fit_data  = power[mask]
    I_fit_data  = intensity[mask]
    err_fit_data = error[mask]

    # Fit the power-law model on the valid subset, weighted by errors.
    # absolute_sigma=True: pcov is in physical units, not rescaled by chi2.
    # This gives correct parameter uncertainties when the supplied sigma values
    # are reliable estimates of the true measurement uncertainty.
    popt, pcov = curve_fit(power_law, P_fit_data, I_fit_data, p0=[1.0, 2.0],
                           sigma=err_fit_data, absolute_sigma=True)
    alpha, n = popt
    alpha_err, n_err = np.sqrt(np.diag(pcov))

    # Reduced chi-square on the fit data only
    residuals = I_fit_data - power_law(P_fit_data, alpha, n)
    chi2 = np.sum((residuals / err_fit_data) ** 2)
    dof = len(P_fit_data) - len(popt)      # degrees of freedom
    chi2_r = chi2 / dof

    # R² (coefficient of determination) — does not depend on error estimates
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((I_fit_data - I_fit_data.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot

    print(f"Fit results (on {int(mask.sum())} points): I = α · P^n")
    print(f"  α = {alpha:.4e} ± {alpha_err:.4e}")
    print(f"  n = {n:.4f} ± {n_err:.4f}")
    print(f"  Reduced χ² = {chi2_r:.4f}  (dof = {dof})")
    print(f"  R²          = {r2:.4f}")
    #if chi2_r < 0.1:
        #print(f"  Note: χ²_r << 1 suggests the error bars are overestimated.")
    #elif chi2_r > 10:
        #print(f"  Note: χ²_r >> 1 suggests the model is a poor fit or errors are underestimated.")

    # --- smooth curve over the full power range ---
    P_smooth = np.linspace(power.min(), power.max(), 300)
    I_smooth = power_law(P_smooth, alpha, n)

    fig, ax = plt.subplots(figsize=(10, 6))

    # valid points
    ax.errorbar(P_fit_data, I_fit_data, yerr=err_fit_data,
                fmt='o', markersize=6, capsize=4, label='Data')

    ax.plot(P_smooth, I_smooth, '-', linewidth=1.5,
            label=(
                f'Fit: $I = \\alpha \\cdot P^n$\n'
                f'$\\alpha = {alpha:.2e} \\pm {alpha_err:.2e}$\n'
                f'$n = {n:.3f} \\pm {n_err:.3f}$\n'
                f'$\\chi^2_r = {chi2_r:.4f}$'
            ))

    ax.set_xlabel('Power (mW)', fontsize=12)
    ax.set_ylabel('Mean Intensity (a.u.)', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=11)
    fig.tight_layout()

    if output_filename:
        save_dir = DEFAULT_OUTPUT_DIR
        os.makedirs(save_dir, exist_ok=True)
        output_path = os.path.join(save_dir, output_filename)
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_path}")

    #plt.show()

if __name__ == "__main__":
    # Check if filename was provided as command-line argument
    if len(sys.argv) < 2:
        print("Usage: python plot_power_intensity.py <data_file.xlsx>")
        print("Example: python plot_power_intensity.py my_data.xlsx")
        sys.exit(1)
    
    # Get the filename from command-line argument
    data_file = sys.argv[1]
    
    try:
        # Read the data
        power, intensity, error = read_data(data_file)
        
        # Print full data table with significance flags
        print(f"\nData read from '{data_file}':")
        print_data_table(power, intensity, error)
        
        # Plot the data
        plot_data(power, intensity, error, output_filename="power_intensity_plot.svg")
        
        # Plot with power-law fit
        plot_data_with_fit(power, intensity, error, output_filename="power_intensity_fit.svg")
        
    except FileNotFoundError:
        print(f"Error: File '{data_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
