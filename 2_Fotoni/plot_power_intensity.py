import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

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

def plot_data(power, intensity, output_filename=None):
    """
    Plot Power vs Mean Intensity.
    
    Parameters:
    -----------
    power : array
        Power values (mW)
    intensity : array
        Mean Intensity values (a.u.)
    output_filename : str, optional
        If provided, save the plot to this filename
    """
    plt.figure(figsize=(10, 6))
    plt.plot(power, intensity, 'o-', markersize=6, linewidth=1.5, label='Data')
    plt.xlabel('Power (mW)', fontsize=12)
    plt.ylabel('Mean Intensity (a.u.)', fontsize=12)
    plt.title('Mean Intensity vs Power', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    if output_filename:
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_filename}")
    
    plt.show()

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
        
        # Print summary
        print(f"Read {len(power)} data points from '{data_file}'")
        print(f"Power range: {power.min():.2f} - {power.max():.2f} mW")
        print(f"Intensity range: {intensity.min():.2f} - {intensity.max():.2f} a.u.")
        
        # Plot the data
        plot_data(power, intensity, output_filename="power_intensity_plot.svg")
        
    except FileNotFoundError:
        print(f"Error: File '{data_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
