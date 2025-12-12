# Script Reference Documentation

## Overview
This document provides detailed technical information about each script in the MSD Analysis Toolkit, including all command-line parameters, input/output specifications, and advanced usage examples.

---

## Core Analysis Scripts

### `data_reader.py`
**Module Type**: Utility library (not executable)
**Purpose**: Robust CSV reading and trajectory data parsing

**Key Functions**:
```python
from data_reader import read_trajectories_from_csv, estimate_global_time_step

# Read trajectories
trajectories = read_trajectories_from_csv("file.csv")  
# Returns: Dict[track_id, Trajectory]

# Estimate time step
dt = estimate_global_time_step(trajectories)
# Returns: float (median dt across all tracks)
```

**Trajectory Object**:
```python
@dataclass
class Trajectory:
    track_id: str
    time: np.ndarray     # Time points [s]
    x: np.ndarray        # X positions [μm]
    y: np.ndarray        # Y positions [μm]  
    dt: float            # Time step [s]
```

**Supported Column Headers**:
- Track ID: `Track ID`, `track_id`, `trackid`, `ID Track`
- X Position: `X position`, `x`, `Position X`, `Pos X`, `position_x`
- Y Position: `Y position`, `y`, `Position Y`, `Pos Y`, `position_y`
- Time: `Time`, `t`, `Position T`, `Pos T`, `position_t`

---

### `msd_analyzer.py`
**Module Type**: Utility library (not executable)
**Purpose**: MSD calculation algorithms

**Key Functions**:

#### Ensemble-Averaged MSD (initial displacement only)
```python
from msd_analyzer import calculate_ensemble_msd

result = calculate_ensemble_msd(
    trajectories,           # Dict[track_id, Trajectory]
    max_lag_steps=None,     # int or None (auto-determine)
    global_dt=None          # float or None (auto-estimate)
)
# Returns: MSDResult(tau, msd, n_trajectories)
```

#### Time-Averaged MSD (single trajectory)
```python
from msd_analyzer import calculate_time_averaged_msd_per_track

tau, msd = calculate_time_averaged_msd_per_track(
    trajectory,             # Single Trajectory object
    max_lag_steps=None,     # int or None
    dt_override=None        # float or None
)
# Returns: (tau_array, msd_array)
```

---

## Plotting Scripts

### `eamsd_plot.py`
**Purpose**: Generate ensemble-averaged MSD plots

```bash
python eamsd_plot.py <csv_file> [OPTIONS]
```

**Required Arguments**:
- `csv`: Path to trajectory CSV file

**Optional Arguments**:
- `--max-lag-fraction FLOAT`: Maximum lag fraction (0.0-1.0, default: auto-determined)
- `--output PATH`: Output file path (default: auto-generated .svg)

**Examples**:
```powershell
# Basic usage
python eamsd_plot.py "Data/experiment1.csv"

# Use only first 30% of lag times
python eamsd_plot.py "Data/experiment1.csv" --max-lag-fraction 0.3

# Custom output file
python eamsd_plot.py "Data/experiment1.csv" --output "my_eamsd.svg"

# Full path example
python eamsd_plot.py "C:\Data\trajectories\sample.csv" --max-lag-fraction 0.5 --output "results\eamsd_analysis.svg"
```

**Output**:
- SVG plot file with linear axes
- X-axis: Time Lag (τ) [seconds]
- Y-axis: Ensemble-Averaged MSD
- Grid, legend, and publication-quality formatting

---

### `tamsd_plot.py`
**Purpose**: Generate time-averaged MSD plots for single trajectory

```bash
python tamsd_plot.py <csv_file> --track-id <ID> [OPTIONS]
```

**Required Arguments**:
- `csv`: Path to trajectory CSV file
- `--track-id INT`: ID of the trajectory to analyze

**Optional Arguments**:
- `--max-lag-fraction FLOAT`: Maximum lag fraction (0.0-1.0, default: auto-determined)
- `--output PATH`: Output file path (default: auto-generated .svg)

**Examples**:
```powershell
# Basic usage
python tamsd_plot.py "Data/experiment1.csv" --track-id 5

# Limit to first 25% of lag times
python tamsd_plot.py "Data/experiment1.csv" --track-id 5 --max-lag-fraction 0.25

# Custom output
python tamsd_plot.py "Data/experiment1.csv" --track-id 3 --output "track3_tamsd.svg"
```

**How to find Track IDs**:
```python
# Method 1: Examine CSV file directly
import pandas as pd
df = pd.read_csv("your_file.csv")
track_ids = df['Track ID'].unique()  # Adjust column name as needed
print("Available Track IDs:", track_ids)

# Method 2: Use data_reader
from data_reader import read_trajectories_from_csv
trajectories = read_trajectories_from_csv("your_file.csv")
print("Available Track IDs:", list(trajectories.keys()))
```

**Output**:
- SVG plot file with linear axes
- X-axis: Time Lag (τ) [seconds]
- Y-axis: Time-Averaged MSD
- Title includes track ID for identification

---

## Fitting Scripts

### `msd_fitting.py`
**Purpose**: Linear MSD fitting (MSD = 4D·τ)

```bash
python msd_fitting.py <csv_file> [OPTIONS]
```

**Required Arguments**:
- `csv`: Path to trajectory CSV file

**Optional Arguments**:
- `--manual-fraction FLOAT`: Fraction of MSD data to use for fitting (0.1-1.0)
- `--output-dir PATH`: Directory to save results (default: current directory)

**Examples**:
```powershell
# Automatic fitting range selection
python msd_fitting.py "Data/experiment1.csv"

# Use first 30% of MSD data for fitting
python msd_fitting.py "Data/experiment1.csv" --manual-fraction 0.3

# Save to specific directory
python msd_fitting.py "Data/experiment1.csv" --output-dir "linear_fits"

# Combined options
python msd_fitting.py "Data/experiment1.csv" --manual-fraction 0.25 --output-dir "results\linear"
```

**Fitting Algorithm**:
- Tests multiple lag fractions (10%-90% in 5% steps)
- Selects optimal range based on R² and residual analysis
- Manual fraction overrides automatic selection
- Uses scipy.optimize.curve_fit with Levenberg-Marquardt

**Output Files**:
- `<filename>_frac<fraction>_linear_fit.svg`: Plot with fit line
- `<filename>_frac<fraction>_fit_results.txt`: Detailed results

**Output Format (text file)**:
```
Linear MSD Fitting Results
=========================
File: experiment1.csv
Model: MSD = 4D·τ
Fitting fraction: 0.30

Parameters:
-----------
Diffusion Coefficient (D): (6.78 ± 0.12) × 10⁻³ μm²/s
R²: 0.9856

Fitting Details:
---------------
Lag range: 0.50 - 7.50 seconds
Data points used: 25
```

---

### `msd_fitting_anomalous.py`
**Purpose**: Anomalous diffusion fitting (MSD = 4D_α·τ^α)

```bash
python msd_fitting_anomalous.py <csv_file> [OPTIONS]
```

**Required Arguments**:
- `csv`: Path to trajectory CSV file

**Optional Arguments**:
- `--manual-fraction FLOAT`: Fraction of MSD data to use for fitting (0.1-1.0)
- `--output-dir PATH`: Directory to save results (default: current directory)

**Examples**:
```powershell
# Automatic fitting
python msd_fitting_anomalous.py "Data/anomalous_data.csv"

# Use first 40% of data
python msd_fitting_anomalous.py "Data/anomalous_data.csv" --manual-fraction 0.4

# Save to specific directory
python msd_fitting_anomalous.py "Data/anomalous_data.csv" --output-dir "anomalous_fits"
```

**Fitting Algorithm**:
- Power-law model: MSD(τ) = 4D_α·τ^α
- Tests multiple lag fractions for optimal R²
- Handles both subdiffusion (α < 1) and superdiffusion (α > 1)
- Robust parameter estimation with bounds

**Physical Interpretation**:
- **α < 1**: Subdiffusion (confined motion, crowded environments)
- **α = 1**: Normal Brownian diffusion
- **α > 1**: Superdiffusion (active transport, ballistic components)

**Output Files**:
- `<filename>_frac<fraction>_anomalous_fit.svg`: Plot with power-law fit
- `<filename>_frac<fraction>_anomalous_results.txt`: Detailed results

**Output Format (text file)**:
```
Anomalous Diffusion Fitting Results
==================================
File: anomalous_data.csv
Model: MSD = 4D_α·τ^α
Fitting fraction: 0.40

Parameters:
-----------
Diffusion Coefficient (D_α): (2.34 ± 0.08) × 10⁻³ μm²/s^α
Anomalous exponent (α): (0.78 ± 0.03)
R²: 0.9923

Interpretation:
--------------
α < 1: Subdiffusion detected
Physical mechanism: Confined or hindered motion

Fitting Details:
---------------
Lag range: 0.50 - 12.0 seconds
Data points used: 48
```

---

### `msd_fitting_anomalous_drift.py`
**Purpose**: Anomalous diffusion with drift (MSD = 4D_α·τ^α + v²·τ²)

```bash
python msd_fitting_anomalous_drift.py <csv_file> [OPTIONS]
```

**Required Arguments**:
- `csv`: Path to trajectory CSV file

**Optional Arguments**:
- `--manual-fraction FLOAT`: Fraction of MSD data to use for fitting (0.1-1.0)
- `--output-dir PATH`: Directory to save results (default: current directory)

**Examples**:
```powershell
# Automatic fitting
python msd_fitting_anomalous_drift.py "Data/drift_data.csv"

# Use first 35% of data  
python msd_fitting_anomalous_drift.py "Data/drift_data.csv" --manual-fraction 0.35

# Save results to specific directory
python msd_fitting_anomalous_drift.py "Data/drift_data.csv" --output-dir "drift_analysis"
```

**Model**: MSD(τ) = 4D_α·τ^α + v²·τ²
- **Diffusive component**: 4D_α·τ^α (anomalous diffusion)
- **Ballistic component**: v²·τ² (directed motion)

**Fitting Algorithm**:
- Three-parameter nonlinear fit (D_α, α, v)
- Automatic velocity bounds estimation from trajectory data
- Handles complex transport mechanisms
- Robust optimization with physical constraints

**Output Files**:
- `<filename>_frac<fraction>_anomalous_drift_fit.svg`: Plot with combined fit
- `<filename>_frac<fraction>_drift_results.txt`: Detailed results

**Output Format (text file)**:
```
Anomalous Diffusion + Drift Fitting Results
==========================================
File: drift_data.csv
Model: MSD = 4D_α·τ^α + v²·τ²
Fitting fraction: 0.35

Parameters:
-----------
Diffusion Coefficient (D_α): (1.89 ± 0.05) × 10⁻³ μm²/s^α
Anomalous exponent (α): (0.85 ± 0.04)
Drift Velocity (v): (0.12 ± 0.01) μm/s
R²: 0.9945

Interpretation:
--------------
α < 1: Subdiffusion with directed motion
Significant drift detected (v > 0.05 μm/s)

Physical mechanism: Combined anomalous diffusion and active transport

Fitting Details:
---------------
Lag range: 0.50 - 14.0 seconds
Data points used: 56
```

---

## Comparison Scripts

### `compare_msd.py`
**Purpose**: Interactive comparison of multiple MSD datasets

```bash
python compare_msd.py
```

**No command-line arguments** - fully interactive interface

**Workflow**:
1. **Analysis type selection**:
   ```
   Choose analysis type:
   1. Compare EA-MSD across different experiments/timesteps
   2. Compare TA-MSD with different truncation lengths (single trajectory)
   ```

2. **File selection** (for EA-MSD comparison):
   - File browser interface
   - Select multiple CSV files
   - Automatic timestep extraction from filenames

3. **Parameters**:
   - Max lag fraction: What portion of lag times to include
   - Output filename customization

4. **Automatic processing**:
   - Computes MSD for each selected dataset
   - Creates overlay plots with legends
   - Saves with parameter-encoded filenames

**Output Files**:
- EA-MSD comparison: `eamsd_plots/compare_eamsd_<timesteps>_f<frac>_<timestamp>.svg`
- TA-MSD comparison: `tamsd_plots/compare_tamsd_ts<timestep>_track<ID>_L<lengths>_f<frac>_<timestamp>.svg`

**Example Output Filename**:
`compare_eamsd_40-50-60min_f0.5_20241121_143052.svg`
- Timesteps: 40, 50, 60 minutes
- Lag fraction: 0.5 (50% of data)
- Timestamp: 2024-11-21 14:30:52

---

## Batch Analysis Scripts

### `batch_analysis_no_anomalous.py`
**Purpose**: Automated analysis of normal diffusion datasets

```bash
python batch_analysis_no_anomalous.py
```

**No command-line arguments** - processes all CSV files in `Data/31_10_no_anomalous/`

**Processing Steps**:
1. **TAMSD Analysis** (longest track per file):
   - Lag fractions: [0.1, 0.25, 0.5, 1.0]
   - Output: `tamsd_plots/<filename>_track<ID>_lag<frac>.svg`

2. **EA-MSD Analysis**:
   - Lag fractions: [0.1, 0.25, 0.5, 1.0] 
   - Output: `eamsd_plots/<filename>_lag<frac>.svg`

3. **Linear Fitting**:
   - Fit fractions: [0.1, 0.25, 0.5]
   - Output: `linear_fits/` directory

4. **Nonlinear Fitting** (drift detection):
   - Fit fractions: [0.1, 0.25, 0.5]
   - Output: `nonlinear_fits/` directory

5. **Summary Report**: `Docu/05_BATCH_ANALYSIS_LOG.md`

**Configuration (editable in script)**:
```python
# Data location
DATA_DIR = Path("Data/31_10_no_anomalous")

# Test parameters
TAMSD_FRACTIONS = [0.1, 0.25, 0.5, 1.0]
EAMSD_FRACTIONS = [0.1, 0.25, 0.5, 1.0]
LINEAR_FIT_FRACTIONS = [0.1, 0.25, 0.5]
```

**Output Summary Table**:
| File | Model | Fraction | D (μm²/s) | v (μm/s) | R² |
|------|-------|----------|-----------|----------|-----|
| file1.csv | Linear | 0.1 | 0.0067 | N/A | 0.985 |
| file1.csv | Nonlinear | 0.1 | 0.0065 | 0.023 | 0.987 |

---

### `batch_analysis_anomalous.py`
**Purpose**: Automated analysis of anomalous diffusion datasets

```bash
python batch_analysis_anomalous.py
```

**No command-line arguments** - processes all CSV files in `Data/17_11_anomalous/`

**Processing Steps**:
1. **TAMSD Analysis** (longest track per file):
   - Lag fractions: [0.1, 0.25, 0.5, 1.0]
   - Output: `anomalous_analysis/tamsd_plot/<filename>_track<ID>_lag<frac>.svg`

2. **EA-MSD Analysis**:
   - Lag fractions: [0.1, 0.25, 0.5, 1.0]
   - Output: `anomalous_analysis/eamsd_plot/<filename>_lag<frac>.svg`

3. **Anomalous Fitting** (no drift):
   - Fit fractions: [0.1, 0.25, 0.5]
   - Output: `anomalous_analysis/linear_fits/` directory

4. **Anomalous + Drift Fitting**:
   - Fit fractions: [0.1, 0.25, 0.5]
   - Output: `anomalous_analysis/nonlinear_fits/` directory

5. **Summary Report**: `Docu/06_BATCH_ANALYSIS_ANOMALOUS_LOG.md`

**Output Summary Table**:
| File | Model | Fraction | D_alpha | alpha | v (μm/s) | R² |
|------|-------|----------|---------|-------|----------|-----|
| file1.csv | Anomalous | 0.25 | 0.0034 | 0.78 | N/A | 0.992 |
| file1.csv | Anomalous+Drift | 0.25 | 0.0031 | 0.81 | 0.087 | 0.995 |

---

## Utility Scripts

### `test_msd_simulation.py`
**Purpose**: Validation of analysis pipeline using simulated data

```bash
python test_msd_simulation.py
```

**Validation Tests**:
1. **Linear Fitting Test**:
   - Simulates Brownian motion with known D
   - Tests if fitting recovers input D within error bounds
   - Reports accuracy statistics

2. **Anomalous Fitting Test**:
   - Simulates anomalous diffusion with known D_α and α
   - Tests parameter recovery accuracy
   - Validates R² values

**Example Output**:
```
--- Testing Linear Fit (Brownian Motion) ---
Simulating 50 tracks with D=0.5 um^2/s...
Calculating Ensemble MSD...
Fitting with automatic range selection...

Results:
True D: 0.500 um^2/s
Fitted D: (0.498 ± 0.008) um^2/s
Error: 0.4%
R²: 0.9987
✓ PASS: D within 5% of true value

--- Testing Anomalous Fit ---
Simulating anomalous diffusion: D_α=0.3, α=0.75...
Results:
True D_α: 0.300, True α: 0.75
Fitted D_α: (0.302 ± 0.005), Fitted α: (0.751 ± 0.008)  
Errors: D_α: 0.7%, α: 0.1%
R²: 0.9995
✓ PASS: Both parameters within 5% of true values
```

### `plot_msd_fit.py`
**Purpose**: Advanced plotting with customizable fit visualization

```bash
python plot_msd_fit.py <csv_file> [OPTIONS]
```

**Required Arguments**:
- `csv`: Path to trajectory CSV file

**Optional Arguments**:
- `--fit-type TYPE`: Fitting model (`linear`, `nonlinear`, `anomalous`, `anomalous_drift`)
- `--tau-max FLOAT`: Maximum tau for fitting (seconds)
- `--output PATH`: Output plot file
- `--show-residuals`: Include residual plot (subplot)

**Examples**:
```powershell
# Linear fit with default settings
python plot_msd_fit.py "Data/sample.csv" --fit-type linear

# Anomalous fit with custom tau range
python plot_msd_fit.py "Data/sample.csv" --fit-type anomalous --tau-max 20.0

# Include residual analysis
python plot_msd_fit.py "Data/sample.csv" --fit-type linear --show-residuals

# Custom output file
python plot_msd_fit.py "Data/sample.csv" --fit-type anomalous_drift --output "custom_fit.svg"
```

**Advanced Features**:
- Multiple fit model support
- Residual analysis plots
- Custom fitting ranges
- Publication-quality formatting
- Automatic parameter annotation

---

## Error Handling & Debugging

### Common Command-Line Errors

#### Invalid File Path
```bash
# Error
python eamsd_plot.py "nonexistent_file.csv"
# Output: FileNotFoundError: [Errno 2] No such file or directory

# Solution: Use absolute paths or verify file exists
python eamsd_plot.py "C:\full\path\to\file.csv"
```

#### Missing Track ID
```bash
# Error  
python tamsd_plot.py "Data/file.csv" --track-id 999
# Output: Error: Track ID 999 not found in file

# Solution: Check available track IDs first
python -c "from data_reader import *; t=read_trajectories_from_csv('Data/file.csv'); print(list(t.keys()))"
```

#### Invalid Fraction Values
```bash
# Error
python msd_fitting.py "Data/file.csv" --manual-fraction 1.5
# Output: Error: fraction must be between 0.1 and 1.0

# Solution: Use valid range
python msd_fitting.py "Data/file.csv" --manual-fraction 0.5
```

### Verbose Mode for Debugging

Most scripts support verbose output for troubleshooting:

```python
# Add to beginning of script for debugging
import logging
logging.basicConfig(level=logging.DEBUG)
```

This will provide detailed information about:
- File reading progress
- Trajectory parsing
- MSD calculation steps
- Fitting algorithm iterations
- Output file creation

### Performance Optimization

#### For Large Files
```bash
# Use smaller lag fractions to reduce computation time
python eamsd_plot.py "large_file.csv" --max-lag-fraction 0.3

# For fitting, start with smaller fractions
python msd_fitting.py "large_file.csv" --manual-fraction 0.2
```

#### For Many Files
```bash
# Batch scripts process sequentially - consider parallel processing
# Run multiple terminals with different file subsets
```

---

## Advanced Configuration

### Modifying Batch Scripts

Edit configuration variables in batch scripts:

```python
# In batch_analysis_no_anomalous.py or batch_analysis_anomalous.py

# Change data directory
DATA_DIR = Path("Data/your_custom_directory")

# Modify test parameters
TAMSD_FRACTIONS = [0.1, 0.2, 0.3, 0.4, 0.5]  # More granular testing
EAMSD_FRACTIONS = [0.1, 0.25, 0.5, 0.75, 1.0]  # Include 75% fraction
FIT_FRACTIONS = [0.1, 0.15, 0.2, 0.25, 0.3]    # Fine-grained fitting tests

# Change output directories
TAMSD_DIR = Path("custom_tamsd_output")
LINEAR_DIR = Path("custom_linear_output")
```

### Custom Analysis Scripts

Template for creating custom analysis:

```python
#!/usr/bin/env python3
"""Custom analysis script template."""

from pathlib import Path
from data_reader import read_trajectories_from_csv
from msd_analyzer import calculate_ensemble_msd
from msd_fitting_anomalous import fit_msd_anomalous

def analyze_file(csv_path):
    """Analyze single CSV file with custom parameters."""
    
    # Load data
    trajectories = read_trajectories_from_csv(csv_path)
    print(f"Loaded {len(trajectories)} trajectories")
    
    # Calculate MSD
    msd_result = calculate_ensemble_msd(trajectories, max_lag_steps=50)
    
    # Fit anomalous model
    fit_result = fit_msd_anomalous(
        msd_result.tau, 
        msd_result.msd, 
        manual_fraction=0.4
    )
    
    # Print results
    print(f"D_α = {fit_result.D_alpha:.3e} μm²/s^α")
    print(f"α = {fit_result.alpha:.3f}")
    print(f"R² = {fit_result.R_squared:.4f}")
    
    return fit_result

# Usage
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python custom_analysis.py <csv_file>")
        sys.exit(1)
    
    result = analyze_file(sys.argv[1])
```

This comprehensive reference should enable users to fully utilize all scripts with complete understanding of parameters and options.