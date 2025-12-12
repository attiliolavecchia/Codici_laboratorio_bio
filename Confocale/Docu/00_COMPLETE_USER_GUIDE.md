# Complete User Guide for MSD Analysis Toolkit

## Table of Contents
1. [Quick Start](#quick-start)
2. [Prerequisites & Setup](#prerequisites--setup)
3. [Script Reference Guide](#script-reference-guide)
4. [Data Requirements](#data-requirements)
5. [Common Workflows](#common-workflows)
6. [Troubleshooting](#troubleshooting)
7. [Output Interpretation](#output-interpretation)

---

## Quick Start

### For Complete Beginners

1. **Install Python dependencies**:
   ```powershell
   pip install pandas numpy scipy matplotlib
   ```

2. **Prepare your data**: Place CSV files with trajectory data in the appropriate folders:
   - `Data/17_11_anomalous/` for anomalous diffusion data
   - `Data/31_10_no_anomalous/` for normal diffusion data

3. **Run automated analysis**:
   ```powershell
   # For normal diffusion data
   python batch_analysis_no_anomalous.py
   
   # For anomalous diffusion data  
   python batch_analysis_anomalous.py
   ```

4. **Check results**: Look in the generated folders (`eamsd_plots/`, `tamsd_plots/`, `linear_fits/`, `nonlinear_fits/`) and the documentation file in `Docu/`

---

## Prerequisites & Setup

### Required Python Packages
```powershell
pip install pandas numpy scipy matplotlib
```

### Data Format Requirements
Your CSV files must contain these columns (case-insensitive, various naming supported):
- **Track ID**: `Track ID`, `track_id`, `trackid`, `ID Track`
- **X Position**: `X position`, `x`, `Position X`, `Pos X`, `position_x` 
- **Y Position**: `Y position`, `y`, `Position Y`, `Pos Y`, `position_y`
- **Time**: `Time`, `t`, `Position T`, `Pos T`, `position_t`

### Directory Structure
```
Confocale/
├── Data/
│   ├── 17_11_anomalous/     # For anomalous diffusion data
│   └── 31_10_no_anomalous/  # For normal diffusion data
├── Docu/                    # Documentation (auto-generated logs)
├── eamsd_plots/            # Output: Ensemble-averaged MSD plots
├── tamsd_plots/            # Output: Time-averaged MSD plots  
├── linear_fits/            # Output: Linear/anomalous fit results
├── nonlinear_fits/         # Output: Nonlinear/drift fit results
└── anomalous_analysis/     # Output: Anomalous data analysis results
```

---

## Script Reference Guide

### 1. Automated Batch Analysis Scripts

#### `batch_analysis_no_anomalous.py`
**Purpose**: Automated analysis of normal diffusion data

**Input**: All CSV files in `Data/31_10_no_anomalous/`

**What it does**:
1. TAMSD analysis on longest track for multiple lag fractions (0.1, 0.25, 0.5, 1.0)
2. EA-MSD analysis for multiple lag fractions (0.1, 0.25, 0.5, 1.0)
3. Linear fitting (normal diffusion) for multiple fractions (0.1, 0.25, 0.5)
4. Nonlinear fitting (drift detection) for multiple fractions (0.1, 0.25, 0.5)
5. Generates summary report in `Docu/05_BATCH_ANALYSIS_LOG.md`

**Usage**:
```powershell
python batch_analysis_no_anomalous.py
```

**Output folders**:
- `tamsd_plots/` - Time-averaged MSD plots
- `eamsd_plots/` - Ensemble-averaged MSD plots
- `linear_fits/` - Linear fit results and plots
- `nonlinear_fits/` - Nonlinear fit results and plots

#### `batch_analysis_anomalous.py`
**Purpose**: Automated analysis of anomalous diffusion data

**Input**: All CSV files in `Data/17_11_anomalous/`

**What it does**:
1. TAMSD analysis on longest track for multiple lag fractions (0.1, 0.25, 0.5, 1.0)
2. EA-MSD analysis for multiple lag fractions (0.1, 0.25, 0.5, 1.0)
3. Anomalous fitting (no drift) for multiple fractions (0.1, 0.25, 0.5)
4. Anomalous + drift fitting for multiple fractions (0.1, 0.25, 0.5)
5. Generates summary report in `Docu/06_BATCH_ANALYSIS_ANOMALOUS_LOG.md`

**Usage**:
```powershell
python batch_analysis_anomalous.py
```

**Output folders**:
- `anomalous_analysis/tamsd_plot/` - Time-averaged MSD plots
- `anomalous_analysis/eamsd_plot/` - Ensemble-averaged MSD plots
- `anomalous_analysis/linear_fits/` - Anomalous fit results (no drift)
- `anomalous_analysis/nonlinear_fits/` - Anomalous + drift fit results

---

### 2. Individual Analysis Scripts

#### `eamsd_plot.py`
**Purpose**: Generate ensemble-averaged MSD plots

**Usage**:
```powershell
# Basic usage
python eamsd_plot.py "Data/your_file.csv"

# With custom lag fraction (use only 50% of data)
python eamsd_plot.py "Data/your_file.csv" --max-lag-fraction 0.5

# Save to specific file
python eamsd_plot.py "Data/your_file.csv" --output "my_eamsd_plot.svg"
```

**Parameters**:
- `csv` (required): Path to CSV file
- `--max-lag-fraction`: Maximum lag fraction to use (0.1-1.0, default: auto)
- `--output`: Output filename (default: auto-generated)

#### `tamsd_plot.py`
**Purpose**: Generate time-averaged MSD plots for single trajectory

**Usage**:
```powershell
# Analyze specific track ID
python tamsd_plot.py "Data/your_file.csv" --track-id 5

# With custom lag fraction
python tamsd_plot.py "Data/your_file.csv" --track-id 5 --max-lag-fraction 0.3

# Save to specific file
python tamsd_plot.py "Data/your_file.csv" --track-id 5 --output "tamsd_track5.svg"
```

**Parameters**:
- `csv` (required): Path to CSV file
- `--track-id` (required): ID of the track to analyze
- `--max-lag-fraction`: Maximum lag fraction to use (default: auto)
- `--output`: Output filename (default: auto-generated)

#### `msd_fitting.py`
**Purpose**: Linear MSD fitting (normal diffusion: MSD = 4D·τ)

**Usage**:
```powershell
# Basic fitting
python msd_fitting.py "Data/your_file.csv"

# Use specific fraction of data for fitting
python msd_fitting.py "Data/your_file.csv" --manual-fraction 0.3

# Specify output directory
python msd_fitting.py "Data/your_file.csv" --output-dir "my_fits"
```

**Parameters**:
- `csv` (required): Path to CSV file
- `--manual-fraction`: Fraction of MSD data to use for fitting (0.1-1.0)
- `--output-dir`: Directory to save results (default: current directory)

**Output**:
- Diffusion coefficient D [μm²/s]
- R² value
- Plot with fit line
- Text file with detailed results

#### `msd_fitting_anomalous.py`
**Purpose**: Anomalous diffusion fitting (MSD = 4D_α·τ^α)

**Usage**:
```powershell
# Basic anomalous fitting
python msd_fitting_anomalous.py "Data/your_file.csv"

# Use specific fraction for fitting
python msd_fitting_anomalous.py "Data/your_file.csv" --manual-fraction 0.25

# Specify output directory
python msd_fitting_anomalous.py "Data/your_file.csv" --output-dir "anomalous_fits"
```

**Parameters**:
- `csv` (required): Path to CSV file
- `--manual-fraction`: Fraction of MSD data to use for fitting (0.1-1.0)
- `--output-dir`: Directory to save results (default: current directory)

**Output**:
- Generalized diffusion coefficient D_α [μm²/s^α]
- Anomalous exponent α (α<1: subdiffusion, α=1: normal, α>1: superdiffusion)
- R² value
- Plot with fit curve
- Text file with detailed results

#### `msd_fitting_anomalous_drift.py`
**Purpose**: Anomalous diffusion with drift fitting (MSD = 4D_α·τ^α + v²·τ²)

**Usage**:
```powershell
# Basic anomalous + drift fitting
python msd_fitting_anomalous_drift.py "Data/your_file.csv"

# Use specific fraction for fitting
python msd_fitting_anomalous_drift.py "Data/your_file.csv" --manual-fraction 0.4

# Specify output directory
python msd_fitting_anomalous_drift.py "Data/your_file.csv" --output-dir "drift_fits"
```

**Parameters**:
- `csv` (required): Path to CSV file
- `--manual-fraction`: Fraction of MSD data to use for fitting (0.1-1.0)
- `--output-dir`: Directory to save results (default: current directory)

**Output**:
- Generalized diffusion coefficient D_α [μm²/s^α]
- Anomalous exponent α
- Drift velocity v [μm/s]
- R² value
- Plot with fit curve
- Text file with detailed results

#### `compare_msd.py`
**Purpose**: Interactive comparison of multiple MSD datasets

**Usage**:
```powershell
python compare_msd.py
```

**Interactive features**:
1. **EA-MSD comparison**: Select multiple files with different timesteps
2. **TA-MSD comparison**: Select one file and track to compare different truncation lengths
3. **Lag fraction selection**: Choose what fraction of data to include
4. **Automatic output**: Saves comparison plots with parameter-encoded filenames

**Output**:
- Comparison plots in `eamsd_plots/` or `tamsd_plots/`
- Overlay of multiple datasets with legends
- Statistical information in plots

---

### 3. Utility Scripts

#### `test_msd_simulation.py`
**Purpose**: Validate the analysis pipeline using simulated data

**Usage**:
```powershell
python test_msd_simulation.py
```

**What it does**:
- Simulates Brownian motion with known parameters
- Tests linear and anomalous fitting algorithms
- Validates that known input parameters are recovered
- Reports accuracy metrics

#### `plot_msd_fit.py`
**Purpose**: Advanced plotting with multiple fit models

**Usage**:
```powershell
# Linear fit plotting
python plot_msd_fit.py "Data/your_file.csv" --fit-type linear

# Nonlinear fit plotting  
python plot_msd_fit.py "Data/your_file.csv" --fit-type nonlinear

# Custom tau max for fitting
python plot_msd_fit.py "Data/your_file.csv" --tau-max 15.0
```

---

## Data Requirements

### CSV File Format
Your CSV files should have columns for:
- **Track/Trajectory ID**: Unique identifier for each particle trajectory
- **X Position**: X coordinate (preferably in micrometers)
- **Y Position**: Y coordinate (preferably in micrometers)  
- **Time**: Time point (preferably in seconds)

### Supported Column Names
The software automatically recognizes these column names (case-insensitive):

| Data Type | Accepted Names |
|-----------|----------------|
| Track ID | `Track ID`, `track_id`, `trackid`, `ID Track` |
| X Position | `X position`, `x`, `Position X`, `Pos X`, `position_x` |
| Y Position | `Y position`, `y`, `Position Y`, `Pos Y`, `position_y` |
| Time | `Time`, `t`, `Position T`, `Pos T`, `position_t` |

### Data Quality Requirements
- **Minimum trajectory length**: At least 10 time points per track
- **Time uniformity**: Approximately uniform time intervals (small variations are handled)
- **Position units**: Micrometers preferred (affects diffusion coefficient units)
- **Time units**: Seconds preferred (affects diffusion coefficient units)

---

## Common Workflows

### Workflow 1: Quick Analysis of New Data

1. **Prepare data**: Put CSV files in `Data/31_10_no_anomalous/` (normal) or `Data/17_11_anomalous/` (anomalous)

2. **Run automated analysis**:
   ```powershell
   python batch_analysis_no_anomalous.py    # or batch_analysis_anomalous.py
   ```

3. **Check results**: 
   - View summary in `Docu/05_BATCH_ANALYSIS_LOG.md` (or `06_BATCH_ANALYSIS_ANOMALOUS_LOG.md`)
   - Browse plots in output folders
   - Check fitted parameters in the summary table

### Workflow 2: Detailed Analysis of Single File

1. **Generate MSD plots**:
   ```powershell
   python eamsd_plot.py "Data/your_file.csv" --max-lag-fraction 0.5
   python tamsd_plot.py "Data/your_file.csv" --track-id 1 --max-lag-fraction 0.5
   ```

2. **Fit different models**:
   ```powershell
   # Try linear fitting first
   python msd_fitting.py "Data/your_file.csv" --manual-fraction 0.3
   
   # If α ≠ 1, try anomalous fitting
   python msd_fitting_anomalous.py "Data/your_file.csv" --manual-fraction 0.3
   
   # If drift suspected, try anomalous + drift
   python msd_fitting_anomalous_drift.py "Data/your_file.csv" --manual-fraction 0.3
   ```

3. **Compare results**: Look at R² values and physical reasonableness of parameters

### Workflow 3: Comparing Multiple Experiments

1. **Interactive comparison**:
   ```powershell
   python compare_msd.py
   ```

2. **Follow prompts**:
   - Select comparison type (EA-MSD or TA-MSD)
   - Choose files to compare
   - Select lag fraction
   - Review generated plots

### Workflow 4: Parameter Optimization

1. **Test different lag fractions**: Start with batch analysis (tests 0.1, 0.25, 0.5)

2. **Examine R² values** in the summary file to find optimal fraction

3. **Re-run with optimal fraction**:
   ```powershell
   python msd_fitting_anomalous.py "Data/your_file.csv" --manual-fraction 0.25
   ```

---

## Troubleshooting

### Common Errors

#### "Could not find required columns"
**Problem**: CSV headers don't match expected names
**Solution**: 
- Check your CSV file has Track ID, X position, Y position, Time columns
- Verify column names match supported formats (see [Data Requirements](#data-requirements))
- Try opening CSV in text editor to check header formatting

#### "No tracks found" or "Empty trajectory"
**Problem**: Data reading failed
**Solution**:
- Check CSV file is properly formatted
- Verify data types (numbers in position/time columns)
- Ensure Track IDs are consistent

#### "Fitting failed" or "Poor R² values"
**Problem**: Model doesn't fit data well
**Solution**:
- Try different lag fractions (--manual-fraction)
- Check if data shows anomalous behavior (use anomalous fitting)
- Verify data quality (sufficient points, reasonable noise level)

#### "FileNotFoundError"
**Problem**: File or directory not found
**Solution**:
- Use absolute paths: `python script.py "C:\full\path\to\file.csv"`
- Check file exists and spelling is correct
- Verify working directory is correct

### Performance Issues

#### Slow processing
- **Large files**: Consider using smaller lag fractions
- **Many files**: Batch scripts process sequentially - be patient
- **Memory issues**: Close other applications, process files individually

#### Poor fits
- **Try different models**: Linear → Anomalous → Anomalous+Drift
- **Adjust lag fractions**: Start with 0.1-0.3 for initial fitting range
- **Check data quality**: Remove short tracks, verify time intervals

---

## Output Interpretation

### Diffusion Coefficient (D)
- **Units**: μm²/s (if positions in μm, time in s)
- **Typical values for 240nm particles in glycerol**: 0.001-0.01 μm²/s
- **Physical meaning**: Rate of diffusive spreading

### Anomalous Exponent (α)
- **α < 1**: Subdiffusion (confined, hindered motion)
- **α = 1**: Normal Brownian diffusion
- **α > 1**: Superdiffusion (active transport, ballistic motion)
- **α ≈ 2**: Pure ballistic motion

### Drift Velocity (v)
- **Units**: μm/s (if positions in μm, time in s)
- **Physical meaning**: Average directed motion speed
- **Typical values**: 0.001-1 μm/s depending on system

### R² (Coefficient of Determination)
- **Range**: 0-1 (higher is better)
- **Good fits**: R² > 0.95
- **Acceptable fits**: R² > 0.90
- **Poor fits**: R² < 0.85 (consider different model)

### File Naming Convention
Output files use parameter-encoded names for organization:
- `filename_track<ID>_lag<fraction>.svg` - Plots
- `filename_frac<fraction>_fit_results.txt` - Fit parameters
- `compare_eamsd_<timesteps>_f<frac>_<timestamp>.svg` - Comparison plots

This naming prevents overwriting and makes it easy to identify analysis parameters.