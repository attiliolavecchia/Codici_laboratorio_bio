# MSD Analysis Toolkit - Project Overview

## Summary

This toolkit provides comprehensive analysis of particle tracking data from confocal microscopy experiments. It calculates Mean Squared Displacement (MSD) using multiple approaches and extracts physical parameters like diffusion coefficients and drift velocities.

## Documentation Guide

**New to the project?** Start here:
- **[📖 Complete User Guide](00_COMPLETE_USER_GUIDE.md)** - Comprehensive documentation with step-by-step instructions
- **[🚀 Getting Started for Beginners](08_GETTING_STARTED_BEGINNERS.md)** - Step-by-step setup and first analysis
- **[📚 Script Reference](07_SCRIPT_REFERENCE.md)** - Detailed technical documentation for all scripts

**Existing users:**
- **[⚙️ Core Analysis](02_DATA_ANALYSIS_CORE.md)** - MSD calculation modules
- **[📊 Comparison Analysis](03_COMPARISON_ANALYSIS.md)** - Multi-experiment studies  
- **[📈 MSD Fitting](04_MSD_FITTING.md)** - Quantitative parameter extraction

## Key Capabilities

1. **Data Processing**: Robust CSV reading with automatic header recognition (TrackMate compatible)
2. **MSD Calculation**: Both ensemble-averaged and time-averaged MSD computation
3. **Visualization**: Publication-quality plots for different analysis types
4. **Comparison Tools**: Multi-experiment and multi-condition analysis
5. **Automated Analysis**: Batch processing with comprehensive reports
6. **Physical Parameter Extraction**: 
   - Linear fitting: Diffusion coefficient (D) from MSD = 4D·τ
   - Nonlinear fitting: D + drift velocity (v) from MSD = 4D·τ + v²·τ²
   - Anomalous diffusion fitting: Generalized D and exponent α from MSD = 4D_α·τ^α
   - Anomalous + drift fitting: D_α, α, and v from MSD = 4D_α·τ^α + v²·τ²

## Experimental Context

**Target System**: 240nm particles in 85% glycerol/water solution
- **Expected D**: 0.006-0.009 μm²/s (Stokes-Einstein)
- **Typical conditions**: High viscosity, slow diffusion
- **Data format**: TrackMate CSV exports with Track ID, X/Y positions, Time

## Module Organization

### Core Analysis (→ `02_DATA_ANALYSIS_CORE.md`)
- `data_reader.py` - CSV reading and trajectory parsing
- `msd_analyzer.py` - MSD computation algorithms  
- `eamsd_plot.py` - Ensemble-averaged MSD plotting
- `tamsd_plot.py` - Time-averaged MSD plotting

### Comparative Analysis (→ `03_COMPARISON_ANALYSIS.md`)
- `compare_msd.py` - Multi-experiment and multi-condition comparisons

### Quantitative Fitting (→ `04_MSD_FITTING.md`)
- `msd_fitting.py` - Linear model fitting (pure normal diffusion)
- `msd_fitting_nonlinear.py` - Nonlinear model fitting (normal diffusion + drift)
- `msd_fitting_anomalous.py` - Anomalous diffusion fitting (no drift)
- `msd_fitting_anomalous_drift.py` - Anomalous diffusion with drift fitting
- `plot_msd_fit.py` - Fitted model visualization

---

## Quick Start for New Users

### 1. Install Dependencies
```powershell
pip install pandas numpy scipy matplotlib
```

### 2. Prepare Your Data
- Place CSV files in `Data/31_10_no_anomalous/` (normal diffusion) or `Data/17_11_anomalous/` (anomalous diffusion)
- CSV must have columns: Track ID, X position, Y position, Time

### 3. Run Automated Analysis
```powershell
# For normal diffusion data
python batch_analysis_no_anomalous.py

# For anomalous diffusion data  
python batch_analysis_anomalous.py
```

### 4. Check Results
- View summary reports in `Docu/` folder
- Browse plots in output folders (`eamsd_plots/`, `tamsd_plots/`, etc.)

### 5. Individual Analysis (Optional)
```powershell
# Create MSD plot
python eamsd_plot.py "Data/your_file.csv"

# Extract diffusion coefficient
python msd_fitting.py "Data/your_file.csv"

# Compare multiple experiments
python compare_msd.py
```

**Need more help?** → [🚀 Getting Started Guide](08_GETTING_STARTED_BEGINNERS.md)

---

## Script Quick Reference

| Script | Purpose | Usage |
|--------|---------|-------|
| `batch_analysis_no_anomalous.py` | Complete automated analysis (normal diffusion) | `python batch_analysis_no_anomalous.py` |
| `batch_analysis_anomalous.py` | Complete automated analysis (anomalous diffusion) | `python batch_analysis_anomalous.py` |
| `eamsd_plot.py` | Ensemble-averaged MSD plot | `python eamsd_plot.py "file.csv"` |
| `tamsd_plot.py` | Time-averaged MSD plot | `python tamsd_plot.py "file.csv" --track-id 5` |
| `msd_fitting.py` | Linear fitting (normal diffusion) | `python msd_fitting.py "file.csv"` |
| `msd_fitting_anomalous.py` | Anomalous diffusion fitting | `python msd_fitting_anomalous.py "file.csv"` |
| `msd_fitting_anomalous_drift.py` | Anomalous + drift fitting | `python msd_fitting_anomalous_drift.py "file.csv"` |
| `compare_msd.py` | Interactive multi-experiment comparison | `python compare_msd.py` |
| `test_msd_simulation.py` | Validate analysis pipeline | `python test_msd_simulation.py` |

**Full details** → [📚 Script Reference](07_SCRIPT_REFERENCE.md)

---
pip install numpy pandas matplotlib scipy

# Verify installation
python -c "import numpy, pandas, matplotlib, scipy; print('All packages installed')"
```

### Basic Workflow

1. **Ensemble Analysis** (population-level behavior):
   ```bash
   python eamsd_plot.py Data/Experiment8_spots_40minstep.csv
   ```

2. **Extract Diffusion Coefficient** (linear model):
   ```bash
   python plot_msd_fit.py Data/Experiment8_spots_40minstep.csv
   ```

3. **Advanced Analysis** (normal diffusion + drift):
   ```bash
   python msd_fitting_nonlinear.py Data/Experiment8_spots_40minstep.csv --compare-linear
   python plot_msd_fit.py --nonlinear Data/Experiment8_spots_40minstep.csv
   ```

4. **Anomalous Diffusion Analysis**:
   ```bash
   # Test for anomalous behavior
   python msd_fitting_anomalous.py Data/240nm_glicerolo50_serie004_251minstep.csv --compare-normal
   
   # Add drift term if needed
   python msd_fitting_anomalous_drift.py Data/240nm_glicerolo50_serie004_251minstep.csv --compare-no-drift --plot-velocity
   ```

5. **Multi-Experiment Comparison**:
   ```bash
   python compare_msd.py  # Interactive selection
   ```

### Example Results

**Linear Fitting** (Normal Diffusion):
- D = (1.437 ± 0.002) × 10⁻¹ μm²/s
- R² = 0.9912

**Nonlinear Fitting** (Normal Diffusion + Drift):  
- D = (6.306 ± 0.003) × 10⁻² μm²/s
- v = (5.002 ± 0.008) × 10⁻² μm/s
- R² = 0.9954

**Anomalous Diffusion Fitting**:
- D_α = (4.414 ± 0.259) × 10⁻³ μm²/s^α
- α = 1.592 ± 0.022 (superdiffusion)
- R² = 0.9965

**Anomalous Diffusion + Drift Fitting**:
- D_α = (7.642 ± 1.67) × 10⁻³ μm²/s^α
- α = 0.980 ± 0.158 (normal diffusion)
- v = (6.238 ± 0.341) × 10⁻² μm/s
- R² = 0.9972

## File Structure

```
Confocale/
├── Data/                           # Experimental data (CSV files)
├── Docu/                          # Documentation (this folder)
├── Core Analysis Modules:
│   ├── data_reader.py
│   ├── msd_analyzer.py  
│   ├── eamsd_plot.py
│   └── tamsd_plot.py
├── Comparison Module:
│   └── compare_msd.py
├── Fitting Modules:
│   ├── msd_fitting.py
│   ├── msd_fitting_nonlinear.py
│   └── plot_msd_fit.py
└── Output Directories (auto-created):
    ├── eamsd_plots/               # Ensemble MSD plots
    ├── tamsd_plots/               # Time-averaged MSD plots
    ├── linear_fits/               # Linear fit results
    └── nonlinear_fits/            # Nonlinear fit results
```

## Data Requirements

### Input CSV Format
Required columns (case-insensitive, various formats accepted):
- **Track ID**: `track_id`, `TrackId`, `TRACK_ID`, etc.
- **X Position**: `X position`, `POSITION_X`, `x`, etc.  
- **Y Position**: `Y position`, `POSITION_Y`, `y`, etc.
- **Time**: `Time`, `T`, `POSITION_T`, etc.

### Units
- **Positions**: μm (micrometers)
- **Time**: seconds
- **Output MSD**: μm²
- **Diffusion coefficient**: μm²/s
- **Velocity**: μm/s

## Analysis Approaches

### 1. Ensemble-Averaged MSD (eaMSD)
- **Purpose**: Population-level diffusion behavior
- **Method**: Initial displacement only (compare to t=0)
- **Best for**: Clean diffusion coefficient extraction
- **Formula**: MSD(τ) = ⟨[r(τ) - r(0)]²⟩

### 2. Time-Averaged MSD (TAMSD)  
- **Purpose**: Individual trajectory analysis
- **Method**: Sliding window within single trajectory
- **Best for**: Heterogeneity studies, single-particle behavior
- **Formula**: TAMSD(τ) = ⟨[r(t+τ) - r(t)]²⟩_t

### 3. Linear Fitting
- **Model**: MSD(τ) = 2D·τ
- **Assumption**: Pure Brownian diffusion
- **Output**: Diffusion coefficient D

### 4. Nonlinear Fitting
- **Model**: MSD(τ) = 4D·τ + v²·τ²  
- **Captures**: Diffusion + directed motion/drift
- **Output**: D (diffusion) + v (drift velocity)

## Next Steps

1. Read `02_DATA_ANALYSIS_CORE.md` for basic MSD analysis
2. Read `03_COMPARISON_ANALYSIS.md` for multi-experiment studies
3. Read `04_MSD_FITTING.md` for quantitative parameter extraction

## Support

Each module includes:
- Comprehensive error handling
- Built-in help (`--help` flag)
- Example usage in documentation
- Troubleshooting guides

For specific technical details, consult the relevant documentation section.

---

## Module Details

### data_reader.py
- Reads the input CSV and builds, for each track, a Trajectory object with NumPy arrays: time, x, y.
- Recognizes TrackMate-style headers (uppercase/lowercase and common variants) and enforces correct data types.
- Sorts points by time, removes temporal duplicates, and estimates the time step Δt for each track as the median of time differences.
- Key functions:
	- `read_trajectories_from_csv(path)` → dict {track_id → Trajectory}
	- `estimate_global_time_step(trajectories)` → Global Δt (median of Δt per track)

When to use it: always, because all other scripts rely on this module to read data robustly.

#### How particles are recognized and grouped (Track ID)
- The reader searches for a column equivalent to "Track ID" in the CSV (also recognizes common variants like "track id", "track_id", "TrackId", etc.).
- Track ID column values are converted to strings and stripped of spaces; this avoids problems if some rows have numeric IDs while others have text IDs.
- All points with the same Track ID are grouped (group-by) to form a single trajectory.
- Within each group:
	- points are sorted by increasing time;
	- any time duplicates are removed (keeping the first one);
	- a Δt per track is estimated as the median of temporal differences.
- The output is a dictionary: {track_id → Trajectory}. Each trajectory contains time, x, y, dt, and n_points.


---

### 2) `msd_analyzer.py`
- Calculates 2D MSD in two ways:
	1. Ensemble-averaged MSD (eaMSD): for each lag n compares position at time n with initial position (no averaging over starting points), then averages across available tracks.
	2. Time-averaged MSD (TAMSD) for **one** single track: averages over all possible windows within the track.
- Units:
	- τ (time lag) in seconds, using global Δt (eaMSD) or track Δt (TAMSD).
	- MSD in squared position units (e.g., μm², if X and Y are in microns).
- Parameter "how much lag to use": now completely user's choice. If `--max-lag-fraction` is not specified, **the entire available interval** is used (up to N_max − 1).
- Execution as script: prints Δt, n_max, τ_max and first τ–MSD values (useful for quick check).

Main functions (also usable from other scripts):
- `calculate_ensemble_msd(trajectories, max_lag_fraction=None, global_dt=None)`
- `calculate_time_averaged_msd_per_track(track, max_lag_fraction=None, dt_override=None)`

---

#### Internal structure and operation of `msd_analyzer.py`
`msd_analyzer.py` is a module that exposes a data-class for results and multiple functions for MSD calculation. The logic is divided into reusable units, not a single monolithic class.

Main components:
- `MSDResult` (data-class): collects tau, msd, tracks_per_lag, dt_used, n_max, total number of tracks, and maximum length.
- Support functions:
	- `determine_maximum_lag_steps(N_max, fraction)`: converts the chosen fraction to a maximum number of steps; always limits to 1..N_max-1. If the fraction is absent or invalid, 1.0 is used (complete interval).
	- `build_tau_array(K, dt)`: builds temporal lags tau = [1..K] * dt.
	- `average_across_trajectories(list_of_arrays)`: averages (ignoring NaN) the lag contributions among various tracks and counts how many tracks contributed to each lag.

eaMSD calculation (ensemble, displacement relative to start):
1) For each track and for each lag n, a single displacement relative to the initial point is calculated: (x[n] - x[0])² + (y[n] - y[0])².
2) Per-track vectors (dimension K) are stacked and averaged across tracks for each lag (nan-safe).
3) Steps n are converted to times with tau = n * Δt, where Δt is the track's time-step.

TAMSD calculation (time-averaged, single track):
- For each lag n, the average is calculated over all temporal windows within the track: average of (x[i+n]-x[i])² + (y[i+n]-y[i])² for i = 0..N-n-1.
- The Δt used is the track's time-step
- The output is an MSDResult relative to a single track (tracks_per_lag = 1).

Important notes:
- `max_lag_fraction` is optional; if omitted, the entire possible interval (N_max-1) is used. If provided, it must be 0 < f ≤ 1 and is always clamped to that range.
- Averages ignore NaN values (tracks too short at high lags), so each lag uses only available contributions.
- A positive intercept of the MSD line is normal in the presence of localization noise or acquisition blur.

### 3) `eamsd_plot.py`
- Reads the CSV, calculates the **eaMSD** and saves a PNG plot on linear axes (MSD vs τ).
- Also displays main information on terminal (Δt, n_max, τ_max, number of tracks).
- Options:
	- `--max-lag-fraction`: fraction of the longest track's span to use (0 < f ≤ 1). If omitted, uses everything (N_max − 1).
	- `--output`: output image file name (PNG).

---

### 4) `tamsd_plot.py`
- Calculates and plots the **TAMSD** for **one** track.
- If no track is selected with `--track-id`, uses the first one in order.
- Options:
	- `--track-id`: ID of the track to analyze (as it appears in the CSV).
	- `--max-lag-fraction`: fraction of the selected track's span (0 < f ≤ 1). If omitted, uses everything (N − 1).
	- `--output`: output image file name (PNG).
---

## How to Use (Example Commands)

Below are some ready-to-copy examples for PowerShell (Windows). Replace the paths with your CSV files.

### A) Quick eaMSD calculation (text only in terminal)
```powershell
python .\msd_analyzer.py .\sferette240nm_spots.csv
```

to use only 25% of the lag interval (in points):
```powershell
python .\msd_analyzer.py .\sferette240nm_spots.csv --max-lag-fraction 0.25
```

Suggestion: to save text output to file
```powershell
python .\msd_analyzer.py .\sferette240nm_spots.csv > msd_results.txt
```

### B) eaMSD plot (PNG)
```powershell
python .\eamsd_plot.py .\sferette240nm_spots.csv --output ea_plot.png
```

With custom lag fraction (e.g., 0.3):
```powershell
python .\eamsd_plot.py .\sferette240nm_spots.csv --max-lag-fraction 0.3 --output ea_plot.png
```

### C) TAMSD plot for a track
```powershell
python .\tamsd_plot.py .\sferette240nm_spots.csv --output tamsd_plot.png
```

If you know the track ID (e.g., "1"):
```powershell
python .\tamsd_plot.py .\sferette240nm_spots.csv --track-id 1 --output tamsd_track1.png
```

With custom lag fraction (e.g., 0.2):
```powershell
python .\tamsd_plot.py .\sferette240nm_spots.csv --track-id 1 --max-lag-fraction 0.2 --output tamsd_track1.png
```

## Quick Results Interpretation
- τ (tau): delay time (lag) in seconds → X-axis.
- MSD: mean of squared displacements in μm² (if coordinates are in microns) → Y-axis.
- Global Δt: robust estimate of the common time step across tracks (median of Δt per track).
- `--max-lag-fraction`: controls how far to push with lags. If omitted, uses the entire available interval. Increasing it too much can reduce the number of tracks that contribute to the averages at longer lags.




