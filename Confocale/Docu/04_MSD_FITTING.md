# MSD Fitting - Quantitative Parameter Extraction

This document covers the five quantitative fitting modules for extracting physical parameters (diffusion coefficient, drift velocity, anomalous exponent) from MSD data.

## Overview

The fitting toolkit provides five complementary approaches:

1. **`msd_fitting.py`** - Linear model fitting (pure normal diffusion)
2. **`msd_fitting_nonlinear.py`** - Nonlinear model fitting (normal diffusion + drift)  
3. **`msd_fitting_anomalous.py`** - Anomalous diffusion fitting (no drift)
4. **`msd_fitting_anomalous_drift.py`** - Anomalous diffusion with drift fitting
5. **`plot_msd_fit.py`** - Visualization for all fitting models

These modules extract quantitative parameters from MSD curves to characterize particle transport mechanisms, from simple Brownian motion to complex anomalous diffusion with drift.

---

## Physical Models

### Linear Model (Pure Normal Diffusion)
**Equation**: `MSD(τ) = 2D·τ`

**Assumptions**:
- Pure Brownian diffusion
- No directed motion or drift
- Linear relationship between MSD and time lag
- Normal diffusion exponent α = 1

**Output**: Diffusion coefficient D [μm²/s]

**Best for**: Simple diffusion systems, initial analysis

### Nonlinear Model (Normal Diffusion + Drift)
**Equation**: `MSD(τ) = 4D·τ + v²·τ²`

**Assumptions**:
- Brownian diffusion + constant drift velocity
- Normal diffusion exponent α = 1
- Additive contributions from diffusive and ballistic motion
- Valid for combined transport mechanisms

**Output**: 
- Diffusion coefficient D [μm²/s]  
- Drift velocity v [μm/s]

**Best for**: Systems with active transport, flow, or directed motion

### Anomalous Diffusion Model (No Drift)
**Equation**: `MSD(τ) = 4D_α·τ^α`

**Assumptions**:
- Anomalous diffusion behavior
- No directed motion or drift
- Power-law relationship between MSD and time lag
- Exponent α characterizes diffusion type

**Output**:
- Generalized diffusion coefficient D_α [μm²/s^α]
- Anomalous exponent α (dimensionless)

**Interpretation of α**:
- **α < 1**: Subdiffusion (confined/hindered motion, crowded environments)
- **α = 1**: Normal Brownian diffusion
- **α > 1**: Superdiffusion (active transport, enhanced mobility)
- **α = 2**: Pure ballistic motion

**Best for**: Complex systems with non-Brownian behavior, viscoelastic media, crowded environments

### Anomalous Diffusion with Drift Model
**Equation**: `MSD(τ) = 4D_α·τ^α + v²·τ²`

**Assumptions**:
- Anomalous diffusion + constant drift velocity
- Combined power-law diffusive and ballistic behavior
- Most general model combining all transport mechanisms

**Output**:
- Generalized diffusion coefficient D_α [μm²/s^α]
- Anomalous exponent α (dimensionless)
- Drift velocity v [μm/s]

**Best for**: Complex systems with both anomalous diffusion and directed motion

---

## Module 1: msd_fitting.py - Linear Fitting

**Purpose**: Extract diffusion coefficient assuming pure Brownian motion.

### Key Features
- **Automatic 10% rule**: Uses first 10% of lag steps (optimal linear regime)
- **Bounded fitting**: Physically realistic parameter bounds
- **Trust Region algorithm**: Robust optimization with bound constraints
- **Uncertainty quantification**: Standard errors from covariance matrix

### Usage

**Command Line**:
```bash
# Basic usage
python msd_fitting.py Data/Experiment8_spots_40minstep.csv

# Custom fit fraction (15% instead of 10%)
python msd_fitting.py Data/Experiment8_spots_40minstep.csv --fit-fraction 0.15

# Custom output directory  
python msd_fitting.py Data/experiment.csv --output-dir custom_fits/

# Adjust parameter bounds
python msd_fitting.py Data/experiment.csv --D-bounds 1e-4 1e-1
```

**As Module**:
```python
from msd_fitting import fit_linear_msd
from data_reader import read_trajectory_data

# Load data and fit
data = read_trajectory_data('Data/Experiment8_spots_40minstep.csv')
result = fit_linear_msd(data, fit_fraction=0.1)

print(f"D = ({result['D']:.3e} ± {result['D_error']:.2e}) μm²/s")
print(f"R² = {result['R_squared']:.4f}")
```

### Command Line Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `csv_file` | str | (required) | Path to trajectory CSV file |
| `--fit-fraction` | float | 0.1 | Fraction of lag steps for fitting (10% rule) |
| `--D-initial` | float | 0.01 | Initial guess for D [μm²/s] |
| `--D-bounds` | float float | 6e-4 1e-1 | Lower/upper bounds for D [μm²/s] |
| `--output-dir` | str | linear_fits | Output directory for results |
| `--min-tracks` | int | 10 | Minimum trajectories required |
| `--max-lag-fraction` | float | None | Max lag as fraction of longest trajectory |

### Output

**Console Output**:
```
Linear MSD Fitting Results:
==========================
  Fit fraction: 10% (14/145 lag steps)
  Diffusion Coeff (D): (1.437e-01 ± 1.82e-03) μm²/s  
  R²: 0.9911
  Fit range: τ ∈ [1.29, 18.10] s
  Number of trajectories: 189
```

**Files Created**:
- `{experiment}_linear_fit_results.txt` - Detailed results summary
- Console output saved for record-keeping

### Physical Interpretation

**For 240nm Particles in 85% Glycerol**:
- **Theoretical D**: 0.006-0.009 μm²/s (Stokes-Einstein)
- **Measured D**: ~0.14 μm²/s (typical result)
- **Discrepancy**: ~15× larger than theoretical

**Possible Explanations**:
1. **Lower glycerol concentration** (e.g., 50% instead of 85%)
2. **Temperature effects** (higher than assumed 25°C)
3. **Particle aggregation** (effective size changes)
4. **Active transport** (need nonlinear model)
5. **Flow/convection** (systematic drift)

**Quality Metrics**:
- **R² > 0.99**: Excellent linear fit
- **R² = 0.95-0.99**: Good fit, minor deviations
- **R² < 0.95**: Poor linear model, consider nonlinear fitting

---

## Module 2: msd_fitting_nonlinear.py - Advanced Fitting

**Purpose**: Extract both diffusion coefficient and drift velocity from MSD data.

### Key Features
- **Automatic interval optimization**: Tests 10%-90% ranges to find optimal fitting interval
- **Velocity pre-analysis**: Calculates trajectory velocities to set intelligent parameter bounds
- **Model comparison**: Optional linear model comparison to quantify improvement
- **Manual control**: Override automatic optimization with specific intervals
- **Comprehensive output**: Statistical summaries, plots, and detailed results

### Physical Model

**Equation**: `MSD(τ) = 4D·τ + v²·τ²`

**Components**:
- **4D·τ**: Diffusive contribution (factor of 4 for 2D independent diffusion)
- **v²·τ²**: Ballistic contribution (constant velocity motion)
- **Combined**: Total displacement variance from both mechanisms

### Usage

**Basic Analysis**:
```bash
# Nonlinear fitting only
python msd_fitting_nonlinear.py Data/Experiment8_spots_40minstep.csv

# Compare with linear model
python msd_fitting_nonlinear.py Data/experiment.csv --compare-linear

# Generate velocity distribution plot
python msd_fitting_nonlinear.py Data/experiment.csv --plot-velocity
```

**Manual Interval Control**:
```bash
# Fit using first 30% of data
python msd_fitting_nonlinear.py Data/experiment.csv --manual-fraction 0.3

# Manual fraction with comparison
python msd_fitting_nonlinear.py Data/experiment.csv --manual-fraction 0.4 --compare-linear
```

**Advanced Options**:
```bash
# Finer interval search (5% steps instead of 10%)
python msd_fitting_nonlinear.py Data/experiment.csv --interval-step 0.05

# Stricter trajectory filtering (≥50 points)
python msd_fitting_nonlinear.py Data/experiment.csv --min-points 50

# Custom parameter bounds
python msd_fitting_nonlinear.py Data/experiment.csv --D-bounds 1e-4 5e-1
```

### Command Line Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `csv_file` | str | (required) | Path to trajectory CSV file |
| `--max-lag-fraction` | float | None | Fraction of longest track for max lag (0<f≤1] |
| `--min-points` | int | 30 | Minimum trajectory length for velocity analysis |
| `--D-initial` | float | 1e-2 | Initial guess for D [μm²/s] |
| `--D-bounds` | float float | 9e-4 1.5e-1 | Lower and upper bounds for D [μm²/s] |
| `--interval-step` | float | 0.10 | Step size for interval search (10% increments) |
| `--manual-fraction` | float | None | Manual fit interval (0.1-0.9), skips optimization |
| `--plot-velocity` | flag | False | Generate velocity distribution histogram |
| `--output-dir` | str | nonlinear_fits | Output directory for results |
| `--compare-linear` | flag | False | Also run linear fit and compare D values |

### Output Analysis

**Velocity Analysis Section**:
```
============================================================
Velocity Analysis (filtering trajectories with ≥30 points)  
============================================================
  Trajectories used: 189/189
  Mean velocity:     0.050 μm/s
  Median velocity:   0.048 μm/s  
  Std deviation:     0.012 μm/s
  Initial guess (v): 0.048 μm/s (median, robust to outliers)
  Bounds for v:      [0.000, 0.086] μm/s (0 to mean+3σ)
============================================================
```

**Interval Optimization**:
```
Testing intervals from 10% to 90% in 10% steps...
   10%: R² = 0.9954, RSS = 1.53e-01, D = 6.31e-02, v = 5.00e-02
   20%: R² = 0.9979, RSS = 1.27e+00, D = 4.78e-02, v = 7.81e-02  
   30%: R² = 0.9990, RSS = 3.12e+00, D = 4.74e-02, v = 7.80e-02
   ...
Optimal interval: 10% (R² = 0.9954, RSS = 1.53e-01)
```

**Final Results**:
```
============================================================
Nonlinear Fit Results:
============================================================
  Optimal interval:          10% (14 steps)
  Diffusion Coeff (D):       (6.31 ± 0.28) × 10⁻² μm²/s
  Drift Velocity (v):        (5.00 ± 0.76) × 10⁻² μm/s  
  R²:                        0.9954
  RSS:                       1.53e-01
  Fit τ range:               [1.29, 18.10] s
============================================================
```

**Model Comparison** (with `--compare-linear`):
```
============================================================
Model Comparison:
============================================================
Linear model (MSD = 2D·τ):
  D = (1.437 ± 0.018) × 10⁻¹ μm²/s, R² = 0.9911

Nonlinear model (MSD = 4D·τ + v²·τ²):  
  D = (6.31 ± 0.28) × 10⁻² μm²/s, R² = 0.9954
  v = (5.00 ± 0.76) × 10⁻² μm/s

Improvement: ΔR² = +0.0043, D_linear/D_nonlinear = 2.28
============================================================
```

### Interpretation Guidelines

**v ≈ 0 (v < 2×v_error)**:
- Pure diffusion dominates  
- Linear and nonlinear models should give similar D (accounting for 2× vs 4× factor)
- Brownian motion behavior

**v > 0 (v ≥ 2×v_error)**:
- Significant drift/directed motion
- Could indicate: active transport, flow, chemical gradients
- Linear model overestimates D by absorbing ballistic term

**R² Improvement**:
- **ΔR² > 0.01**: Strong evidence for nonlinear behavior
- **ΔR² = 0.001-0.01**: Modest improvement, check v significance  
- **ΔR² < 0.001**: Linear model adequate

**Optimal Interval Analysis**:
- **Always 10%**: Well-described by model across time scales
- **Variable optimal**: Time-dependent behavior (confinement, caging)
- **Always 90%**: Possible systematic errors, check data quality

---

## Module 3: plot_msd_fit.py - Visualization

**Purpose**: Generate publication-quality plots for both linear and nonlinear MSD fits.

### Key Features
- **Dual mode support**: Linear (default) or nonlinear (`--nonlinear`) fitting
- **Automatic parameter detection**: Reads from fitting modules or computes directly
- **High-resolution output**: 300 DPI default, customizable
- **Publication ready**: Professional formatting, clear legends, proper axis labels
- **Flexible naming**: Automatic or custom output filenames

### Usage

**Linear Fitting Plots**:
```bash
# Basic linear fit plot
python plot_msd_fit.py Data/Experiment8_spots_40minstep.csv

# Custom fit fraction
python plot_msd_fit.py Data/experiment.csv --fit-fraction 0.15

# Custom output file
python plot_msd_fit.py Data/experiment.csv --output my_linear_fit.png
```

**Nonlinear Fitting Plots**:
```bash
# Basic nonlinear fit plot  
python plot_msd_fit.py --nonlinear Data/Experiment8_spots_40minstep.csv

# Manual fraction specification
python plot_msd_fit.py --nonlinear Data/experiment.csv --manual-fraction 0.3

# Custom output directory
python plot_msd_fit.py --nonlinear Data/experiment.csv --output-dir publication_plots/
```

### Command Line Arguments

**Common Arguments**:
| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `csv_file` | str | (required) | Path to trajectory CSV file |
| `--output` | str | None | Custom output filename (auto-generated if omitted) |
| `--output-dir` | str | linear_fits/ | Output directory for plots |
| `--dpi` | int | 300 | Plot resolution (DPI) |
| `--figsize` | float float | 8 6 | Figure size in inches (width height) |

**Linear Mode Arguments**:
| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--fit-fraction` | float | 0.1 | Fraction of lag steps for fitting |
| `--D-bounds` | float float | 6e-4 1e-1 | Bounds for D parameter |

**Nonlinear Mode Arguments** (with `--nonlinear`):
| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--manual-fraction` | float | None | Manual fit interval (0.1-0.9) |
| `--interval-step` | float | 0.10 | Step size for interval optimization |
| `--D-bounds` | float float | 9e-4 1.5e-1 | Bounds for D parameter |

### Plot Output

**Linear Fit Plot**:
- **Blue circles**: MSD data points (first 10% of lag steps)
- **Red line**: Linear fit `MSD = 2D·τ`
- **Text box**: D value ± error, R², fit parameters
- **Axes**: Log-log scale, proper units (μm², seconds)
- **Title**: Experiment identifier and model type

**Nonlinear Fit Plot**:
- **Blue circles**: MSD data points (optimal interval)
- **Red curve**: Nonlinear fit `MSD = 4D·τ + v²·τ²`
- **Text box**: D ± error, v ± error, R², fit parameters  
- **Axes**: Log-log scale with appropriate range
- **Title**: Experiment identifier and optimal interval

---

## Complete Analysis Workflows

### Workflow 1: Linear Analysis (Simple Diffusion)
```bash
# Step 1: Basic MSD visualization
python eamsd_plot.py Data/Experiment8_spots_40minstep.csv

# Step 2: Linear fitting
python msd_fitting.py Data/Experiment8_spots_40minstep.csv

# Step 3: Publication plot
python plot_msd_fit.py Data/Experiment8_spots_40minstep.csv
```

### Workflow 2: Comprehensive Analysis (Diffusion + Drift)
```bash
# Step 1: Nonlinear analysis with comparison
python msd_fitting_nonlinear.py Data/Experiment8_spots_40minstep.csv --compare-linear --plot-velocity

# Step 2: Generate both plots
python plot_msd_fit.py --nonlinear Data/Experiment8_spots_40minstep.csv
python plot_msd_fit.py Data/Experiment8_spots_40minstep.csv

# Step 3: Review results
# Compare D values, assess v significance, examine R² improvement
```

### Workflow 3: Parameter Exploration  
```bash
# Test different intervals manually
python msd_fitting_nonlinear.py Data/experiment.csv --manual-fraction 0.1
python msd_fitting_nonlinear.py Data/experiment.csv --manual-fraction 0.3
python msd_fitting_nonlinear.py Data/experiment.csv --manual-fraction 0.5

# Generate plots for optimal interval
python plot_msd_fit.py --nonlinear Data/experiment.csv --manual-fraction 0.2
```

---

## Expected Results for 240nm Particles

### Theoretical Framework

**Stokes-Einstein Equation**:
```
D = kT/(6πηr)
```

**For 85% Glycerol Solution**:
- k = 1.381×10⁻²³ J/K (Boltzmann constant)
- T = 298 K (room temperature)  
- η = 125 mPa·s (85% glycerol viscosity)
- r = 120 nm (particle radius)
- **Predicted D**: 0.006-0.009 μm²/s

### Typical Experimental Results

**Linear Fitting**:
- **D_linear**: ~0.14 μm²/s (15× larger than theoretical)
- **R²**: 0.991-0.999
- **Interpretation**: Apparent diffusion includes ballistic contribution

**Nonlinear Fitting**:
- **D_nonlinear**: ~0.063 μm²/s (7× larger than theoretical)  
- **v**: ~0.05 μm/s
- **R²**: 0.995-0.999
- **Improvement**: ΔR² ≈ 0.004

### Physical Interpretation

**D Discrepancy Sources**:
1. **Different glycerol concentration** (actual vs nominal)
2. **Temperature effects** (viscosity decreases ~3% per °C)
3. **Particle aggregation** (effective size changes)
4. **Active processes** (cellular activity, convection)

**Velocity Significance**:
- **v ~ 0.05 μm/s**: Typical for microfluidic flows
- **v > v_error**: Statistically significant drift detected  
- **Timescale analysis**: D/v² ≈ 25 s (crossover from diffusive to ballistic)

---

## Troubleshooting and Best Practices

### Common Issues

**1. Fitting Convergence Failures**:
```
Error: Optimization did not converge
```
**Solutions**:
- Check parameter bounds are realistic
- Verify data quality (no NaN values, reasonable scales)
- Try different initial guesses
- Reduce fitting interval (use smaller fraction)

**2. Unrealistic Parameter Values**:
```
Warning: D value at boundary
```
**Solutions**:
- Expand parameter bounds
- Check experimental conditions match assumptions
- Verify data units (micrometers, seconds)

**3. Poor Fit Quality (R² < 0.95)**:
```
Warning: Low R² detected
```
**Solutions**:
- Inspect MSD curve for anomalies
- Consider alternative models (subdiffusion, confined diffusion)
- Check trajectory length and quality
- Examine residuals for systematic patterns

### Best Practices

1. **Always Start with Comparison**: Use `--compare-linear` to understand model differences
2. **Validate Velocity Significance**: Ensure v > 2×v_error for meaningful drift
3. **Check Time-Scale Consistency**: Test multiple intervals with `--manual-fraction`
4. **Quality Control**: Require R² > 0.99 for reliable parameter extraction
5. **Physical Validation**: Compare with theoretical expectations (Stokes-Einstein)

---

## Integration with Other Modules

### Data Pipeline
```
Raw CSV Data
     ↓ (data_reader.py)
Processed Trajectories  
     ↓ (msd_analyzer.py)
MSD Curves
     ↓ (msd_fitting.py / msd_fitting_nonlinear.py)
Physical Parameters
     ↓ (plot_msd_fit.py)
Publication Plots
```

### Cross-Module Consistency
- **Parameter bounds**: Consistent across fitting and plotting modules
- **Algorithm choices**: Same optimization methods and interval selection
- **Output formats**: Compatible file naming and directory structure
- **Error propagation**: Consistent uncertainty calculations

### Next Steps After Fitting
1. **Compare with literature values** for similar systems
2. **Investigate discrepancies** using time-scale analysis
3. **Validate with control experiments** (known D values)
4. **Consider advanced models** (subdiffusion, confined motion) if needed

---

## Module 4: msd_fitting_anomalous.py - Anomalous Diffusion (No Drift)

**Purpose**: Extract generalized diffusion coefficient and anomalous exponent from MSD data exhibiting non-Brownian behavior.

### Key Features
- **Automatic interval optimization**: Tests 10%-90% ranges to find optimal fitting interval
- **Power-law fitting**: Captures subdiffusion, normal diffusion, and superdiffusion
- **Conservative α bounds**: Default range [0.01, 2.0] covers all physical regimes
- **Model comparison**: Optional comparison with normal diffusion (α=1)
- **Automatic plot generation**: Publication-quality figures with fit parameters
- **Interpretation guidance**: Automatic classification of diffusion type

### Physical Model

**Equation**: `MSD(τ) = 4D_α·τ^α`

**Components**:
- **D_α**: Generalized diffusion coefficient with units [μm²/s^α] (dimension depends on α)
- **α**: Anomalous diffusion exponent (dimensionless)
- **4**: Factor for 2D diffusion (2 independent dimensions)

**Physical Regimes**:
- **α < 0.95**: Subdiffusion - confined motion, crowding, obstacles, viscoelastic media
- **α ≈ 1.0**: Normal diffusion - standard Brownian motion
- **α > 1.05**: Superdiffusion - active transport, flow, enhanced mobility
- **α ≈ 2.0**: Ballistic motion - constant velocity transport

### Usage

**Basic Analysis**:
```bash
# Automatic optimization (tests 10%-90%)
python msd_fitting_anomalous.py Data/240nm_glicerolo50_serie004_251minstep.csv

# Compare with normal diffusion (α=1)
python msd_fitting_anomalous.py Data/experiment.csv --compare-normal
```

**Manual Interval Control**:
```bash
# Fit using first 30% of data
python msd_fitting_anomalous.py Data/experiment.csv --manual-fraction 0.3

# Custom α bounds for subdiffusion regime
python msd_fitting_anomalous.py Data/experiment.csv --alpha-bounds 0.1 0.99
```

**Advanced Options**:
```bash
# Finer interval search (5% steps)
python msd_fitting_anomalous.py Data/experiment.csv --interval-step 0.05

# Custom D_α and α parameters
python msd_fitting_anomalous.py Data/experiment.csv --D-alpha-initial 0.005 --alpha-initial 0.8
```

### Command Line Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `csv_file` | str | (required) | Path to trajectory CSV file |
| `--max-lag-fraction` | float | None | Fraction of longest track for max lag (0<f≤1] |
| `--D-alpha-initial` | float | 1e-2 | Initial guess for D_α [μm²/s^α] |
| `--D-alpha-bounds` | float float | 1e-6 1e2 | Lower/upper bounds for D_α [μm²/s^α] |
| `--alpha-initial` | float | 1.0 | Initial guess for α (default: normal diffusion) |
| `--alpha-bounds` | float float | 0.01 2.0 | Lower/upper bounds for α (subdiffusion to ballistic) |
| `--interval-step` | float | 0.10 | Step size for interval search (10% increments) |
| `--manual-fraction` | float | None | Manual fit interval (0.1-0.9), skips optimization |
| `--output-dir` | str | anomalous_fits | Output directory for results and plots |
| `--compare-normal` | flag | False | Compare with normal diffusion (α fixed at 1.0) |

### Output Analysis

**Interval Optimization**:
```
Testing intervals from 10% to 90% in 10% steps...
   10%: R² = 0.9965, RSS = 5.04e-02, D_α = 4.41e-03, α = 1.5921
   20%: R² = 0.9957, RSS = 8.29e-01, D_α = 8.74e-03, α = 1.3540
   30%: R² = 0.9935, RSS = 6.19e+00, D_α = 5.95e-03, α = 1.4628
   ...
Optimal interval: 10% (R² = 0.9965, RSS = 5.04e-02)
```

**Final Results**:
```
============================================================
Anomalous Diffusion Fit Results:
============================================================
  Optimal interval:                10% (49 steps)
  Diffusion Coeff (D_α):         (4.414131e-03 ± 2.59e-04) μm²/s^α
  Anomalous exponent (α):        (1.5921 ± 0.0219)
  R²:                            0.996474
  RSS:                           5.042357e-02
  Fit τ range:                   [0.37, 18.20] s
  Interpretation:                Superdiffusion (active transport)
============================================================

Plot saved to: anomalous_fits/240nm_glicerolo50_serie004_251minstep_anomalous_fit_10pct.png
```

**Model Comparison** (with `--compare-normal`):
```
============================================================
Normal Diffusion Comparison (α fixed at 1.0)
============================================================

Normal diffusion (α=1, MSD = 2D·τ):
  D = (4.135305e-02 ± 1.01e-03) μm²/s
  R² = 0.924352

Anomalous diffusion (α free, MSD = 4D_α·τ^α):
  D_α = (4.414131e-03 ± 2.59e-04) μm²/s^α
  α = (1.5921 ± 0.0219)
  R² = 0.996474

Anomalous model provides better fit (ΔR² = +0.0721)
============================================================
```

### Interpretation Guidelines

**α Significance**:
- **|α - 1| > 2×α_error**: Statistically significant deviation from normal diffusion
- **|α - 1| < α_error**: Consistent with normal Brownian diffusion
- Check if anomalous model provides significant R² improvement over normal model

**R² Improvement**:
- **ΔR² > 0.05**: Strong evidence for anomalous diffusion
- **ΔR² = 0.01-0.05**: Moderate anomalous behavior
- **ΔR² < 0.01**: Normal diffusion model may be sufficient

**Physical Interpretation**:

**Subdiffusion (α < 1)**:
- **Causes**: Crowding, obstacles, binding, caging, viscoelastic media
- **Examples**: Proteins in cells, particles in polymer networks, confined geometries
- **Time-scale**: D_α/τ^(1-α) increases with time → slowing down

**Superdiffusion (α > 1)**:
- **Causes**: Active transport, flow, ballistic motion, Lévy flights
- **Examples**: Motor-driven transport, convection, enhanced diffusion
- **Time-scale**: D_α·τ^(α-1) increases with time → speeding up

### Generated Files

- **`{experiment}_anomalous_fit_{XX}pct.png`** - Publication-quality plot showing:
  - MSD data points (blue circles)
  - Anomalous fit curve (red line)
  - Fit parameters in text box (D_α, α, R²)
  - Proper axis labels and formatting

---

## Module 5: msd_fitting_anomalous_drift.py - Anomalous Diffusion with Drift

**Purpose**: Extract generalized diffusion coefficient, anomalous exponent, and drift velocity from complex MSD data.

### Key Features
- **Three-parameter fitting**: Simultaneously fits D_α, α, and v
- **Velocity pre-analysis**: Calculates trajectory velocities for intelligent v bounds
- **Most general model**: Captures anomalous diffusion + directed motion
- **Model comparison**: Optional comparison with non-drifting anomalous model
- **Velocity histogram**: Optional visualization of trajectory velocity distribution
- **Automatic plot generation**: Publication-quality figures with all fit parameters

### Physical Model

**Equation**: `MSD(τ) = 4D_α·τ^α + v²·τ²`

**Components**:
- **4D_α·τ^α**: Anomalous diffusion term (power-law in time)
- **v²·τ²**: Ballistic drift term (quadratic in time)
- **Combined**: Total displacement variance from both mechanisms

**Physical Interpretation**:
- Generalizes all previous models
- Can capture subdiffusion + drift (e.g., active transport in crowded environments)
- Can capture superdiffusion + drift (e.g., enhanced diffusion with flow)
- Most flexible model for complex biological systems

### Usage

**Basic Analysis**:
```bash
# Automatic optimization
python msd_fitting_anomalous_drift.py Data/240nm_glicerolo50_serie004_251minstep.csv

# With velocity histogram
python msd_fitting_anomalous_drift.py Data/experiment.csv --plot-velocity

# Compare with non-drifting anomalous model
python msd_fitting_anomalous_drift.py Data/experiment.csv --compare-no-drift
```

**Manual Interval Control**:
```bash
# Fit using first 20% of data
python msd_fitting_anomalous_drift.py Data/experiment.csv --manual-fraction 0.2

# Custom α bounds with comparison
python msd_fitting_anomalous_drift.py Data/experiment.csv --alpha-bounds 0.1 1.5 --compare-no-drift
```

**Advanced Options**:
```bash
# Stricter trajectory filtering
python msd_fitting_anomalous_drift.py Data/experiment.csv --min-points 50

# Custom interval search
python msd_fitting_anomalous_drift.py Data/experiment.csv --interval-step 0.05

# All features enabled
python msd_fitting_anomalous_drift.py Data/experiment.csv --plot-velocity --compare-no-drift --interval-step 0.1
```

### Command Line Arguments

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `csv_file` | str | (required) | Path to trajectory CSV file |
| `--max-lag-fraction` | float | None | Fraction of longest track for max lag (0<f≤1] |
| `--min-points` | int | 30 | Minimum trajectory length for velocity analysis |
| `--D-alpha-initial` | float | 1e-2 | Initial guess for D_α [μm²/s^α] |
| `--D-alpha-bounds` | float float | 1e-6 1e2 | Lower/upper bounds for D_α [μm²/s^α] |
| `--alpha-initial` | float | 1.0 | Initial guess for α (default: normal diffusion) |
| `--alpha-bounds` | float float | 0.01 2.0 | Lower/upper bounds for α |
| `--interval-step` | float | 0.10 | Step size for interval search (10% increments) |
| `--manual-fraction` | float | None | Manual fit interval (0.1-0.9), skips optimization |
| `--plot-velocity` | flag | False | Generate velocity distribution histogram |
| `--output-dir` | str | anomalous_drift_fits | Output directory for results and plots |
| `--compare-no-drift` | flag | False | Compare with non-drifting anomalous model |

### Output Analysis

**Velocity Analysis Section**:
```
============================================================
Velocity Analysis (filtering trajectories with ≥30 points)
============================================================
  Trajectories used: 36/36
  Mean velocity:     0.320430 μm/s
  Median velocity:   0.316114 μm/s
  Std deviation:     0.027451 μm/s
  Initial guess (v): 0.316114 μm/s
  Bounds for v:      [0.000000, 0.402783] μm/s
============================================================
```

**Interval Optimization**:
```
Testing intervals from 10% to 90% in 10% steps...
   10%: R² = 0.9972, RSS = 4.06e-02, D_α = 7.64e-03, α = 0.9802, v = 6.24e-02
   20%: R² = 0.9957, RSS = 8.29e-01, D_α = 8.74e-03, α = 1.3540, v = 6.51e-13
   30%: R² = 0.9948, RSS = 4.93e+00, D_α = 1.78e-02, α = 0.9545, v = 4.16e-02
   ...
Optimal interval: 10% (R² = 0.9972, RSS = 4.06e-02)
```

**Final Results**:
```
============================================================
Anomalous Diffusion with Drift Fit Results:
============================================================
  Optimal interval:                10% (49 steps)
  Diffusion Coeff (D_α):         (7.642423e-03 ± 1.67e-03) μm²/s^α
  Anomalous exponent (α):        (0.9802 ± 0.1583)
  Drift velocity (v):            (6.238181e-02 ± 3.41e-03) μm/s
  R²:                            0.997162
  RSS:                           4.058596e-02
  Fit τ range:                   [0.37, 18.20] s
  Interpretation:                Normal Brownian diffusion with drift
============================================================

Plot saved to: anomalous_drift_fits/240nm_glicerolo50_serie004_251minstep_anomalous_drift_fit_10pct.png
```

**Model Comparison** (with `--compare-no-drift`):
```
============================================================
Non-Drifting Anomalous Model Comparison
============================================================

Non-drifting model (MSD = 4D_α·τ^α):
  D_α = (4.414131e-03 ± 2.59e-04) μm²/s^α
  α = (1.5921 ± 0.0219)
  R² = 0.996474

Drifting model (MSD = 4D_α·τ^α + v²·τ²):
  D_α = (7.642423e-03 ± 1.67e-03) μm²/s^α
  α = (0.9802 ± 0.1583)
  v = (6.238181e-02 ± 3.41e-03) μm/s
  R² = 0.997162

Both models fit similarly well (ΔR² = +0.0007)
============================================================
```

### Interpretation Guidelines

**Model Selection**:
1. **Compare R² values**: Does drift term significantly improve fit?
2. **Check v significance**: Is v > 2×v_error?
3. **Examine α change**: Does α change significantly between models?
4. **Physical reasoning**: Does drift make sense in your system?

**Drift Significance**:
- **v ≥ 3×v_error**: Strong evidence for directed motion
- **v = 1-3×v_error**: Moderate drift, check R² improvement
- **v < v_error**: Drift term not significant, use non-drifting model

**α Interpretation with Drift**:
- **α ≈ 1, v > 0**: Normal diffusion + drift (most common)
- **α < 1, v > 0**: Subdiffusion + drift (e.g., crowded environment with flow)
- **α > 1, v ≈ 0**: Apparent superdiffusion may be due to missing drift term
  - Refit with drift model to separate contributions

**Common Scenarios**:

**Scenario 1: Apparent Superdiffusion**
- Non-drifting model: α = 1.59 ± 0.02, R² = 0.996
- Drifting model: α = 0.98 ± 0.16, v = 0.062 ± 0.003, R² = 0.997
- **Interpretation**: Superdiffusion was due to unmodeled drift

**Scenario 2: True Subdiffusion with Flow**
- Non-drifting model: α = 0.65 ± 0.03, R² = 0.985
- Drifting model: α = 0.68 ± 0.05, v = 0.045 ± 0.008, R² = 0.993
- **Interpretation**: Genuine subdiffusion in presence of drift

**Scenario 3: Drift Dominates**
- Non-drifting model: α = 1.95 ± 0.05, R² = 0.980
- Drifting model: α = 1.02 ± 0.25, v = 0.180 ± 0.012, R² = 0.998
- **Interpretation**: Ballistic behavior captured by v²τ² term

### Generated Files

- **`{experiment}_anomalous_drift_fit_{XX}pct.png`** - Plot with:
  - MSD data points (blue circles)
  - Anomalous drift fit curve (red line)
  - All fit parameters (D_α, α, v, R²)
  
- **`{experiment}_velocity_distribution.png`** (with `--plot-velocity`) - Histogram showing:
  - Velocity distribution from trajectories
  - Mean and median markers
  - Number of trajectories used

---

## Complete Analysis Workflows for Anomalous Diffusion

### Workflow 1: Initial Characterization
```bash
# Step 1: Quick visual check
python eamsd_plot.py Data/experiment.csv

# Step 2: Test normal diffusion
python msd_fitting.py Data/experiment.csv

# Step 3: Test anomalous diffusion
python msd_fitting_anomalous.py Data/experiment.csv --compare-normal

# Decision: If ΔR² > 0.01, proceed with anomalous models
```

### Workflow 2: Anomalous Diffusion Analysis
```bash
# Step 1: Non-drifting anomalous fit
python msd_fitting_anomalous.py Data/experiment.csv --compare-normal

# Step 2: Add drift term
python msd_fitting_anomalous_drift.py Data/experiment.csv --compare-no-drift --plot-velocity

# Step 3: Compare all models
# Review R², parameter errors, physical interpretation
```

### Workflow 3: Comprehensive Parameter Study
```bash
# Test multiple intervals for non-drifting model
python msd_fitting_anomalous.py Data/experiment.csv --manual-fraction 0.1
python msd_fitting_anomalous.py Data/experiment.csv --manual-fraction 0.3
python msd_fitting_anomalous.py Data/experiment.csv --manual-fraction 0.5

# Test with drift at optimal interval
python msd_fitting_anomalous_drift.py Data/experiment.csv --manual-fraction 0.3 --plot-velocity

# Compare stability of α and v across intervals
```

### Workflow 4: Model Selection Decision Tree
```
Start: MSD curve visualization
  ↓
Test normal diffusion (msd_fitting.py)
  ↓
R² > 0.99? → YES → Linear behavior confirmed
  ↓ NO
Test anomalous (msd_fitting_anomalous.py --compare-normal)
  ↓
ΔR² > 0.01? → NO → Stick with normal diffusion
  ↓ YES
Anomalous behavior detected
  ↓
Add drift term (msd_fitting_anomalous_drift.py --compare-no-drift)
  ↓
v significant AND ΔR² > 0.005? → YES → Use anomalous + drift
  ↓ NO
Use anomalous without drift
```

---

## Model Comparison Summary

| Model | Equation | Parameters | Best For | When to Use |
|-------|----------|------------|----------|-------------|
| **Linear** | 2D·τ | D | Simple Brownian motion | R² > 0.99, no drift |
| **Nonlinear** | 4D·τ + v²·τ² | D, v | Normal diffusion + drift | α ≈ 1, significant drift |
| **Anomalous** | 4D_α·τ^α | D_α, α | Non-Brownian behavior | α ≠ 1, no drift |
| **Anomalous + Drift** | 4D_α·τ^α + v²·τ² | D_α, α, v | Complex transport | α ≠ 1 AND drift |

### Decision Criteria

**Use Linear Model When**:
- Visual inspection shows linear MSD vs τ
- R² > 0.99 with linear fit
- No evidence of drift or anomalous behavior
- Simple system characterization needed

**Use Nonlinear Model When**:
- Significant drift detected (v > 2×v_error)
- MSD shows upward curvature at longer times
- Linear model R² < 0.95
- α ≈ 1 in anomalous fits

**Use Anomalous Model When**:
- α significantly different from 1 (|α-1| > 2×α_error)
- Improves R² by > 0.01 over normal diffusion
- No significant drift component
- Studying subdiffusion/superdiffusion mechanisms

**Use Anomalous + Drift Model When**:
- Both α ≠ 1 AND v > 0 significantly
- Improves R² over anomalous-only model
- Complex biological or active matter systems
- Need to separate diffusive and drift contributions

---

## Troubleshooting Anomalous Diffusion Fits

### Common Issues

**1. α at Boundary Values**:
```
Warning: α = 0.01 (at lower bound)
Warning: α = 2.00 (at upper bound)
```
**Solutions**:
- Expand α bounds if physically justified
- Check if model is appropriate for data
- Examine residuals for systematic errors
- May indicate need for different model

**2. Large α Uncertainty**:
```
α = 1.2 ± 0.8 (large error)
```
**Solutions**:
- Use longer fitting interval (increase fraction)
- Filter for longer trajectories (increase --min-points)
- Reduce interval step for finer optimization
- Check data quality and trajectory lengths

**3. Inconsistent α Between Models**:
```
Non-drifting: α = 1.59 ± 0.02
Drifting: α = 0.98 ± 0.16
```
**Interpretation**:
- **Not a problem!** This indicates drift was causing apparent superdiffusion
- Use drifting model result
- Check if v is significant (v > 2×v_error)

**4. Poor Convergence**:
```
RuntimeError: Fitting failed at fraction 0.5
```
**Solutions**:
- Reduce maxfev parameter (already set to 5000-10000)
- Try different initial guesses
- Use manual fraction with smaller values (0.1-0.3)
- Check data for NaN or infinite values

### Best Practices for Anomalous Diffusion

1. **Always Compare Models**: Use `--compare-normal` and `--compare-no-drift` flags
2. **Check Parameter Significance**: Require parameter > 2×error for significance
3. **Validate Physically**: Does α value make sense for your system?
4. **Test Multiple Intervals**: Use both automatic and manual fractions
5. **Examine Residuals**: Look for systematic patterns in MSD - fit curve
6. **Consider Time Scales**: Does α vary with fitting interval? May indicate time-dependent behavior

---

## Advanced Topics

### Time-Dependent Anomalous Exponent

If α changes significantly with fitting interval:
- **Early times different from late times**: May indicate multiple regimes
- **Solution**: Fit different time windows separately
- **Example**: Early subdiffusion (caging) → late normal diffusion (free)

### Confined Anomalous Diffusion

For confined geometries with anomalous behavior:
- MSD may plateau at long times
- Use shorter fitting intervals (10-30%)
- Consider specialized confined diffusion models

### Distinguishing Superdiffusion from Drift

**Test**:
1. Fit without drift: If α > 1, could be superdiffusion OR hidden drift
2. Fit with drift: If α drops to ~1 and v becomes significant → drift
3. If α stays > 1 with drift → genuine superdiffusion + drift

**Example**:
```bash
# Test 1
python msd_fitting_anomalous.py Data/experiment.csv
# Result: α = 1.6 ± 0.03

# Test 2
python msd_fitting_anomalous_drift.py Data/experiment.csv
# Result: α = 1.0 ± 0.15, v = 0.08 ± 0.004

# Conclusion: Apparent superdiffusion was due to drift
```

---

## Integration with Analysis Pipeline

### Complete Pipeline
```
Raw Trajectory CSV
     ↓ (data_reader.py)
Processed Trajectories
     ↓ (msd_analyzer.py)
Ensemble MSD Curve
     ↓ (Try all models)
  ├─→ msd_fitting.py (normal diffusion)
  ├─→ msd_fitting_nonlinear.py (normal + drift)
  ├─→ msd_fitting_anomalous.py (anomalous)
  └─→ msd_fitting_anomalous_drift.py (anomalous + drift)
     ↓ (Compare R², parameters)
Select Best Model
     ↓ (Physical interpretation)
Publication Plots + Parameter Extraction
```

This comprehensive fitting toolkit provides the quantitative foundation for extracting meaningful physical parameters from particle tracking experiments, enabling direct comparison with theoretical predictions and literature values.
