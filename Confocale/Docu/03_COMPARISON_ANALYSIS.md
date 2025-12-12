# Comparison Analysis - Multi-Experiment Studies

This document covers `compare_msd.py`, the module for comparing MSD results across different experiments, conditions, and time intervals.

## Overview

The comparison module enables sophisticated analysis of multiple datasets to:

1. **Compare different experiments** (e.g., Experiment 2 vs Experiment 8)
2. **Analyze time-step effects** (40min vs 50min vs 60min intervals)
3. **Investigate temporal evolution** (different acquisition conditions)
4. **Statistical validation** (reproducibility across conditions)

## Module: compare_msd.py

**Purpose**: Interactive and automated comparison of ensemble-averaged MSD across multiple CSV files.

### Key Features

**Interactive File Selection**:
- Graphical file browser for easy dataset selection
- Support for multiple file selection
- Automatic file validation and preview

**Flexible Comparison Modes**:
- **Experiment comparison**: Different experimental conditions
- **Time-step comparison**: Various acquisition intervals
- **Batch analysis**: Process multiple file sets automatically

**Advanced Visualization**:
- **Multi-curve plots**: Up to 10 datasets on single plot
- **Color coding**: Automatic legend and curve differentiation
- **Error bars**: Standard error across trajectories
- **Log-log scaling**: Optimal for MSD analysis

**Statistical Analysis**:
- **Diffusion coefficient extraction**: Linear fitting for each dataset
- **Statistical comparison**: Mean values, confidence intervals
- **Reproducibility metrics**: Inter-experiment variability

---

## Usage Modes

### 1. Interactive Mode (Recommended)

```bash
python compare_msd.py
```

**Workflow**:
1. **File Selection**: Graphical browser opens automatically
2. **Multi-selection**: Hold Ctrl/Cmd to select multiple CSV files
3. **Processing**: Automatic MSD calculation for all selected files
4. **Visualization**: Combined plot with all datasets
5. **Analysis**: Statistical summary and comparison metrics

**Interactive Features**:
- **Live preview**: See file contents before confirmation
- **Error handling**: Automatic validation and user feedback
- **Custom labels**: Option to modify legend labels
- **Export options**: Save plots and statistical summaries

### 2. Command Line Mode

```bash
# Compare specific experiments
python compare_msd.py --files "Data/Experiment2_spots_40minstep.csv" "Data/Experiment8_spots_40minstep.csv"

# Compare time intervals for single experiment
python compare_msd.py --pattern "Data/Experiment8_spots_*minstep.csv"

# Batch processing with custom output
python compare_msd.py --batch --output-dir comparison_analysis/
```

### 3. Automated Pattern Matching

```bash
# All Experiment 8 files with different time steps
python compare_msd.py --pattern "*Experiment8*"

# All 40-minute interval files across experiments
python compare_msd.py --pattern "*40minstep*"

# All files in Data directory
python compare_msd.py --pattern "Data/*.csv"
```

---

## Command Line Arguments

### Required Arguments
- **No arguments**: Launches interactive mode with file browser

### Optional Arguments
- `--files FILE1 FILE2 ...` - Specific files to compare
- `--pattern PATTERN` - Glob pattern for file matching
- `--output-dir DIR` - Custom output directory (default: comparison_plots/)
- `--batch` - Non-interactive batch processing mode
- `--no-plot` - Skip plotting, only generate statistics
- `--format {png,pdf,svg}` - Output plot format (default: png)
- `--dpi DPI` - Plot resolution (default: 300)
- `--figsize WIDTH HEIGHT` - Figure size in inches (default: 10x8)

### Advanced Options
- `--fit-range FRACTION` - Fraction of data to use for fitting (default: 0.1)
- `--min-tracks N` - Minimum tracks required per file (default: 10)
- `--exclude-pattern PATTERN` - Files to exclude from analysis
- `--legend-location {best,upper,lower,right}` - Legend position

---

## Output Files

### 1. Comparison Plots
**Location**: `comparison_plots/comparison_YYYYMMDD_HHMMSS.png`

**Content**:
- **Multi-curve MSD plot**: All datasets on log-log scale
- **Error bars**: Standard error across trajectories
- **Legend**: File names and fitted diffusion coefficients
- **Fit lines**: Linear fits in appropriate regime (first 10% of data)
- **Statistics box**: Summary of D values and R² metrics

### 2. Statistical Summary
**Location**: `comparison_plots/comparison_statistics_YYYYMMDD_HHMMSS.txt`

**Content**:
```
MSD Comparison Analysis Summary
Generated: 2024-11-14 14:30:45

Files Analyzed:
1. Data/Experiment2_spots_40minstep.csv (245 trajectories)
2. Data/Experiment8_spots_40minstep.csv (189 trajectories)

Diffusion Coefficients:
File 1: D = (1.23 ± 0.02) × 10⁻¹ μm²/s, R² = 0.9934
File 2: D = (1.44 ± 0.03) × 10⁻¹ μm²/s, R² = 0.9912

Statistical Comparison:
Mean D: (1.34 ± 0.15) × 10⁻¹ μm²/s
Coefficient of Variation: 11.2%
Relative Difference: 17.1%

Quality Metrics:
All R² > 0.99 (excellent fits)
Average trajectory length: 145 points
Total analysis time: 2.34 seconds
```

### 3. Data Export (Optional)
**Location**: `comparison_plots/comparison_data_YYYYMMDD_HHMMSS.csv`

**Content**: Combined MSD data with lag times, MSD values, and errors for each file

---

## Analysis Scenarios

### Scenario 1: Experiment Comparison
**Objective**: Compare diffusion behavior between different experimental conditions

**Files**: 
- `Experiment2_spots_40minstep.csv` 
- `Experiment8_spots_40minstep.csv`

**Expected Results**:
- Similar D values → consistent experimental conditions
- Different D values → condition-dependent effects
- Similar error bars → good reproducibility

**Command**:
```bash
python compare_msd.py --files "Data/Experiment2_spots_40minstep.csv" "Data/Experiment8_spots_40minstep.csv"
```

### Scenario 2: Time-Step Analysis
**Objective**: Investigate effect of acquisition interval on measured diffusion

**Files**: 
- `Experiment8_spots_40minstep.csv` (40-minute intervals)
- `Experiment8_spots_50minstep.csv` (50-minute intervals)  
- `Experiment8_spots_60minstep.csv` (60-minute intervals)

**Expected Results**:
- Consistent D values → robust measurement
- Systematic differences → time-sampling artifacts
- Larger errors at longer intervals → reduced statistics

**Command**:
```bash
python compare_msd.py --pattern "Data/Experiment8_spots_*minstep.csv"
```

### Scenario 3: Complete Dataset Analysis
**Objective**: Comprehensive analysis of all available data

**Files**: All CSV files in Data directory

**Expected Results**:
- Population of D values → experimental variability
- Outliers → identify problematic datasets
- Trends → systematic effects across conditions

**Command**:
```bash
python compare_msd.py --pattern "Data/*.csv" --batch
```

---

## Statistical Interpretation

### Diffusion Coefficient Variability

**Excellent Reproducibility**: CV < 10%
- **Interpretation**: Highly consistent experimental conditions
- **Example**: D = (1.40 ± 0.05) × 10⁻¹ μm²/s across 5 experiments

**Good Reproducibility**: CV = 10-20%
- **Interpretation**: Minor experimental variations
- **Action**: Check for systematic trends, consider experimental factors

**Poor Reproducibility**: CV > 20%
- **Interpretation**: Significant experimental variability
- **Action**: Investigate experimental conditions, check data quality

### R² Quality Assessment

**Excellent Fits**: R² > 0.995
- **Interpretation**: Clean diffusion regime identified
- **Confidence**: High reliability in D extraction

**Good Fits**: R² = 0.99-0.995
- **Interpretation**: Minor deviations from linear model
- **Action**: Consider nonlinear fitting if systematic

**Poor Fits**: R² < 0.99
- **Interpretation**: Significant non-linear behavior
- **Action**: Investigate subdiffusion, confinement, or directed motion

### Inter-Experiment Comparison

**Statistical Tests Available**:
- **Mean comparison**: t-test for significant differences
- **Variance comparison**: F-test for experimental consistency
- **Trend analysis**: Linear regression across conditions

---

## Advanced Features

### 1. Custom Analysis Ranges
```bash
# Use first 20% of data for fitting (instead of default 10%)
python compare_msd.py --fit-range 0.2 --files file1.csv file2.csv
```

### 2. Quality Filtering
```bash
# Require minimum 20 trajectories per file
python compare_msd.py --min-tracks 20 --pattern "Data/*.csv"
```

### 3. Export Control
```bash
# Generate PDF plots and CSV data export
python compare_msd.py --format pdf --export-data --files file1.csv file2.csv
```

### 4. Publication Mode
```bash
# High-resolution plots with custom sizing
python compare_msd.py --dpi 600 --figsize 12 8 --format svg --files file1.csv file2.csv
```

---

## Troubleshooting

### Common Issues

**1. File Selection Problems**:
```
Error: Could not open file browser
```
- **Cause**: Display/GUI issues on server systems
- **Solution**: Use command line mode with `--files` or `--pattern`

**2. Inconsistent File Formats**:
```
Error: Column mismatch between files
```
- **Cause**: Different CSV formats across files
- **Solution**: Verify all files use same TrackMate export format

**3. Memory Issues with Large Datasets**:
```
MemoryError: Cannot load all files simultaneously
```
- **Cause**: Too many large files selected
- **Solution**: Process subsets or use `--batch` mode

**4. Poor Statistical Quality**:
```
Warning: Low R² values detected (R² < 0.99)
```
- **Cause**: Non-linear diffusion behavior or poor data quality
- **Solution**: Check individual files, consider nonlinear fitting

### Performance Optimization

**For Large Datasets**:
- Use pattern matching instead of manual selection
- Enable batch mode for automated processing
- Process subsets if memory constraints exist

**For Publication Quality**:
- Use high DPI settings (600+)
- Choose vector formats (SVG, PDF) for scalability
- Customize figure size for target publication

---

## Integration with Other Modules

### Workflow Integration

**1. Basic Analysis → Comparison**:
```bash
# First: individual analysis
python eamsd_plot.py Data/Experiment8_spots_40minstep.csv

# Then: comparative analysis
python compare_msd.py --pattern "Data/Experiment8_*"
```

**2. Comparison → Detailed Fitting**:
```bash
# First: identify interesting datasets
python compare_msd.py --pattern "Data/*.csv"

# Then: detailed fitting on selected files
python msd_fitting_nonlinear.py Data/selected_file.csv --compare-linear
```

### Next Steps

1. **Quantitative Analysis**: Use `04_MSD_FITTING.md` for detailed parameter extraction
2. **Individual Analysis**: Return to `02_DATA_ANALYSIS_CORE.md` for specific trajectories
3. **Advanced Modeling**: Consider nonlinear fitting if systematic deviations observed

---

## Examples and Case Studies

### Case Study 1: Time-Step Sensitivity
**Dataset**: Experiment 8 with varying intervals

**Results**:
- 40min: D = (1.44 ± 0.03) × 10⁻¹ μm²/s
- 50min: D = (1.39 ± 0.04) × 10⁻¹ μm²/s  
- 60min: D = (1.42 ± 0.05) × 10⁻¹ μm²/s

**Interpretation**: Consistent D values across time intervals → robust measurement

### Case Study 2: Experiment-to-Experiment Variation
**Dataset**: Experiment 2 vs Experiment 8 (40min intervals)

**Results**:
- Exp 2: D = (1.23 ± 0.02) × 10⁻¹ μm²/s
- Exp 8: D = (1.44 ± 0.03) × 10⁻¹ μm²/s

**Interpretation**: 17% difference → investigate experimental conditions

This comparison functionality provides essential tools for validating experimental reproducibility and identifying systematic trends across datasets.