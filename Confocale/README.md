# Confocale — MSD Analysis

Python toolkit for Mean Squared Displacement analysis of particle trajectories from confocal microscopy.

## Requirements

```
pip install numpy scipy pandas matplotlib
```

## Data format

CSV files with columns: **Track ID**, **X position**, **Y position**, **Time**.  
Header names are flexible (TrackMate-style `POSITION_X/Y/T` and `TRACK_ID` also work).

Place data in:
- `Data/31_10_no_anomalous/` — standard Brownian motion
- `Data/14_11_anomalous/`    — anomalous diffusion (glycerol)

## Scripts

### `run_all.py` — batch analysis (start here)

Runs all MSD plots and fits on every CSV in both datasets:

```
python run_all.py
```

Output:
- `Results/<dataset>/eamsd/` and `tamsd/`  — MSD plot SVGs at lag fractions 10 %, 25 %, 50 %, 100 %
- `Results/<dataset>/linear_fits/` and `nonlinear_fits/` — fit plot SVGs
- `Docu/fits_<dataset>_results.csv` and `.md` — summary tables

### `plot_msd.py` — single MSD plot

```
python plot_msd.py <csv> --mode eamsd [--max-lag-fraction 0.25] [--output out.svg]
python plot_msd.py <csv> --mode tamsd [--track-id ID] [--max-lag-fraction 0.25]
```

### `fit_msd.py` — single fit

```
python fit_msd.py <csv> --model linear           [--fit-fraction 0.10]
python fit_msd.py <csv> --model nonlinear
python fit_msd.py <csv> --model anomalous
python fit_msd.py <csv> --model anomalous_drift   [--output-dir fits/]
```

### `compare_msd.py` — interactive multi-experiment comparison

```
python compare_msd.py
```

Scans `Data/` for files matching `*_spots_<N>minstep.csv`, lets you pick 4 timesteps for EA-MSD overlay and one trajectory for TA-MSD truncation study.

### `check_ergodicity.py` — ergodicity hypothesis check

```
python check_ergodicity.py [--lag-fraction 0.25] [--compare-fraction 0.10]
```

For each CSV in the non-anomalous dataset, overlays EA-MSD and ensemble-averaged TA-MSD on the same plot and computes relative deviation metrics in the linear region. Output in `Results/no_anomalous/ergodicity/`.

### `test_msd_simulation.py` — validation tests

```
python test_msd_simulation.py
```

Simulates Brownian motion and anomalous diffusion to verify that MSD computation and fitting recover known parameters.

### `diffusion_statistics.py` — rigorous D extraction and statistics

```
python diffusion_statistics.py
```

Extracts the diffusion coefficient D from the non-anomalous dataset using three estimators:

1. **Per-track taMSD** — fits each individual track's time-averaged MSD → distribution of D_i values
2. **Per-file eaMSD** — fits the ensemble-averaged MSD of each CSV file → one D per file
3. **Per-file ⟨taMSD⟩** — fits the ensemble-averaged taMSD of each CSV file → one D per file (ergodicity check)

Each estimator is fitted with both linear (`MSD = 4Dτ`) and nonlinear/drift (`MSD = 4Dτ + v²τ²`) models at 10 % and 25 % lag fractions. Compares measured D against the Einstein–Stokes prediction via one-sample t-tests.

Output:
- `Results/no_anomalous/diffusion_statistics/` — D histograms (SVG) and comparison plot
- `Results/no_anomalous/linear_fits/` and `nonlinear_fits/` — fit plot SVGs
- `Docu/diffusion_statistics_results.csv` — summary table
- `Docu/diffusion_statistics_tamsd_per_track.csv` — per-track D values
- `Docu/diffusion_statistics_eamsd_per_file.csv` — per-file eaMSD D values
- `Docu/diffusion_statistics_ens_tamsd_per_file.csv` — per-file ⟨taMSD⟩ D values

### `comp_glycerol_viscosity.py` — theoretical D (Einstein–Stokes)

Computes the theoretical diffusion coefficient using the Cheng 2008 viscosity model for water–glycerol mixtures. Used internally by `diffusion_statistics.py`.

## Diffusion models

| Model | Equation | Parameters | Dataset |
|---|---|---|---|
| Linear | MSD = 4Dτ | D | no_anomalous |
| Nonlinear (drift) | MSD = 4Dτ + v²τ² | D, v | no_anomalous |
| Anomalous | MSD = 4D_α τ^α | D_α, α | anomalous |
| Anomalous + drift | MSD = 4D_α τ^α + v²τ² | D_α, α, v | anomalous |

Models 2–4 automatically test fit intervals from 10 % to 90 % and select the interval with minimum reduced χ².

## Module structure

| File | Role |
|---|---|
| `data_reader.py` | CSV → `Trajectory` objects |
| `msd_analyzer.py` | EA-MSD and TA-MSD computation |
| `msd_fitting.py` | All four fitting models + shared utilities |
| `plot_msd.py` | MSD plotting CLI |
| `fit_msd.py` | Fitting CLI with plot output |
| `run_all.py` | Batch orchestrator |
| `compare_msd.py` | Interactive multi-experiment comparison |
| `check_ergodicity.py` | Ergodicity check (EA-MSD vs TA-MSD) |
| `diffusion_statistics.py` | Rigorous D extraction, histograms, and t-tests |
| `comp_glycerol_viscosity.py` | Einstein–Stokes D_theory (Cheng 2008) |
