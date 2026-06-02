# Confocale - MSD Analysis Toolkit

Python tools for Mean Squared Displacement (MSD) analysis of particle trajectories from confocal microscopy.

## Goal
This repository contains the full analysis workflow used for doctoral exam results:
- non-anomalous diffusion analysis
- anomalous diffusion analysis
- fit/diagnostic plots and summary tables
- LaTeX-ready figures used in the report

## Folder layout
- `Data/31_10_no_anomalous/`: CSV trajectories for normal diffusion
- `Data/14_11_anomalous/`: CSV trajectories for anomalous diffusion
- `Results/no_anomalous/`: generated figures and fit outputs (normal diffusion)
- `Results/anomalous/`: generated figures and fit outputs (anomalous diffusion)

## Requirements
Install dependencies in your Python environment:

```bash
pip install numpy scipy pandas matplotlib
```

## Input format
CSV files must contain trajectory coordinates and time. Accepted headers are flexible:
- `Track ID`, `X position`, `Y position`, `Time`
- or TrackMate-style names such as `TRACK_ID`, `POSITION_X`, `POSITION_Y`, `POSITION_T`

## Main scripts
### 1) Full batch run
Run all analyses and figure generation for both datasets:

```bash
python run_all.py
```

### 2) Single-file MSD plotting

```bash
python plot_msd.py <csv> --mode eamsd [--max-lag-fraction 0.25] [--output out.svg]
python plot_msd.py <csv> --mode tamsd [--track-id ID] [--max-lag-fraction 0.25]
```

### 3) Single-file fitting

```bash
python fit_msd.py <csv> --model linear [--fit-fraction 0.10]
python fit_msd.py <csv> --model nonlinear
python fit_msd.py <csv> --model anomalous
python fit_msd.py <csv> --model anomalous_drift [--output-dir fits/]
```

### 4) Ergodicity check (normal diffusion)

```bash
python check_ergodicity.py [--lag-fraction 0.25] [--compare-fraction 0.10]
```

### 5) Diffusion statistics summary (normal diffusion)

```bash
python diffusion_statistics.py
```

Produces:
- diffusion-coefficient histograms
- per-file and per-track summary tables
- linear vs nonlinear model comparison

## Models used
- Linear: `MSD = 4 D tau`
- Linear + drift: `MSD = 4 D tau + v^2 tau^2`
- Anomalous: `MSD = 4 D_alpha tau^alpha`
- Anomalous + drift correction workflows (eaMSD variance correction and ensemble drift subtraction for taMSD)

## Reproducibility notes
- Use the same Python version and package set across runs.
- Keep raw data unchanged in `Data/`.
- Regenerate figures with scripts instead of manual edits.
