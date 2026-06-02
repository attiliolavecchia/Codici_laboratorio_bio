# Figure Coverage Checklist (Exam Delivery)

Date: 2026-06-02
Scope: verification of SVG files referenced by LaTeX sources against generated outputs in `Results/`.

Checked LaTeX files:
- `Latex/Data_Analysis.tex` (non-anomalous diffusion)
- `Latex/Anomalous_Data_Analysis.tex` (anomalous diffusion)

## Summary
- Non-anomalous references found: 10
- Anomalous references found: 32
- Missing references after fix: 0

## Key figure groups requested by exam delivery

### Non-anomalous diffusion
- eaMSD plots: present
- taMSD plots: present
- Linear and nonlinear fits: present
- Residuals and fit diagnostics: present
- Diffusion coefficient histogram(s): present
- Ergodicity comparison (eaMSD vs <taMSD>): present

### Anomalous diffusion
- eaMSD and taMSD plots: present
- Drift-corrected comparisons (eaMSD and taMSD): present in `Results/anomalous/drift_comparison/`
- Anomalous fits (eaMSD and ensemble taMSD): present
- Residual plots for anomalous fits: present
- Histograms of `D_alpha`: present
- Histograms of apparent `D` with fixed alpha=1: present
- eaMSD vs <taMSD> alpha comparison plots: present
- Van Hove distributions and beta-vs-tau plots: present

## Fix applied
In `Latex/Anomalous_Data_Analysis.tex` two references used filenames with `_2` suffix:
- `beta_vs_tau_glic50_2.svg`
- `beta_vs_tau_glic200_2.svg`

Actual files in `Results/` are:
- `beta_vs_tau_glic50.svg`
- `beta_vs_tau_glic200.svg`

The LaTeX references were updated to match existing files.

## Suggested Drive upload structure
- `01_non_anomalous/eamsd_tamsd/`
- `01_non_anomalous/fits/`
- `01_non_anomalous/histograms/`
- `01_non_anomalous/ergodicity/`
- `02_anomalous/eamsd_tamsd/`
- `02_anomalous/drift_comparison/`
- `02_anomalous/fits_and_residuals/`
- `02_anomalous/histograms_Dalpha_D/`
- `02_anomalous/ergodicity_and_alpha/`
- `02_anomalous/van_hove_and_beta/`

This structure mirrors the analysis logic and makes review easier for examiners.
