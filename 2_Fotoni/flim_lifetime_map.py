"""
FLIM Per-Pixel Lifetime Map
============================

Builds a spatial τ map for each pixel (or rolling NxN superpixel) inside the
Otsu-masked region of a FLIM TIFF stack.  Mono-exponential and amplitude-
weighted bi-exponential (⟨τ⟩) maps are saved as separate PNG files.

Algorithm per pixel (y, x)
--------------------------
1. Sum the photon-count decay curves of the (2r+1)×(2r+1) neighbourhood,
   restricting to pixels inside the Otsu mask (no background contamination).
2. Fit the summed trace with the mono- and bi-exponential models using the
   same Levenberg-Marquardt code as the global fit pipeline.
3. Physical validity: τ outside [0.05, 12] ns is discarded (→ NaN).

Fit window (fixed for all pixels in an image)
---------------------------------------------
  fit_start = global_peak_bin + peak_offset_bins   (default +7 ≈ 0.34 ns,
              to skip the IRF-dominated rising edge)
  fit_end   = last bin in the global trace with count ≥ 1

Output per TIFF  →  img/Prepuzio/lifetime_maps/   (or --output-dir)
  {stem}_tau_mono_{tstart}-{tend}ns_color.png        τ_mono,  inferno colormap
  {stem}_tau_biexp_mean_{tstart}-{tend}ns_color.png  ⟨τ⟩,      inferno colormap
  {stem}_tau_distribution_{tstart}-{tend}ns.png      τ histograms (mean ± σ)
  lifetime_map_summary.csv                           one row per processed file

IMPORTANT – working directory
------------------------------
  The script accepts exactly ONE positional argument: the path to a TIFF file
  or (with --batch) to the folder containing the TIFF files.

  Run from the repository root (Codici/):
      python .\2_Fotoni\flim_lifetime_map.py .\2_Fotoni\data\Prepuzio --batch

  Or cd into 2_Fotoni/ first:
      cd 2_Fotoni
      python flim_lifetime_map.py data/Prepuzio --batch

  Do NOT pass two positional paths – only one input is accepted.

Example commands
----------------
  # Single file (from Codici/)
  python .\2_Fotoni\flim_lifetime_map.py `
      .\2_Fotoni\data\Prepuzio\exp7_50mW_200frame_RAW_Intensity.tif

  # All Prepuzio files, default 3×3 binning (from Codici/)
  python .\2_Fotoni\flim_lifetime_map.py .\2_Fotoni\data\Prepuzio --batch

  # All Prepuzio files, single-pixel (1×1, lower threshold) (from Codici/)
  python .\2_Fotoni\flim_lifetime_map.py .\2_Fotoni\data\Prepuzio `
      --batch --bin-radius 0 --min-counts 15

  # Custom output folder
  python .\2_Fotoni\flim_lifetime_map.py .\2_Fotoni\data\Prepuzio `
      --batch --output-dir .\2_Fotoni\img\my_maps

Key options
-----------
  --bin-radius INT    Superpixel half-width: 0=1×1 (single px, fast but noisy),
                      1=3×3 (default, good SNR), 2=5×5 (smooth, lower resolution)
  --min-counts INT    Min photons in fit window to attempt a fit (default 30;
                      lower to ~15 when using --bin-radius 0)
  --peak-offset-bins  Bins right of global peak to skip for IRF (default 7)
  --colormap STR      Matplotlib colormap for colour maps (default: inferno)
  --laser-rate FLOAT  Laser rep. rate in MHz (default: 80)
  --no-biexp          Skip bi-exponential fitting (only τ_mono; ~2× faster)
  --fit-start-ns NS   Fit window start in ns (default: auto from peak + offset)
  --fit-end-ns NS     Fit window end in ns (default: auto from last bin with signal)
  --batch             Process every .tif/.tiff in the input folder
  --output-dir PATH   Override default output directory
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from tqdm import tqdm

# ─── Make sure sibling modules are importable ─────────────────────────────────
_HERE = Path(__file__).parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from flim_tiff_fit import ( 
    DEFAULT_LASER_REP_RATE_MHZ,
    _build_mask,
    _remove_hot_pixels_2d,
    _remove_hot_pixels_stack,
    detect_trailing_artifact_cutoff,
)
from flim_exponential_fit import ( 
    fit_bi_exponential,
    fit_mono_exponential,
)

# ─── Constants ────────────────────────────────────────────────────────────────
_DEFAULT_OUTPUT_DIR = _HERE / "img" / "Prepuzio" / "lifetime_maps"
_TAU_MIN_NS: float = 0.05   # lower physical validity bound (ns)
_TAU_MAX_NS: float = 12.0   # upper bound ≈ laser period (ns)
_TAU_BOUND_TOL: float = 0.1  # values within this tolerance of either bound are
                              # treated as boundary-hitting (degenerate) fits


# =============================================================================
# Helper: colormap compatible with all matplotlib ≥ 3.2
# =============================================================================

def _get_cmap(name: str):
    """Return a *copy* of the named colormap (compatible with mpl ≥ 3.2)."""
    try:
        # matplotlib ≥ 3.5
        cmap = matplotlib.colormaps[name].copy()
    except AttributeError:
        import copy
        cmap = copy.copy(plt.cm.get_cmap(name))
    return cmap


# =============================================================================
# Superpixel trace extraction
# =============================================================================

def extract_superpixel_trace(
    stack: np.ndarray,
    y: int,
    x: int,
    radius: int,
    mask: np.ndarray,
) -> np.ndarray:
    """
    Sum the photon-count decay curves for a rolling (2r+1)×(2r+1) neighbourhood.

    Only pixels that lie inside ``mask`` contribute to the sum so that
    background pixels never dilute the cell signal.

    Parameters
    ----------
    stack  : (N_bins, H, W) float64 array
    y, x   : centre pixel coordinates
    radius : half-width (1 → 3×3 window, 2 → 5×5, …)
    mask   : (H, W) bool, True = inside cell mask

    Returns
    -------
    trace : (N_bins,) float64 — summed decay curve of contributing pixels
    """
    H, W = stack.shape[1], stack.shape[2]
    y0, y1 = max(0, y - radius), min(H, y + radius + 1)
    x0, x1 = max(0, x - radius), min(W, x + radius + 1)

    patch_m = mask[y0:y1, x0:x1]           # (ny, nx) bool
    if not patch_m.any():
        return np.zeros(stack.shape[0], dtype=np.float64)

    patch_s = stack[:, y0:y1, x0:x1]       # (N_bins, ny, nx)
    return patch_s[:, patch_m].sum(axis=1)  # (N_bins,)


# =============================================================================
# Global fit-window determination
# =============================================================================

def compute_global_fit_window(
    stack: np.ndarray,
    mask: np.ndarray,
    time_ns: np.ndarray,
    peak_offset_bins: int = 7,
    min_counts_global: float = 1.0,
) -> Tuple[int, int, float, float]:
    """
    Determine a single fit window from the spatially-aggregated decay trace.

    The start is placed ``peak_offset_bins`` to the *right* of the global
    intensity peak to skip the IRF-dominated rising edge (the same shift used
    in the existing Prepuzio fits which start at t ≈ 2.7–2.9 ns).

    The end is the last bin whose aggregated count is ≥ ``min_counts_global``.
    Both indices are fixed and reused for every per-pixel fit in the image.

    Returns
    -------
    (start_idx, end_idx, start_ns, end_ns)
        ``end_idx`` is exclusive (Python-slice convention).
    """
    if not mask.any():
        raise RuntimeError("Empty mask – no pixels to aggregate.")

    global_trace = stack[:, mask].astype(np.float64).sum(axis=1)
    peak_idx  = int(np.argmax(global_trace))
    start_idx = min(peak_idx + peak_offset_bins, len(global_trace) - 2)

    # Last bin with enough signal, starting from start_idx
    valid = np.where(global_trace[start_idx:] >= min_counts_global)[0]
    end_idx = (start_idx + int(valid[-1]) + 1) if valid.size else len(global_trace)
    end_idx = min(end_idx, len(global_trace))

    # Guarantee ≥ 2 fitting points
    if end_idx - start_idx < 2:
        end_idx = min(len(global_trace), start_idx + 2)

    return start_idx, end_idx, float(time_ns[start_idx]), float(time_ns[end_idx - 1])


# =============================================================================
# Per-pixel fitting loop
# =============================================================================

def compute_lifetime_maps(
    stack: np.ndarray,
    mask: np.ndarray,
    time_ns: np.ndarray,
    start_idx: int,
    end_idx: int,
    bin_radius: int = 1,
    min_counts: int = 30,
    fit_biexp: bool = True,
    tau_max_ns: float = 12.0,
) -> Dict[str, np.ndarray]:
    """
    Fit mono- and (optionally) bi-exponential decays for every masked pixel.

    For each pixel the trace used for fitting is the sum of photon counts over
    the (2·bin_radius+1)² neighbourhood intersected with the mask.  Pixels with
    fewer than ``min_counts`` photons in the fit window are skipped (→ NaN).

    Parameters
    ----------
    stack        : (N_bins, H, W) float64, hot-pixel cleaned
    mask         : (H, W) bool
    time_ns      : (N_bins,) time axis in ns
    start_idx, end_idx : fit window (exclusive end), from compute_global_fit_window
    bin_radius   : rolling-window half-width (1 → 3×3, 2 → 5×5)
    min_counts   : minimum total counts in the fit window
    fit_biexp    : whether to run the bi-exponential fit

    Returns
    -------
    dict with keys 'tau_mono', 'tau_avg', 'tau1', 'tau2':
        each an (H, W) float64 array, NaN where the fit is absent / failed /
        returned a non-physical value.
    """
    H, W = stack.shape[1], stack.shape[2]
    tau_mono = np.full((H, W), np.nan, dtype=np.float64)
    tau_avg  = np.full((H, W), np.nan, dtype=np.float64)
    tau1_map = np.full((H, W), np.nan, dtype=np.float64)
    tau2_map = np.full((H, W), np.nan, dtype=np.float64)

    time_fit = time_ns[start_idx:end_idx]
    ys, xs   = np.where(mask)

    for y, x in tqdm(zip(ys, xs), total=len(ys), desc="Fitting pixels", unit="px",
                     dynamic_ncols=True):
        trace     = extract_superpixel_trace(stack, y, x, bin_radius, mask)
        trace_fit = trace[start_idx:end_idx]

        if trace_fit.sum() < min_counts:
            continue

        # ── Mono-exponential ─────────────────────────────────────────────
        mono_res = None
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mono_res = fit_mono_exponential(time_fit, trace_fit, counts=trace_fit)
            tau_m = float(mono_res.tau)
            if _TAU_MIN_NS <= tau_m <= _TAU_MAX_NS:
                tau_mono[y, x] = tau_m
            else:
                mono_res = None   # non-physical → don't seed bi-exp with it
        except Exception:
            pass

        if not fit_biexp:
            continue

        # ── Bi-exponential (seeded from mono result when available) ───────
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                bi_res = fit_bi_exponential(
                    time_fit, trace_fit,
                    mono_result=mono_res,
                    counts=trace_fit,
                    tau_max_ns=tau_max_ns,
                )
            t1    = float(bi_res.tau1)
            t2    = float(bi_res.tau2)
            t_avg = float(bi_res.tau_avg)
            # Require both individual τs to be physical (tau_avg then is too)
            if (_TAU_MIN_NS <= t1 <= tau_max_ns and
                    _TAU_MIN_NS <= t2 <= tau_max_ns):
                tau_avg[y, x]  = t_avg
                tau1_map[y, x] = t1
                tau2_map[y, x] = t2
        except Exception:
            pass

    return {
        "tau_mono": tau_mono,
        "tau_avg":  tau_avg,
        "tau1":     tau1_map,
        "tau2":     tau2_map,
    }


# =============================================================================
# Plotting helpers
# =============================================================================

def _percentile_range(
    arr: np.ndarray, lo: float = 2.0, hi: float = 98.0
) -> Tuple[float, float]:
    """Return (vmin, vmax) from percentiles of the finite values."""
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return 0.0, 1.0
    return float(np.percentile(finite, lo)), float(np.percentile(finite, hi))


def save_lifetime_plot(
    intensity: np.ndarray,
    tau_map: np.ndarray,
    output_path: Path,
    colormap: str,
    vmin: float,
    vmax: float,
    title: str = "",
    dpi: int = 200,
) -> None:
    """
    Save a single PNG: grayscale intensity background with a τ-map overlay.

    NaN pixels in ``tau_map`` are rendered fully transparent so the
    intensity background shows through wherever no valid fit exists.

    The overlay is built as an RGBA image where:
      - colour encodes τ (mapped through ``colormap``)
      - alpha = 0.65 on valid pixels, 0.0 on NaN
    This approach is compatible with all matplotlib ≥ 3.2 and avoids
    mutating the global colormap registry.
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.axis("off")

    # Grayscale intensity background
    ax.imshow(intensity, cmap="gray", interpolation="nearest")

    # Build RGBA overlay manually so NaN → alpha=0 (transparent)
    norm    = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    cmap_fn = _get_cmap(colormap)
    tau_rgba = cmap_fn(norm(np.nan_to_num(tau_map, nan=vmin)))  # (H, W, 4)
    # Hard edge: alpha = 0.85 on valid pixels, 0.0 outside the mask.
    tau_rgba[:, :, 3] = np.where(np.isfinite(tau_map), 0.85, 0.0)

    im = ax.imshow(tau_rgba, interpolation="nearest")

    # Colorbar — attach a ScalarMappable so the bar shows the τ scale
    sm = plt.cm.ScalarMappable(cmap=_get_cmap(colormap), norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("τ (ns)", fontsize=11)
    cbar.ax.tick_params(labelsize=9)

    if title:
        ax.set_title(title, fontsize=8, pad=4)

    fig.tight_layout(pad=0.5)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


def save_tau_histogram(
    tau_map: np.ndarray,
    stats: dict,
    label: str,
    output_path: Path,
    color: str = "#e25822",
    dpi: int = 200,
    exclude_bounds_tol: float = 0.0,
) -> None:
    """
    Save a single-panel histogram of per-pixel τ values with mean ± σ lines.

    Parameters
    ----------
    exclude_bounds_tol : float
        If > 0, exclude histogram bars for values within this tolerance of the
        physical bounds (_TAU_MIN_NS, _TAU_MAX_NS).  The mean/σ lines from
        ``stats`` already reflect filtered values when the caller passes a
        bound-filtered stats dict.
    """
    v = tau_map[np.isfinite(tau_map)]
    n_total = len(v)
    if exclude_bounds_tol > 0.0:
        v = v[
            (v > _TAU_MIN_NS + exclude_bounds_tol) &
            (v < _TAU_MAX_NS - exclude_bounds_tol)
        ]
    n_excl = n_total - len(v)
    fig, ax = plt.subplots(figsize=(5, 4))
    if v.size == 0:
        ax.text(0.5, 0.5, "No valid fits", ha="center", va="center",
                transform=ax.transAxes)
    else:
        ax.hist(v, bins=40, color=color, alpha=0.75, edgecolor="white", linewidth=0.4)
        median = stats.get("median", np.nan)
        mean   = stats["mean"]
        ax.axvline(median, color="black", lw=1.5,
                   label=f"median = {median:.3f} ns")
        ax.axvline(mean,   color="gray",  lw=1.0, ls="--",
                   label=f"mean = {mean:.3f} ns")
        ax.set_xlabel("τ (ns)", fontsize=11)
        ax.set_ylabel("Pixel count", fontsize=11)
        #ax.set_title(label, fontsize=11)
        ax.legend(fontsize=9)
    fig.tight_layout(pad=1.0)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


# =============================================================================
# Single-file pipeline
# =============================================================================

def process_tiff_lifetime_map(
    tiff_path: Path,
    output_dir: Path,
    peak_offset_bins: int = 7,
    min_counts: int = 30,
    bin_radius: int = 1,
    laser_rate: float = DEFAULT_LASER_REP_RATE_MHZ,
    colormap: str = "inferno",
    hot_pixel_k: float = 5.0,
    artifact_spike_factor: float = 2.5,
    fit_biexp: bool = True,
    fit_start_ns: Optional[float] = None,
    fit_end_ns: Optional[float] = None,
    tau_max_ns: float = 12.0,
) -> dict:
    """Full lifetime-map pipeline for one TIFF stack."""
    sep = "═" * 65
    print(f"\n{sep}")
    print(f"  Processing: {tiff_path.name}")
    print(sep)

    # ── 1. Load raw stack ─────────────────────────────────────────────────
    stack_raw = tifffile.imread(str(tiff_path))
    if stack_raw.ndim != 3:
        raise ValueError(
            f"Expected 3-D TIFF stack (N_bins, H, W), got shape {stack_raw.shape}"
        )
    N0, H, W = stack_raw.shape
    print(f"  Stack shape : {N0} bins × {H} × {W} px")

    # ── 2. Trailing-artefact trimming ─────────────────────────────────────
    mean_per_frame = stack_raw.astype(np.float64).mean(axis=(1, 2))
    cutoff, trim_report = detect_trailing_artifact_cutoff(
        mean_per_frame, spike_factor=artifact_spike_factor
    )
    if trim_report["detected"]:
        stack_raw = stack_raw[:cutoff]
        print(
            f"  Artifact trimming: removed {trim_report['n_removed']} frames "
            f"→ {stack_raw.shape[0]} bins"
        )
    else:
        print("  Artifact trimming: nothing detected")

    # ── 3. Hot-pixel removal (frame-by-frame) ─────────────────────────────
    stack = _remove_hot_pixels_stack(stack_raw.astype(np.float64), k=hot_pixel_k)

    # ── 4. Time axis ──────────────────────────────────────────────────────
    N_bins       = stack.shape[0]
    time_range_ns = 1000.0 / laser_rate        # laser period in ns
    dt            = time_range_ns / N_bins      # ns per bin
    time_ns       = np.arange(N_bins, dtype=np.float64) * dt
    print(f"  Time axis   : {N_bins} bins × {dt:.4f} ns/bin = {time_ns[-1]:.2f} ns total")

    # ── 5. Otsu mask ──────────────────────────────────────────────────────
    intensity_proj  = stack.sum(axis=0)
    intensity_clean = _remove_hot_pixels_2d(intensity_proj, k=hot_pixel_k)
    mask, otsu_thr  = _build_mask(intensity_clean)
    n_mask = int(mask.sum())
    print(f"  Otsu mask   : {n_mask} pixels  (threshold = {otsu_thr:.3g})")

    # ── 5b. Save Otsu mask for visual inspection ──────────────────────────
    output_dir.mkdir(parents=True, exist_ok=True)
    _stem_early = tiff_path.stem.replace(" ", "_")

    # Binary mask (white = masked-in, black = excluded)
    fig_m, ax_m = plt.subplots(figsize=(5, 5))
    ax_m.axis("off")
    ax_m.imshow(mask, cmap="gray", vmin=0, vmax=1, interpolation="nearest")
    #ax_m.set_title(f"Otsu mask  (thr={otsu_thr:.3g},  {n_mask} px)", fontsize=8)
    fig_m.tight_layout(pad=0.3)
    _mask_path = output_dir / f"{_stem_early}_otsu_mask.png"
    fig_m.savefig(_mask_path, dpi=200, bbox_inches="tight")
    plt.close(fig_m)
    print(f"  Saved: {_mask_path.name}")

    # Mask overlaid on intensity projection (red = excluded, green = included)
    fig_o, ax_o = plt.subplots(figsize=(5, 5))
    ax_o.axis("off")
    _int_norm = intensity_clean / (intensity_clean.max() + 1e-9)
    _overlay  = np.stack([_int_norm, _int_norm, _int_norm], axis=-1)
    _red_ch = _overlay[:, :, 0].copy(); _red_ch[mask]  = 0.0
    _grn_ch = _overlay[:, :, 1].copy(); _grn_ch[~mask] = 0.0
    _overlay[:, :, 0] = _red_ch
    _overlay[:, :, 1] = _grn_ch
    ax_o.imshow(np.clip(_overlay, 0, 1), interpolation="nearest")
    #ax_o.set_title("green = fitted  |  red = excluded", fontsize=8)
    fig_o.tight_layout(pad=0.3)
    _overlay_path = output_dir / f"{_stem_early}_otsu_mask_overlay.png"
    fig_o.savefig(_overlay_path, dpi=200, bbox_inches="tight")
    plt.close(fig_o)
    print(f"  Saved: {_overlay_path.name}")

    # ── 6. Global fit window ──────────────────────────────────────────────
    start_idx, end_idx, start_ns, end_ns = compute_global_fit_window(
        stack, mask, time_ns, peak_offset_bins=peak_offset_bins
    )
    # Override with user-specified bounds when provided
    if fit_start_ns is not None:
        start_idx = int(np.clip(np.searchsorted(time_ns, fit_start_ns),
                                0, len(time_ns) - 2))
        start_ns  = float(time_ns[start_idx])
    if fit_end_ns is not None:
        end_idx = int(np.clip(np.searchsorted(time_ns, fit_end_ns, side="right"),
                              start_idx + 2, len(time_ns)))
        end_ns  = float(time_ns[end_idx - 1])
    src = ("user-specified"
           if (fit_start_ns is not None or fit_end_ns is not None)
           else f"IRF offset = {peak_offset_bins} bins = {peak_offset_bins * dt:.3f} ns")
    print(
        f"  Fit window  : [{start_ns:.3f}, {end_ns:.3f}] ns  "
        f"(bins {start_idx}–{end_idx - 1}, {src})"
    )

    # ── 7. Per-pixel fitting ──────────────────────────────────────────────
    maps = compute_lifetime_maps(
        stack, mask, time_ns, start_idx, end_idx,
        bin_radius=bin_radius,
        min_counts=min_counts,
        fit_biexp=fit_biexp,
        tau_max_ns=tau_max_ns,
    )
    n_valid_mono = int(np.isfinite(maps["tau_mono"]).sum())
    n_valid_bi   = int(np.isfinite(maps["tau_avg"]).sum())
    print(
        f"  Valid fits  : mono {n_valid_mono}/{n_mask}  |  "
        f"bi-exp {n_valid_bi}/{n_mask}"
    )

    # ── 8. Colour range (2nd–98th percentile per map) ─────────────────────
    vm_mono = _percentile_range(maps["tau_mono"])
    vm_bi   = _percentile_range(maps["tau_avg"])
    # For τ₁/τ₂ exclude boundary-hitting pixels from the colour range so that
    # a few degenerate fits at 12 ns do not compress the entire colourmap.
    _t1_ok = maps["tau1"][np.isfinite(maps["tau1"]) &
                           (maps["tau1"] < _TAU_MAX_NS - _TAU_BOUND_TOL)]
    vm1 = _percentile_range(_t1_ok) if _t1_ok.size > 0 else (0.0, _TAU_MAX_NS)
    _t2_ok = maps["tau2"][np.isfinite(maps["tau2"]) &
                           (maps["tau2"] < _TAU_MAX_NS - _TAU_BOUND_TOL)]
    vm2 = _percentile_range(_t2_ok) if _t2_ok.size > 0 else (0.0, _TAU_MAX_NS)

    # ── 9. Statistics ─────────────────────────────────────────────────────
    def _stats(arr: np.ndarray, bound_tol: float = 0.0) -> dict:
        v = arr[np.isfinite(arr)]
        n_total = int(v.size)
        if bound_tol > 0.0:
            v = v[
                (v > _TAU_MIN_NS + bound_tol) &
                (v < _TAU_MAX_NS - bound_tol)
            ]
        n_excl = n_total - int(v.size)
        if v.size == 0:
            return {"mean": np.nan, "median": np.nan, "std": np.nan,
                    "mad": np.nan, "p05": np.nan,
                    "p95": np.nan, "n_excl": n_excl}
        _med = float(np.median(v))
        return {
            "mean":   float(v.mean()),
            "median": _med,
            "std":    float(v.std()),
            "mad":    float(np.median(np.abs(v - _med))),
            "p05":    float(np.percentile(v, 5)),
            "p95":    float(np.percentile(v, 95)),
            "n_excl": n_excl,
        }

    sm = _stats(maps["tau_mono"])
    sb = _stats(maps["tau_avg"])
    s1 = _stats(maps["tau1"], bound_tol=_TAU_BOUND_TOL)
    s2 = _stats(maps["tau2"], bound_tol=_TAU_BOUND_TOL)
    print(
        f"  τ_mono      : {sm['mean']:.3f} ± {sm['std']:.3f} ns  "
        f"median = {sm['median']:.3f}  [5–95 %: {sm['p05']:.3f} – {sm['p95']:.3f}]"
    )
    print(
        f"  ⟨τ⟩_bi-exp  : {sb['mean']:.3f} ± {sb['std']:.3f} ns  "
        f"median = {sb['median']:.3f}  [5–95 %: {sb['p05']:.3f} – {sb['p95']:.3f}]"
    )
    if fit_biexp:
        excl1 = s1.get("n_excl", 0)
        excl2 = s2.get("n_excl", 0)
        excl1_str = f"  [{excl1} px at bound excl.]" if excl1 > 0 else ""
        excl2_str = f"  [{excl2} px at bound excl.]" if excl2 > 0 else ""
        print(
            f"  τ₁          : {s1['mean']:.3f} ± {s1['std']:.3f} ns  "
            f"median = {s1['median']:.3f}  [5–95 %: {s1['p05']:.3f} – {s1['p95']:.3f}]{excl1_str}"
        )
        print(
            f"  τ₂          : {s2['mean']:.3f} ± {s2['std']:.3f} ns  "
            f"median = {s2['median']:.3f}  MAD = {s2['mad']:.3f}  "
            f"[5–95 %: {s2['p05']:.3f} – {s2['p95']:.3f}]{excl2_str}"
        )

    # ── 10. Save numpy arrays (for cross-section comparison) ─────────────
    output_dir.mkdir(parents=True, exist_ok=True)
    stem          = tiff_path.stem.replace(" ", "_")
    fit_range_str = f"{start_ns:.2f}-{end_ns:.2f}ns"

    npy_keys = ["tau_mono", "tau_avg", "tau1", "tau2"] if fit_biexp else ["tau_mono"]
    for key in npy_keys:
        np.save(output_dir / f"{stem}_{key}.npy", maps[key])

    # ── 11. Save plots ────────────────────────────────────────────────────

    # (map_key, filename_tag, cmap_name, (vmin, vmax))
    plot_specs = [
        ("tau_mono", "tau_mono",       colormap, vm_mono),
        ("tau_avg",  "tau_biexp_mean", colormap, vm_bi),
        ("tau1",     "tau1",           colormap, vm1),
        ("tau2",     "tau2",           colormap, vm2),
    ]

    for map_key, name_tag, cmap, (vmin, vmax) in plot_specs:
        if not fit_biexp and map_key in ("tau_avg", "tau1", "tau2"):
            continue
        cmap_tag = "gray" if cmap == "gray" else "color"
        fname    = f"{stem}_{name_tag}_{fit_range_str}_{cmap_tag}.png"
        save_lifetime_plot(
            intensity_proj,
            maps[map_key],
            output_dir / fname,
            colormap=cmap,
            vmin=vmin,
            vmax=vmax,
        )

    # ── 11b. Lifetime histograms (one PNG per quantity) ───────────────────
    save_tau_histogram(
        maps["tau_mono"], sm,
        "τ (mono-exp fit)",
        output_dir / f"{stem}_hist_tau_mono_{fit_range_str}.png",
        color="#e25822",
    )
    if fit_biexp and np.isfinite(maps["tau_avg"]).any():
        save_tau_histogram(
            maps["tau_avg"], sb,
            "⟨τ⟩ (bi-exp amplitude-weighted mean)",
            output_dir / f"{stem}_hist_tau_avg_{fit_range_str}.png",
            color="#5577cc",
        )
        save_tau_histogram(
            maps["tau1"], s1,
            "τ₁ (bi-exponential component 1)",
            output_dir / f"{stem}_hist_tau1_{fit_range_str}.png",
            color="#2ca02c",
            exclude_bounds_tol=_TAU_BOUND_TOL,
        )
        save_tau_histogram(
            maps["tau2"], s2,
            "τ₂ (bi-exponential component 2)",
            output_dir / f"{stem}_hist_tau2_{fit_range_str}.png",
            color="#9467bd",
            exclude_bounds_tol=_TAU_BOUND_TOL,
        )

    # ── 12. Summary row ───────────────────────────────────────────────────
    return {
        "file":               str(tiff_path),
        "n_bins":             N_bins,
        "time_bin_width_ns":  dt,
        "fit_start_ns":       start_ns,
        "fit_end_ns":         end_ns,
        "peak_offset_bins":   peak_offset_bins,
        "bin_radius":         bin_radius,
        "n_mask_pixels":      n_mask,
        "n_valid_mono":       n_valid_mono,
        "n_valid_biexp":      n_valid_bi,
        "tau_mono_mean_ns":   sm["mean"],
        "tau_mono_median_ns": sm["median"],
        "tau_mono_std_ns":    sm["std"],
        "tau_mono_p05_ns":    sm["p05"],
        "tau_mono_p95_ns":    sm["p95"],
        "tau_avg_mean_ns":    sb["mean"],
        "tau_avg_median_ns":  sb["median"],
        "tau_avg_std_ns":     sb["std"],
        "tau_avg_p05_ns":     sb["p05"],
        "tau_avg_p95_ns":     sb["p95"],
        "tau1_mean_ns":        s1["mean"],
        "tau1_median_ns":      s1["median"],
        "tau1_std_ns":         s1["std"],
        "tau1_mad_ns":         s1["mad"],
        "tau2_mean_ns":        s2["mean"],
        "tau2_median_ns":      s2["median"],
        "tau2_std_ns":         s2["std"],
        "tau2_mad_ns":         s2["mad"],
        "otsu_threshold":     float(otsu_thr),
        "artifact_detected":  bool(trim_report["detected"]),
    }


# =============================================================================
# CLI
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build per-pixel FLIM lifetime maps from TCSPC TIFF stacks.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "input", type=str,
        help="Path to a TIFF stack, or a folder when --batch is used",
    )
    parser.add_argument(
        "--batch", action="store_true",
        help="Process every .tif/.tiff inside the input folder",
    )
    parser.add_argument(
        "--output-dir", type=str, default=None,
        help=f"Output folder (default: {_DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--laser-rate", type=float, default=DEFAULT_LASER_REP_RATE_MHZ,
        help=f"Laser repetition rate in MHz (default: {DEFAULT_LASER_REP_RATE_MHZ})",
    )
    parser.add_argument(
        "--peak-offset-bins", type=int, default=7,
        help="Bins to skip to the right of the global peak to avoid the IRF "
             "rising edge (default: 7)",
    )
    parser.add_argument(
        "--min-counts", type=int, default=30,
        help="Min photon count in the fit window before attempting a fit "
             "(default: 30)",
    )
    parser.add_argument(
        "--bin-radius", type=int, default=1,
        help="Rolling-window half-width: 0 → 1×1 (single pixel), 1 → 3×3, 2 → 5×5 (default: 1)",
    )
    parser.add_argument(
        "--colormap", type=str, default="inferno",
        help="Matplotlib colormap for the colour plots (default: inferno)",
    )
    parser.add_argument(
        "--hot-pixel-k", type=float, default=5.0,
        help="Hot-pixel detection threshold factor k (default: 5.0)",
    )
    parser.add_argument(
        "--artifact-spike-factor", type=float, default=2.5,
        help="Trailing-artefact spike sensitivity factor (default: 2.5)",
    )
    parser.add_argument(
        "--no-biexp", action="store_true",
        help="Skip bi-exponential fitting (produces only τ_mono maps; faster)",
    )
    parser.add_argument(
        "--fit-start-ns", type=float, default=None,
        help="Fit window start in ns (default: auto from global peak + offset)",
    )
    parser.add_argument(
        "--fit-end-ns", type=float, default=None,
        help="Fit window end in ns (default: auto from last bin with signal)",
    )
    parser.add_argument(
        "--tau-max-ns", type=float, default=12.0,
        help="Upper bound on τ for bi-exponential fitting in ns (default: 12.0). "
             "Tighten to reduce outliers in τ₂, e.g. --tau-max-ns 6 for DAPI in cells",
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir) if args.output_dir else _DEFAULT_OUTPUT_DIR

    common_kw = dict(
        output_dir=output_dir,
        peak_offset_bins=args.peak_offset_bins,
        min_counts=args.min_counts,
        bin_radius=args.bin_radius,
        laser_rate=args.laser_rate,
        colormap=args.colormap,
        hot_pixel_k=args.hot_pixel_k,
        artifact_spike_factor=args.artifact_spike_factor,
        fit_biexp=not args.no_biexp,
        fit_start_ns=args.fit_start_ns,
        fit_end_ns=args.fit_end_ns,
        tau_max_ns=args.tau_max_ns,
    )

    input_path = Path(args.input)

    if args.batch:
        if not input_path.is_dir():
            parser.error(f"--batch requires a folder, got: {input_path}")
        tiff_files = sorted(
            list(input_path.glob("*.tif")) + list(input_path.glob("*.tiff"))
        )
        if not tiff_files:
            parser.error(f"No .tif/.tiff files found in {input_path}")
        print(f"Batch mode: {len(tiff_files)} TIFF file(s) in {input_path}")
        rows = []
        for tp in tiff_files:
            try:
                rows.append(process_tiff_lifetime_map(tp, **common_kw))
            except Exception as exc:
                print(f"\n  !! FAILED {tp.name}: {exc}")
        if rows:
            summary_path = output_dir / "lifetime_map_summary.csv"
            pd.DataFrame(rows).to_csv(summary_path, index=False)
            print(f"\nSummary saved: {summary_path}")

    else:
        if not input_path.is_file():
            parser.error(f"File not found: {input_path}")
        row = process_tiff_lifetime_map(input_path, **common_kw)
        summary_path = output_dir / "lifetime_map_summary.csv"
        pd.DataFrame([row]).to_csv(summary_path, index=False)
        print(f"\nSummary saved: {summary_path}")

    print("\nDone!")


if __name__ == "__main__":
    main()
