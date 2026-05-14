"""
FLIM TIFF Stack Reader and Fitter
==================================

This script reads a FLIM intensity TIFF stack (e.g. 256 frames of 128x128
pixels) where each frame corresponds to one TCSPC time bin, builds a 1-D
intensity-vs-time decay by aggregating the photon counts over all pixels
(or over an optional ROI), and fits it using the same mono-/bi-exponential
machinery already implemented in ``flim_exponential_fit.py``.

The TIFF stack is expected to have shape (N_bins, H, W). For an 80 MHz
laser the laser period is 12.5 ns, so with N_bins = 256 the time bin width
is 12.5 / 256 ~ 0.0488 ns. The number of bins is auto-detected from the
stack length, so it works equally well with 128, 256, ... bins.

Trailing-artifact removal
--------------------------
Some TCSPC software re-injects a delayed copy of the laser pulse at the end
of the acquisition window.  This produces a sharp intensity spike in the last
few non-zero frames that has no physical meaning and would bias the lifetime
fit.  The script detects and removes these frames automatically before any
other processing (trimming is ON by default, use ``--no-trim`` to disable).
The sensitivity is controlled by ``--artifact-spike-factor`` (default 2.5):
a frame is flagged as artifact when its mean intensity exceeds 2.5 times the
local baseline estimated from the 30 frames just before the spike region.

Denoising pipeline
-------------------
When ``--denoise`` is active, before collapsing the stack to a 1-D trace the
script applies:
  1. Hot-pixel removal: pixels with value > k × local 3×3 median are replaced
     by that median (``--hot-pixel-k``, default 5).
  2. Intensity mask: only pixels above the Otsu threshold on the time-projected
     image contribute to the decay trace (requires scikit-image).
     The denoised & masked 3-D stack (N_bins × H × W) is saved as a multi-page
     TIFF (``<name>_stack_denoised.tif``) for inspection in ImageJ / Fiji.

Quick start
-----------
    # Single file, both models, Otsu denoising, no GUI (recommended default)
    python flim_tiff_fit.py "altro/morti decadimento immagine 1s_Intensity.tif" \\
        --denoise --mask otsu --no-show

    # All 'morti/vivi decadimento' files → Batteri folder
    foreach ($f in @(
        "morti decadimento immagine 1s_Intensity.tif",
        "morti decadimento immagine 2s_Intensity.tif",
        "vivi decadimento immagine 1s_Intensity.tif",
        "vivi decadimento immagine 2_Intensity.tif",
        "vivi decadimento immagine 2s_Intensity.tif",
        "vivi decadimento immagine 3_Intensity.tif"
    )) {
        python flim_tiff_fit.py "altro/$f" --denoise --mask otsu \\
            --output-dir img/Batteri --no-show
    }

    # All 'Morti Decay Raw Data' files → Batteri_2 folder
    foreach ($f in @(
        "Morti Decay Raw Data 1_Intensity.tif",
        "Morti Decay Raw Data 3_Intensity.tif",
        "Morti Decay Raw Data 4_Intensity.tif"
    )) {
        python flim_tiff_fit.py "altro/$f" --denoise --mask otsu \\
            --output-dir img/Batteri_2 --no-show
    }

    # Fit a single TIFF stack (mono-exponential only)
    python flim_tiff_fit.py "altro/vivi decadimento immagine 1s_Intensity.tif" --model mono

    # Bi-exponential, sum pixels instead of averaging, no GUI
    python flim_tiff_fit.py "altro/morti decadimento immagine 1s_Intensity.tif" \\
        --model bi --reduce sum --no-show

    # Batch process every .tif in a folder (no denoising)
    python flim_tiff_fit.py altro --batch --model bi --no-show

    # Disable trailing-artifact trimming (not recommended)
    python flim_tiff_fit.py "altro/morti decadimento immagine 1s_Intensity.tif" --no-trim

    # Tune detection sensitivity (lower = more aggressive)
    python flim_tiff_fit.py "altro/morti decadimento immagine 1s_Intensity.tif" \\
        --artifact-spike-factor 2.0

Key options
-----------
    --model {mono,bi,both}          Decay model (default: both)
    --denoise                       Enable hot-pixel removal + Otsu mask
    --mask {otsu,percentile,        Masking method when --denoise is used
            manual,none}            (Otsu threshold; scikit-image required)
    --hot-pixel-k FLOAT             Hot-pixel sensitivity (default: 5.0)
    --smooth-sigma FLOAT            Gaussian pre-blur sigma in px (default: 0)
    --no-trim                       Disable trailing-artifact removal
    --artifact-spike-factor FLOAT   Spike detection sensitivity (default: 2.5)
    --reduce {mean,sum}             Pixel aggregation per frame (default: mean)
    --roi Y0,Y1,X0,X1               Restrict analysis to a rectangular ROI
    --fit-start FLOAT               Manual fit start in ns (default: auto peak)
    --fit-end FLOAT                 Manual fit end in ns (default: last non-zero)
    --output-dir PATH               Output folder (default: img/Batteri)
    --output-name NAME              Custom filename prefix (single-file mode)
    --no-show                       Save plots without opening the window
    --batch                         Process every .tif in the input folder
    --laser-rate FLOAT              Laser repetition rate in MHz (default: 80)
"""

from __future__ import annotations

import argparse
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi

try:
    from skimage.filters import threshold_otsu
    from skimage.morphology import remove_small_objects, dilation as sk_dilation, disk
    _HAS_SKIMAGE = True
except ImportError:  # pragma: no cover - fallback when skimage missing
    threshold_otsu = None
    remove_small_objects = None
    binary_dilation = None
    _HAS_SKIMAGE = False

from flim_exponential_fit import (
    FLIMData,
    DEFAULT_LASER_REP_RATE_MHZ,
    extract_decay_region,
    get_data_summary,
    fit_mono_exponential,
    fit_bi_exponential,
    plot_mono_exponential_fit,
    plot_bi_exponential_fit,
    print_mono_result,
    print_bi_result,
    resolve_fit_indices,
)


# Default output directory for TIFF-based FLIM analysis
OUTPUT_DIR = Path(__file__).parent / "img" / "Batteri"


# =============================================================================
# TIFF reading
# =============================================================================

def _aggregate_stack(
    stack: np.ndarray,
    reduce: str = "mean",
    roi: Optional[Tuple[int, int, int, int]] = None,
) -> np.ndarray:
    """
    Collapse a (N_bins, H, W) stack into a 1-D intensity trace of length N_bins.

    Parameters
    ----------
    stack : np.ndarray
        TIFF stack with shape (N_bins, H, W).
    reduce : {'mean', 'sum'}
        How to aggregate pixels within each frame. 'mean' matches what
        ImageJ's "Plot Z-axis profile" produces (mean grey value per slice);
        'sum' gives the total photon count per time bin.
    roi : (y0, y1, x0, x1), optional
        Optional rectangular region of interest (slice end-exclusive).

    Returns
    -------
    np.ndarray
        1-D float array of length N_bins.
    """
    if stack.ndim != 3:
        raise ValueError(
            f"Expected a 3-D TIFF stack with shape (N_bins, H, W), got shape {stack.shape}"
        )

    if roi is not None:
        y0, y1, x0, x1 = roi
        stack = stack[:, y0:y1, x0:x1]

    arr = stack.astype(np.float64, copy=False)

    if reduce == "mean":
        return arr.mean(axis=(1, 2))
    if reduce == "sum":
        return arr.sum(axis=(1, 2))
    raise ValueError(f"Unknown reduce mode: {reduce!r} (use 'mean' or 'sum')")


# =============================================================================
# Lightweight denoising (hot-pixel removal + intensity mask)
# =============================================================================

@dataclass
class FLIMDenoiseArtifacts:
    """Side products of the optional denoising stage (for saving / inspection)."""
    intensity_raw: np.ndarray          # 2-D, sum of the raw stack along time
    intensity_denoised: np.ndarray     # 2-D, after hot-pixel removal (+ smoothing)
    mask: np.ndarray                   # 2-D bool, True = kept pixels (bacteria)
    intensity_masked: np.ndarray       # 2-D float, denoised image * mask
    stack_denoised: np.ndarray         # 3-D (N_bins, H, W), hot-pixel cleaned + masked
    n_pixels_in_mask: int
    mask_threshold: float
    hot_pixel_k: float


def _remove_hot_pixels_2d(img: np.ndarray, k: float = 5.0) -> np.ndarray:
    """
    Replace isolated hot pixels by the local 3x3 median.

    A pixel is flagged as 'hot' when its value exceeds ``k`` times the local
    median (computed in a 3x3 neighbourhood). This kills the bright single-pixel
    speckles that dominate the background of these 2PE FLIM images without
    touching extended bacterial clusters.
    """
    if img.size == 0:
        return img
    med = ndi.median_filter(img, size=3, mode="reflect")
    # Avoid div-by-zero on background regions (median == 0)
    safe_med = np.where(med > 0, med, 1.0)
    hot = img > k * safe_med
    out = img.copy()
    out[hot] = med[hot]
    return out


def _remove_hot_pixels_stack(stack: np.ndarray, k: float = 5.0) -> np.ndarray:
    """Apply ``_remove_hot_pixels_2d`` frame-by-frame on a (N, H, W) stack."""
    out = np.empty_like(stack, dtype=np.float64)
    for i in range(stack.shape[0]):
        out[i] = _remove_hot_pixels_2d(stack[i].astype(np.float64, copy=False), k=k)
    return out


def _build_mask(
    img: np.ndarray,
    min_object_size: int = 5,
    dilate: bool = True,
) -> Tuple[np.ndarray, float]:
    """
    Build a boolean mask of 'bacterial' pixels using Otsu thresholding.

    Requires scikit-image. Applies morphological cleanup to remove tiny
    isolated bright pixels and dilates the mask by 1 pixel.
    """
    if not _HAS_SKIMAGE:
        raise RuntimeError(
            "scikit-image is required for Otsu masking. "
            "Install it with: pip install scikit-image"
        )

    thr = float(threshold_otsu(img))
    mask = img > thr

    # Morphological cleanup: drop tiny isolated bright pixels (background dots).
    if min_object_size > 1:
        mask = remove_small_objects(mask, max_size=int(min_object_size) - 1)
        if dilate:
            mask = sk_dilation(mask, disk(1))

    return mask.astype(bool), thr


def _aggregate_stack_denoised(
    stack: np.ndarray,
    reduce: str = "mean",
    roi: Optional[Tuple[int, int, int, int]] = None,
    hot_pixel_k: float = 5.0
) -> Tuple[np.ndarray, FLIMDenoiseArtifacts]:
    """
    Hot-pixel cleanup + intensity-mask + temporal aggregation.

    Returns the 1-D decay trace AND a ``FLIMDenoiseArtifacts`` bundle with the
    intermediate 2-D images used to build the mask (handy for saving/QC).

    With ``reduce='mean'``, the average is taken over the masked pixels only,
    so the trace is interpretable as 'mean intensity per bacterial pixel'.
    """
    if stack.ndim != 3:
        raise ValueError(
            f"Expected a 3-D TIFF stack with shape (N_bins, H, W), got shape {stack.shape}"
        )

    if roi is not None:
        y0, y1, x0, x1 = roi
        stack = stack[:, y0:y1, x0:x1]

    arr = stack.astype(np.float64, copy=False)

    # 1) projected intensity image (raw)
    intensity_raw = arr.sum(axis=0)

    # 2) hot-pixel cleanup of the projected image (used to build the mask)
    intensity_denoised = _remove_hot_pixels_2d(intensity_raw, k=hot_pixel_k)

    # 3) mask from the cleaned projected image
    mask, thr = _build_mask(intensity_denoised)

    # 4) clean the time-resolved stack frame-by-frame and apply the mask
    stack_clean = _remove_hot_pixels_stack(arr, k=hot_pixel_k)
    masked_stack = stack_clean * mask[None, :, :]

    n_pixels = int(mask.sum())
    if n_pixels == 0:
        raise RuntimeError(
            "Denoising mask is empty: no pixels survived the Otsu threshold. "
            "Check image quality or disable denoising with --no-denoise."
        )

    if reduce == "mean":
        trace = masked_stack.sum(axis=(1, 2)) / n_pixels
    elif reduce == "sum":
        trace = masked_stack.sum(axis=(1, 2))
    else:
        raise ValueError(f"Unknown reduce mode: {reduce!r} (use 'mean' or 'sum')")

    artifacts = FLIMDenoiseArtifacts(
        intensity_raw=intensity_raw,
        intensity_denoised=intensity_denoised,
        mask=mask,
        intensity_masked=intensity_denoised * mask,
        stack_denoised=masked_stack.astype(np.float32),
        n_pixels_in_mask=n_pixels,
        mask_threshold=thr,
        hot_pixel_k=hot_pixel_k,
    )
    return trace, artifacts


def save_denoised_stack(artifacts: FLIMDenoiseArtifacts, output_path: Path) -> None:
    """
    Save the hot-pixel cleaned and masked 3-D stack to a multi-page TIFF.

    The output has shape (N_bins, H, W) with float32 pixel values and can be
    opened directly in ImageJ / Fiji as a hyperstack.
    """
    tifffile.imwrite(
        str(output_path),
        artifacts.stack_denoised,          # (N_bins, H, W) float32
        photometric="minisblack",
    )


# =============================================================================
# Trailing-artifact detection and trimming
# =============================================================================

def detect_trailing_artifact_cutoff(
    mean_per_frame: np.ndarray,
    eps_frac: float = 1e-4,
    spike_factor: float = 2.5,
    look_back: int = 30,
) -> Tuple[int, dict]:
    """
    Detect the index at which trailing artifact frames begin.

    Some TCSPC acquisition software re-injects a delayed copy of the laser
    pulse at the very end of the time window, producing a sharp intensity
    spike in the last few non-zero frames that has no physical meaning and
    would bias any lifetime fit.

    Strategy
    --------
    1. Identify trailing near-zero frames (detector dark counts).
    2. Estimate a local baseline as the 20th percentile of the
       ``look_back`` frames just before the suspected spike (robust to the
       spike itself being inside the window).
    3. Walk backwards from the last non-zero frame; flag every frame whose
       mean exceeds ``spike_factor * baseline`` as an artifact.
    4. Return the index of the first artifact frame (= the cut point).

    Parameters
    ----------
    mean_per_frame : 1-D array
        Mean pixel intensity for each temporal frame of the stack.
    eps_frac : float
        A frame is considered 'zero' if mean < eps_frac * global_max.
    spike_factor : float
        Ratio above the local baseline that flags a frame as artifact.
    look_back : int
        Size of the window (in frames) used to estimate the baseline.

    Returns
    -------
    cutoff_idx : int
        First frame index to remove.  Use ``stack[:cutoff_idx]`` to keep
        only good data.  Returns ``len(mean_per_frame)`` when no artifact
        is detected.
    report : dict
        Diagnostic fields for logging / CSV export.
    """
    n = len(mean_per_frame)
    global_max = float(mean_per_frame.max())
    eps = global_max * eps_frac

    # Step 1 – find last non-zero frame
    last_nz = n - 1
    while last_nz >= 0 and mean_per_frame[last_nz] <= eps:
        last_nz -= 1

    if last_nz < 0:
        return 0, {"detected": False, "reason": "all-zero stack",
                   "cutoff_idx": 0, "n_removed": n,
                   "n_trailing_zeros": n, "baseline": 0.0, "threshold": 0.0}

    n_trailing_zeros = n - 1 - last_nz

    # Step 2 – estimate local baseline just before the potential spike
    win_end   = last_nz + 1
    win_start = max(0, win_end - look_back)
    baseline  = float(np.percentile(mean_per_frame[win_start:win_end], 20))

    if baseline <= 0.0:
        return n, {"detected": False, "reason": "zero baseline",
                   "cutoff_idx": n, "n_removed": 0,
                   "n_trailing_zeros": n_trailing_zeros,
                   "baseline": 0.0, "threshold": 0.0}

    threshold = baseline * spike_factor

    # Step 3 – walk backwards to find start of contiguous artifact block
    spike_start = n  # sentinel: no spike found
    for i in range(last_nz, -1, -1):
        if mean_per_frame[i] > threshold:
            spike_start = i
        else:
            break  # first normal frame; stop

    if spike_start == n:
        return n, {
            "detected": False, "cutoff_idx": n, "n_removed": 0,
            "n_trailing_zeros": n_trailing_zeros,
            "baseline": baseline, "threshold": threshold,
        }

    n_artifact = last_nz + 1 - spike_start   # non-zero artifact frames
    n_removed  = n - spike_start             # total frames trimmed (artifact + zeros)

    return spike_start, {
        "detected": True,
        "cutoff_idx": spike_start,
        "n_removed": n_removed,
        "n_artifact_frames": n_artifact,
        "n_trailing_zeros": n_trailing_zeros,
        "baseline": baseline,
        "threshold": threshold,
        "spike_mean_max": float(mean_per_frame[spike_start: last_nz + 1].max()),
        "spike_factor_observed": float(
            mean_per_frame[spike_start: last_nz + 1].max() / baseline
        ),
    }


def _print_artifact_report(tiff_path: Path, report: dict) -> None:
    """Print a human-readable summary of the artifact-detection result."""
    sep = "-" * 60
    print(sep)
    print(f"Artifact scan: {tiff_path.name}")
    if not report["detected"]:
        reason = report.get("reason", "intensity within normal range")
        print(f"  No trailing artifact detected ({reason}).")
        if report.get("n_trailing_zeros", 0):
            print(f"  Trailing zero frames: {report['n_trailing_zeros']} (left intact).")
    else:
        print(f"  *** Trailing artifact DETECTED ***")
        print(f"  First artifact frame : {report['cutoff_idx']} "
              f"(1-indexed: {report['cutoff_idx'] + 1})")
        print(f"  Artifact frames      : {report['n_artifact_frames']} "
              f"(peak mean = {report['spike_mean_max']:.4g}, "
              f"{report['spike_factor_observed']:.1f}x baseline)")
        print(f"  Trailing zero frames : {report['n_trailing_zeros']}")
        print(f"  Total frames removed : {report['n_removed']} "
              f"(stack trimmed to {report['cutoff_idx']} frames)")
        print(f"  Baseline / threshold : {report['baseline']:.4g} / "
              f"{report['threshold']:.4g}")
    print(sep)


def read_flim_tiff(
    tiff_path: str | Path,
    laser_rep_rate_mhz: float = DEFAULT_LASER_REP_RATE_MHZ,
    reduce: str = "mean",
    roi: Optional[Tuple[int, int, int, int]] = None,
    denoise: bool = True,
    hot_pixel_k: float = 5.0,
    return_artifacts: bool = False,
    trim_artifacts: bool = True,
    artifact_spike_factor: float = 2.5,
):
    """
    Read a FLIM intensity TIFF stack and convert it to a FLIMData object.

    Each frame in the stack is interpreted as one TCSPC time bin. The stack
    is collapsed to a 1-D intensity trace by averaging (or summing) pixels.

    When ``trim_artifacts=True`` (default), the stack is scanned for trailing
    frames with anomalously high mean intensity before any other processing.
    These are a known TCSPC software artefact where the laser-triggered signal
    is re-injected at the end of the time window; they have no physical meaning
    and would bias the lifetime fit.  Detected frames are removed and a
    diagnostic report is printed.  Use ``artifact_spike_factor`` to tune
    sensitivity (default 3.0 = frame must be >3x the local background).

    When ``denoise=True``, a light-touch cleanup is applied before aggregation:
    isolated hot pixels are replaced by the local 3x3 median (controlled by
    ``hot_pixel_k``), and only the pixels that pass the Otsu mask contribute to
    the decay trace. With ``reduce='mean'`` the average is then taken over the
    masked pixels only.

    If ``return_artifacts=True``, returns ``(FLIMData, FLIMDenoiseArtifacts |
    None)`` so callers can save the projected image / mask alongside the fit.
    """
    tiff_path = Path(tiff_path)
    if not tiff_path.exists():
        raise FileNotFoundError(f"TIFF file not found: {tiff_path}")

    stack = tifffile.imread(str(tiff_path))

    # --- Trailing-artifact detection and trimming (before any other processing) ---
    trim_report: dict = {"detected": False, "cutoff_idx": stack.shape[0] if stack.ndim == 3 else 0, "n_removed": 0}
    if trim_artifacts and stack.ndim == 3:
        mean_per_frame = stack.astype(np.float64).mean(axis=(1, 2))
        cutoff_idx, trim_report = detect_trailing_artifact_cutoff(
            mean_per_frame, spike_factor=artifact_spike_factor
        )
        _print_artifact_report(tiff_path, trim_report)
        if trim_report["detected"]:
            stack = stack[:cutoff_idx]

    artifacts: Optional[FLIMDenoiseArtifacts] = None
    if denoise:
        intensity, artifacts = _aggregate_stack_denoised(
            stack, reduce=reduce, roi=roi,
            hot_pixel_k=hot_pixel_k,
        )
    else:
        intensity = _aggregate_stack(stack, reduce=reduce, roi=roi)

    num_bins = int(intensity.size)
    time_range_ns = 1000.0 / laser_rep_rate_mhz       # MHz -> ns period
    time_bin_width_ns = time_range_ns / num_bins
    time_ns = np.arange(num_bins, dtype=float) * time_bin_width_ns

    peak_index = int(np.argmax(intensity))
    peak_time_ns = float(time_ns[peak_index])

    data = FLIMData(
        time_ns=time_ns,
        intensity=intensity,
        time_bin_width_ns=time_bin_width_ns,
        num_bins=num_bins,
        peak_index=peak_index,
        peak_time_ns=peak_time_ns,
    )
    if return_artifacts:
        return data, artifacts, trim_report
    return data


# =============================================================================
# Driver: read TIFF + fit + plot (re-uses flim_exponential_fit functions)
# =============================================================================


def _parse_roi(spec: Optional[str]) -> Optional[Tuple[int, int, int, int]]:
    if spec is None:
        return None
    parts = [int(p) for p in spec.split(",")]
    if len(parts) != 4:
        raise argparse.ArgumentTypeError(
            "ROI must be 'y0,y1,x0,x1' (got: %r)" % spec
        )
    return tuple(parts)  # type: ignore[return-value]


def process_tiff(
    tiff_path: Path,
    output_dir: Path,
    model: str,
    laser_rate: float,
    reduce: str,
    roi: Optional[Tuple[int, int, int, int]],
    fit_start_ns: Optional[float],
    fit_end_ns: Optional[float],
    show_plot: bool,
    output_name: Optional[str] = None,
    denoise: bool = True,
    hot_pixel_k: float = 5.0,
    trim_artifacts: bool = True,
    artifact_spike_factor: float = 2.5,
) -> dict:
    """Read one TIFF, run fit(s), save outputs, and return extracted lifetimes."""
    print(f"\nReading FLIM TIFF stack: {tiff_path}")
    print(f"Laser repetition rate: {laser_rate} MHz   |   pixel reduce: {reduce}")
    if roi is not None:
        print(f"ROI: y={roi[0]}:{roi[1]}, x={roi[2]}:{roi[3]}")
    if not trim_artifacts:
        print("Trailing-artifact trimming: DISABLED (--no-trim)")
    if denoise:
        print(f"Denoise: ON (hot_pixel_k={hot_pixel_k}, mask=otsu)")

    data, artifacts, trim_report = read_flim_tiff(
        tiff_path, laser_rep_rate_mhz=laser_rate, reduce=reduce, roi=roi,
        denoise=denoise, hot_pixel_k=hot_pixel_k,
        return_artifacts=True,
        trim_artifacts=trim_artifacts,
        artifact_spike_factor=artifact_spike_factor,
    )
    print(get_data_summary(data))
    if artifacts is not None:
        print(f"Mask: {artifacts.n_pixels_in_mask} pixels kept "
              f"(method=otsu, threshold={artifacts.mask_threshold:.3g})")

    output_dir.mkdir(parents=True, exist_ok=True)
    input_stem_clean = tiff_path.stem.replace(" ", "_")
    suffix = "_denoised" if denoise else ""

    # Save artifacts of the denoising stage (raw vs cleaned image, mask)
    if artifacts is not None:
        raw_png = output_dir / f"{input_stem_clean}_intensity_raw.png"
        den_png = output_dir / f"{input_stem_clean}_intensity_denoised.png"
        msk_png = output_dir / f"{input_stem_clean}_mask.png"
        msk_tif = output_dir / f"{input_stem_clean}_intensity_masked.tif"
        plt.imsave(raw_png, artifacts.intensity_raw, cmap="gray")
        plt.imsave(den_png, artifacts.intensity_denoised, cmap="gray")
        plt.imsave(msk_png, artifacts.mask.astype(np.uint8), cmap="gray", vmin=0, vmax=1)
        tifffile.imwrite(str(msk_tif), artifacts.intensity_masked.astype(np.float32))
        stk_tif = output_dir / f"{input_stem_clean}_stack_denoised.tif"
        save_denoised_stack(artifacts, stk_tif)
        print(f"Saved denoising artifacts: {raw_png.name}, {den_png.name}, "
              f"{msk_png.name}, {msk_tif.name}, {stk_tif.name}")

    # Save the 1-D time/intensity trace for reference and repucibility
    csv_output_path = output_dir / f"{input_stem_clean}{suffix}_time_intensity.csv"
    pd.DataFrame({"time_ns": data.time_ns,
                  "intensity": data.intensity}).to_csv(csv_output_path, index=False)
    print(f"Time-intensity trace saved to: {csv_output_path}")

    # Decide the fit window using the same logic as the CSV-based script.
    # min_counts=5 excludes the low-statistics tail (valid only for reduce='sum')
    _min_counts = 5 if reduce == "sum" else 1
    start_idx, end_idx, auto_used = resolve_fit_indices(data, fit_start_ns, fit_end_ns, min_counts=_min_counts)
    if auto_used:
        print(f"Auto-detected peak at bin {start_idx} (t = {data.peak_time_ns:.2f} ns)")

    time_decay, intensity_decay = extract_decay_region(data, start_idx, end_idx)
    print(
        f"Fitting decay window: {time_decay[0]:.2f} - {time_decay[-1]:.2f} ns "
        f"({len(time_decay)} points)"
    )

    fit_range_str = f"t{time_decay[0]:.1f}-{time_decay[-1]:.1f}ns"
    base_name = output_name if output_name else f"{input_stem_clean}{suffix}"

    result_row = {
        "file": str(tiff_path),
        "num_bins": int(data.num_bins),
        "time_bin_width_ns": float(data.time_bin_width_ns),
        "fit_start_ns": float(time_decay[0]),
        "fit_end_ns": float(time_decay[-1]),
        "mono_tau_ns": np.nan,
        "mono_tau_error_ns": np.nan,
        "mono_chi2_red": np.nan,
        "biexp_tau1_ns": np.nan,
        "biexp_tau1_error_ns": np.nan,
        "biexp_tau2_ns": np.nan,
        "biexp_tau2_error_ns": np.nan,
        "biexp_tau_avg_ns": np.nan,
        "biexp_chi2_red": np.nan,
        "denoise": bool(denoise),
        "hot_pixel_k": float(hot_pixel_k) if denoise else np.nan,
        "mask_method": "otsu" if denoise else "",
        "mask_threshold": float(artifacts.mask_threshold) if artifacts is not None else np.nan,
        "n_pixels_in_mask": int(artifacts.n_pixels_in_mask) if artifacts is not None else np.nan,
        "trim_artifact_detected": bool(trim_report.get("detected", False)),
        "trim_cutoff_idx": int(trim_report.get("cutoff_idx", data.num_bins)),
        "trim_n_removed": int(trim_report.get("n_removed", 0)),
    }

    run_mono = model in {"mono", "both", "bi"}
    run_bi = model in {"bi", "both"}

    mono_result = None
    if run_mono:
        print("\n>>> Fitting MONO-EXPONENTIAL model...")
        mono_result = fit_mono_exponential(time_decay, intensity_decay, counts=intensity_decay)
        print_mono_result(mono_result)

        result_row["mono_tau_ns"] = float(mono_result.tau)
        result_row["mono_tau_error_ns"] = float(mono_result.tau_error)
        result_row["mono_chi2_red"] = float(mono_result.chi2_red)

        output_file = output_dir / f"{base_name}_mono_{fit_range_str}.svg"
        plot_mono_exponential_fit(
            data, mono_result, time_decay, intensity_decay,
            output_path=output_file, show_plot=show_plot,
        )

    if run_bi:
        if mono_result is None:
            print("\n>>> Fitting mono-exponential first for initial estimates...")
            mono_result = fit_mono_exponential(time_decay, intensity_decay, counts=intensity_decay)
            print(f"    Mono-exp tau = {mono_result.tau:.3f} ns (used for bi-exp init)")
        else:
            print(f"\n>>> Using mono fit (tau = {mono_result.tau:.3f} ns) for bi-exp initialization...")

        print("\n>>> Fitting BI-EXPONENTIAL model...")
        bi_result = fit_bi_exponential(time_decay, intensity_decay, mono_result=mono_result, counts=intensity_decay)
        print_bi_result(bi_result)

        result_row["biexp_tau1_ns"] = float(bi_result.tau1)
        result_row["biexp_tau1_error_ns"] = float(bi_result.tau1_error)
        result_row["biexp_tau2_ns"] = float(bi_result.tau2)
        result_row["biexp_tau2_error_ns"] = float(bi_result.tau2_error)
        result_row["biexp_tau_avg_ns"] = float(bi_result.tau_avg)
        result_row["biexp_chi2_red"] = float(bi_result.chi2_red)

        output_file = output_dir / f"{base_name}_biexp_{fit_range_str}.svg"
        plot_bi_exponential_fit(
            data, bi_result, time_decay, intensity_decay,
            output_path=output_file, show_plot=show_plot,
        )

    return result_row


# =============================================================================
# CLI
# =============================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Read a FLIM TIFF stack and fit its mean intensity decay "
                    "with mono- or bi-exponential models.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument("input", type=str,
                        help="Path to a TIFF stack, or to a folder when --batch is used")
    parser.add_argument("--batch", action="store_true",
                        help="Treat 'input' as a folder and process every .tif/.tiff inside it")
    parser.add_argument("--model", choices=["mono", "bi", "both"], default="both",
                        help="Decay model: mono, bi, or both (default: both)")
    parser.add_argument("--laser-rate", type=float, default=DEFAULT_LASER_REP_RATE_MHZ,
                        help=f"Laser repetition rate in MHz (default: {DEFAULT_LASER_REP_RATE_MHZ})")
    parser.add_argument("--reduce", choices=["mean", "sum"], default="sum",
                        help="Per-frame pixel aggregation (default: sum)")
    parser.add_argument("--roi", type=str, default=None,
                        help="Optional ROI as 'y0,y1,x0,x1' (slice end-exclusive)")
    parser.add_argument("--fit-start", type=float, default=None,
                        help="Manual start time for fit (ns); default: auto peak")
    parser.add_argument("--fit-end", type=float, default=None,
                        help="Manual end time for fit (ns); default: last non-zero")
    parser.add_argument("--output-dir", type=str, default=None,
                        help=f"Output directory (default: {OUTPUT_DIR})")
    parser.add_argument("--output-name", type=str, default=None,
                        help="Custom output filename prefix (single-file mode only)")
    parser.add_argument("--no-show", action="store_true",
                        help="Do not open the plot window (just save the SVG)")

    # Optional denoising of the (N_bins, H, W) stack before aggregation
    parser.add_argument("--denoise", action="store_true",
                        help="Enable hot-pixel removal + intensity-mask cleanup "
                             "before aggregating the decay (default: off)")
    parser.add_argument("--hot-pixel-k", type=float, default=5.0,
                        help="A pixel is hot when value > k * local 3x3 median (default: 5)")


    # Trailing-artifact trimming
    parser.add_argument("--no-trim", action="store_true",
                        help="Disable automatic detection and removal of trailing "
                             "artifact frames (re-injected laser pulse artefact). "
                             "Trimming is ON by default.")
    parser.add_argument("--artifact-spike-factor", type=float, default=2.5,
                        help="A frame is flagged as artifact when its mean exceeds "
                             "this factor times the local background (default: 2.5)")

    args = parser.parse_args()

    output_dir = Path(args.output_dir) if args.output_dir else OUTPUT_DIR
    roi = _parse_roi(args.roi)
    show_plot = not args.no_show

    common_kwargs = dict(
        denoise=args.denoise,
        hot_pixel_k=args.hot_pixel_k,
        trim_artifacts=not args.no_trim,
        artifact_spike_factor=args.artifact_spike_factor,
    )
    summary_name = "lifetime_summary_denoised.csv" if args.denoise else "lifetime_summary.csv"

    input_path = Path(args.input)

    if args.batch:
        if not input_path.is_dir():
            parser.error(f"--batch requires a folder, got: {input_path}")
        tiff_files = sorted(
            list(input_path.glob("*.tif")) + list(input_path.glob("*.tiff"))
        )
        if not tiff_files:
            parser.error(f"No .tif/.tiff files found in {input_path}")
        print(f"Batch mode: {len(tiff_files)} TIFF files in {input_path}")
        summary_rows = []
        for tp in tiff_files:
            try:
                result_row = process_tiff(
                    tp, output_dir, args.model, args.laser_rate, args.reduce,
                    roi, args.fit_start, args.fit_end, show_plot=show_plot,
                    output_name=None, **common_kwargs,
                )
                summary_rows.append(result_row)
            except Exception as exc:  # keep batch going
                print(f"  !! Failed on {tp.name}: {exc}")

        if summary_rows:
            summary_path = output_dir / summary_name
            pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
            print(f"\nLifetime summary saved to: {summary_path}")
    else:
        if not input_path.is_file():
            parser.error(f"Input file does not exist: {input_path}")
        result_row = process_tiff(
            input_path, output_dir, args.model, args.laser_rate, args.reduce,
            roi, args.fit_start, args.fit_end, show_plot=show_plot,
            output_name=args.output_name, **common_kwargs,
        )

        summary_path = output_dir / summary_name
        pd.DataFrame([result_row]).to_csv(summary_path, index=False)
        print(f"Lifetime summary saved to: {summary_path}")

    print("\nDone!")


if __name__ == "__main__":
    main()
