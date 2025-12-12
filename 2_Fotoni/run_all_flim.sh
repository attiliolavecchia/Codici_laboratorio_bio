#!/usr/bin/env bash
set -euo pipefail

# Batch runner for FLIM fits.
# - Finds all Plot_Values_*.csv under a data root
# - Runs mono and bi fits with auto window (peak -> last non-zero)
# - Writes outputs to an img root mirroring the data subfolders
# - Suppresses plot display (--no-show)
#
# Usage:
#   ./run_all_flim.sh [--data-root PATH] [--output-root PATH] [--laser-rate MHz] [--python PATH]
#
# Defaults assume this file lives next to flim_exponential_fit.py with data/ and img/ siblings.

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHON_BIN="${PYTHON:-python}"
DATA_ROOT="$SCRIPT_DIR/data"
OUTPUT_ROOT="$SCRIPT_DIR/img"
LASER_RATE="80"

usage() {
  cat <<'EOF'
Usage: run_all_flim.sh [--data-root PATH] [--output-root PATH] [--laser-rate MHz] [--python PATH]

Options:
  --data-root PATH    Root folder containing Plot_Values_*.csv files (default: ./data)
  --output-root PATH  Root folder for outputs, mirrors data structure (default: ./img)
  --laser-rate MHz    Laser repetition rate in MHz (default: 80)
  --python PATH       Python executable to use (default: $PYTHON or python)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-root)
      DATA_ROOT="$2"; shift 2;;
    --output-root)
      OUTPUT_ROOT="$2"; shift 2;;
    --laser-rate)
      LASER_RATE="$2"; shift 2;;
    --python)
      PYTHON_BIN="$2"; shift 2;;
    -h|--help)
      usage; exit 0;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 1;;
  esac
done

if [[ ! -d "$DATA_ROOT" ]]; then
  echo "Data root not found: $DATA_ROOT" >&2
  exit 1
fi

# Normalize paths (remove trailing slashes)
DATA_ROOT="${DATA_ROOT%/}"
OUTPUT_ROOT="${OUTPUT_ROOT%/}"

common_args=(--no-show --laser-rate "$LASER_RATE" --x-column '[cm]' --y-column 'Mean')

# Use find with -L to follow symlinks and suppress errors for broken links
find "$DATA_ROOT" -type f -name 'Plot_Values_*.csv' 2>/dev/null | while IFS= read -r csv_path; do
  rel_path="${csv_path#$DATA_ROOT}"
  rel_path="${rel_path#/}"  # Remove leading slash if present
  rel_dir="$(dirname "$rel_path")"
  
  # If rel_dir is ".", just use OUTPUT_ROOT
  if [[ "$rel_dir" == "." ]]; then
    out_dir="$OUTPUT_ROOT"
  else
    out_dir="$OUTPUT_ROOT/$rel_dir"
  fi

  echo "Processing: $csv_path"
  echo "  Output -> $out_dir"

  "$PYTHON_BIN" "$SCRIPT_DIR/flim_exponential_fit.py" "${common_args[@]}" --output-dir "$out_dir" --model mono "$csv_path"
  "$PYTHON_BIN" "$SCRIPT_DIR/flim_exponential_fit.py" "${common_args[@]}" --output-dir "$out_dir" --model bi "$csv_path"

done
