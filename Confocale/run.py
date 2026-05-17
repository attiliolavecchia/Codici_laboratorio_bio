from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def main() -> int:
    this_dir = Path(__file__).resolve().parent
    app_file = this_dir / "andi_demo.py"

    if not app_file.exists():
        print(f"Cannot find app file: {app_file}")
        return 1

    # --wasm  → export as self-contained HTML for GitHub Pages
    if "--wasm" in sys.argv:
        out_dir = this_dir / "dist"
        out_dir.mkdir(exist_ok=True)
        cmd = [
            sys.executable, "-m", "marimo", "export", "html-wasm",
            str(app_file),
            "-o", str(out_dir),
            "--mode", "run",
            "--no-show-code",
            "-f",  # overwrite if exists
        ]
        print("Exporting WASM → ", out_dir)
        print("Run:", " ".join(cmd))
        rc = subprocess.call(cmd, cwd=str(this_dir))
        if rc == 0:
            html = out_dir / "index.html"
            print(f"\n✓ Done!  Open  {html}\n  or deploy the '{out_dir.name}/' folder to GitHub Pages.")
        return rc

    mode = "run"
    if "--edit" in sys.argv:
        mode = "edit"

    cmd = [sys.executable, "-m", "marimo", mode, str(app_file)]
    print("Launching:", " ".join(cmd))
    return subprocess.call(cmd, cwd=str(this_dir))


if __name__ == "__main__":
    raise SystemExit(main())
