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

    cmd = [sys.executable, "-m", "marimo", "run", str(app_file)]
    print("Launching:", " ".join(cmd))
    return subprocess.call(cmd, cwd=str(this_dir))


if __name__ == "__main__":
    raise SystemExit(main())
