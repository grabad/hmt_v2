"""
Runs steps 01-05 in order, forwarding any flags (--quick, --cells, --overwrite).
Stops at the first step that fails.

    python drivers/run_all.py --quick       # smoke test on DATA_DIR
    python drivers/run_all.py               # full run
"""
import subprocess
import sys
from pathlib import Path

STEPS = ["01_preprocess.py", "02_calibrate.py", "03_cluster.py", "04_spatial.py", "05_compare.py"]

here = Path(__file__).resolve().parent
for step in STEPS:
    print(f"\n===== {step} =====", flush=True)
    result = subprocess.run([sys.executable, str(here / step), *sys.argv[1:]])
    if result.returncode != 0:
        sys.exit(f"{step} failed (exit {result.returncode})")
