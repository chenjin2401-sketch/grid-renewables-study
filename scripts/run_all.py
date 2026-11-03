"""Run every phase in order (a few hours on a laptop with 2 worker processes).

    python scripts/download_rts_gmlc.py      # once
    python scripts/run_all.py
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STEPS = ["phase1_rts_baseline.py", "phase2_hosting_capacity.py", "select_project.py",
         "phase3_n1_screening.py", "phase4_mitigations.py", "phase4b_volt_var_map.py",
         "phase5_sensitivities.py"]

if __name__ == "__main__":
    if (ROOT.parent / "data" / "private" / "load_profile_5bus.csv").exists():
        STEPS.insert(0, "phase0_reproduce.py")
    for step in STEPS:
        print(f"\n=== {step} ===", flush=True)
        subprocess.run([sys.executable, str(ROOT / step)], check=True)
