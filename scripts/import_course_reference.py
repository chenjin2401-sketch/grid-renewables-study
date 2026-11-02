"""Copy the Phase 0 validation inputs from the MATLAB coursework folder into data/private/.

Usage: python scripts/import_course_reference.py "<MATLAB project folder>"

Reads the 24 x 5 load profile from project_config.m and copies the MATPOWER
hourly results (results/task*_*.csv) as the reference for the Python port.
"""
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "data" / "private"


def read_load_profile(config_m: Path) -> list[list[float]]:
    text = config_m.read_text(encoding="utf-8")
    m = re.search(r"cfg\.Pload\s*=\s*\[(.*?)\];", text, re.S)
    if not m:
        raise ValueError(f"cfg.Pload not found in {config_m}")
    rows = []
    for line in m.group(1).strip().splitlines():
        vals = [float(v) for v in line.split()]
        if vals:
            rows.append(vals)
    if len(rows) != 24 or any(len(r) != 5 for r in rows):
        raise ValueError("expected a 24 x 5 load profile")
    return rows


def main(src: Path) -> None:
    DEST.mkdir(parents=True, exist_ok=True)
    rows = read_load_profile(src / "project_config.m")
    with open(DEST / "load_profile_5bus.csv", "w", encoding="utf-8") as f:
        f.write("hour,bus1,bus2,bus3,bus4,bus5\n")
        for h, r in enumerate(rows):
            f.write(f"{h}," + ",".join(f"{v:g}" for v in r) + "\n")
    ref = DEST / "matpower_reference"
    ref.mkdir(exist_ok=True)
    n = 0
    for p in sorted((src / "results").glob("task*_*.csv")):
        shutil.copy2(p, ref / p.name)
        n += 1
    print(f"wrote {DEST / 'load_profile_5bus.csv'} and {n} reference files to {ref}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(Path(sys.argv[1]))
