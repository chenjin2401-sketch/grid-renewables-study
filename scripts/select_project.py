"""Pick the case-study solar project from the Phase 2 hosting-capacity map.

Rule (stated up front so the choice is not cherry-picked):
  - sunniest area (highest utility-PV capacity factor);
  - point of interconnection limited by a branch thermal rating (the common case in
    transmission interconnection studies), with firm hosting capacity of 50-300 MW
    (window widened to 500, then 800 MW if no bus qualifies);
  - among those, the bus with the median firm hosting capacity;
  - project size = 2 x firm hosting capacity, rounded to 10 MW (a request well above firm).
Writes results/project.json.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gridstudy.rts import load_rts, load_timeseries   # noqa: E402

RES = ROOT / "results"


def main():
    rts = load_rts()
    ts = load_timeseries(rts)
    cf = {a: float(s.mean()) for a, s in ts.pv_shape_area.items()}
    area = max(cf, key=cf.get)
    hc = pd.read_csv(RES / "phase2_hosting_capacity.csv")
    cand = pd.DataFrame()
    for hi in (300, 500, 800):                     # widen the size window only if needed
        cand = hc[(hc.area == area) & hc.ac_limit.str.startswith("thermal") & hc.ac_firm_mw.between(50, hi)]
        if not cand.empty:
            break
    if cand.empty:
        cand = hc[(hc.area == area) & (hc.ac_firm_mw > 0)]
    cand = cand.sort_values("ac_firm_mw").reset_index(drop=True)
    pick = cand.iloc[len(cand) // 2]
    proj = {"bus": int(pick.bus), "area": int(area), "area_pv_capacity_factor": round(cf[area], 4),
            "firm_mw": float(pick.ac_firm_mw), "dc_firm_mw": float(pick.dc_firm_mw),
            "mw": float(10 * np.round(2 * pick.ac_firm_mw / 10)),
            "limit": pick.ac_limit, "candidates": cand[["bus", "ac_firm_mw", "ac_limit"]].to_dict("records")}
    (RES / "project.json").write_text(json.dumps(proj, indent=2))
    print(json.dumps({k: v for k, v in proj.items() if k != "candidates"}, indent=2))


if __name__ == "__main__":
    main()
