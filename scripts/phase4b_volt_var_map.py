"""Phase 4b: does inverter Volt-VAR raise hosting capacity where it matters?

Reruns the Phase 2 study (same hours, same criteria) with the new plant on the IEEE 1547-2018
Category B Volt-VAR curve instead of unity power factor, at every bus whose Phase 2 limit was
its own voltage (the plant raising the voltage of the bus it connects to).
Outputs results/phase4b_volt_var_hc.csv and results/phase4b_summary.json.
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

RES = ROOT / "results"


def run(b):
    import phase2_hosting_capacity as p2
    from gridstudy.hosting import volt_var_ieee1547
    if not p2._W:
        p2._init()
    return p2.study_bus(b, q_of_v=volt_var_ieee1547)


def main():
    hc = pd.read_csv(RES / "phase2_hosting_capacity.csv")
    own = hc[hc.ac_limit == "voltage " + hc.bus.astype(str)]
    with ProcessPoolExecutor(max_workers=2) as ex:
        rows = list(ex.map(run, own.bus.astype(int).tolist()))
    vv = pd.DataFrame(rows)[["bus", "ac_firm_mw", "ac_limit"]].rename(
        columns={"ac_firm_mw": "volt_var_mw", "ac_limit": "volt_var_limit"})
    df = own[["bus", "kv", "dc_firm_mw", "ac_firm_mw", "ac_limit"]].rename(
        columns={"ac_firm_mw": "unity_pf_mw", "ac_limit": "unity_pf_limit"}).merge(vv, on="bus")
    df["gain_mw"] = df.volt_var_mw - df.unity_pf_mw
    df["gain_pct"] = 100 * df.gain_mw / df.unity_pf_mw
    df.round(1).to_csv(RES / "phase4b_volt_var_hc.csv", index=False)
    s = {"buses": int(len(df)), "buses_improved_gt_10pct": int((df.gain_pct > 10).sum()),
         "median_gain_pct": round(float(df.gain_pct.median()), 1),
         "median_unity_mw": float(df.unity_pf_mw.median()), "median_volt_var_mw": float(df.volt_var_mw.median()),
         "total_gain_mw": round(float(df.gain_mw.sum()), 1),
         "new_limits": df.volt_var_limit.str.extract(r"^(\w+)")[0].value_counts().to_dict()}
    (RES / "phase4b_summary.json").write_text(json.dumps(s, indent=2))
    print(df.round(1).to_string(index=False))
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main()
