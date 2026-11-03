"""Phase 3: automated N-1 contingency analysis on critical hours, with and without the solar project.

For each case (baseline, with project):
  - full AC N-1 (every non-islanding branch outage, every critical hour) as the benchmark
  - DC LODF screening at several thresholds + AC verification of flagged cases
  - speed-up and missed overloads of screening vs. the full AC run

Needs cache/phase1.npz and results/phase2_hosting_capacity.csv.
The project bus/size come from results/project.json (written by phase 4 prep) or the CLI.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pandapower as pp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gridstudy.contingency import ac_base_case, ac_contingencies, critical_hours, dc_screen, project_dispatch  # noqa: E402
from gridstudy.rts import build_ac, load_rts, load_timeseries                                    # noqa: E402

CACHE, RES = ROOT / "cache", ROOT / "results"
THRESHOLDS = [1.0, 0.95, 0.9, 0.85]


def run_case(rts, ts, y, hours, project):
    ac = build_ac(rts)
    sgen = None
    if project is not None:
        sgen = pp.create_sgen(ac.net, project["bus"], p_mw=0.0, name="project solar")
    pg, flow, p_new = project_dispatch(rts, ts, y, hours, project)
    # intact-system AC: hours that do not converge even without an outage are reported, not screened
    base_conv, q_base = ac_base_case(rts, ac, ts, hours, pg, p_new, sgen=sgen)
    excluded = [str(ts.index[h]) for h in hours[~base_conv]]
    hours, pg, flow, p_new, q_base = hours[base_conv], pg[base_conv], flow[base_conv], p_new[base_conv], \
        q_base[base_conv]

    full, n_full, failed, t_full = ac_contingencies(rts, ac, ts, y, hours, pg, p_new, sgen=sgen)
    truth = set(zip(full.hour_pos, full.outage, full.monitored))
    truth_pairs = set(zip(full.hour_pos, full.outage))
    screens = []
    for method in ("DC (LODF)", "hybrid (LODF + base-case Q)"):
        for thr in THRESHOLDS:
            t0 = time.perf_counter()
            flagged, island = dc_screen(rts, flow, threshold=thr, q_base=q_base if method.startswith("hybrid") else None)
            pairs = sorted(set(zip(flagged.hour_pos, flagged.outage)))
            t_dc = time.perf_counter() - t0
            ver, n_ver, _, t_ver = ac_contingencies(rts, ac, ts, y, hours, pg, p_new, pairs=pairs, sgen=sgen)
            found = set(zip(ver.hour_pos, ver.outage, ver.monitored))
            screens.append({"method": method, "threshold": thr, "flagged_cases": len(pairs), "ac_runs": n_ver,
                            "runtime_s": round(t_dc + t_ver, 2), "overloads_found": len(found & truth),
                            "overloads_missed": len(truth - found),
                            "recall_pct": round(100 * len(found & truth) / max(len(truth), 1), 1),
                            "nonconverged_cases_flagged": len(set(failed) & set(pairs)),
                            "speedup_vs_full_ac": round(t_full / max(t_dc + t_ver, 1e-9), 1)})
    lab = rts.branch.uid.values
    full = full.assign(hour=[str(ts.index[hours[j]]) for j in full.hour_pos],
                       outage_uid=lab[full.outage], monitored_uid=lab[full.monitored])
    summary = {
        "hours_screened": int(len(hours)), "hours_excluded_base_not_converged": excluded,
        "contingencies_per_hour": int((~island).sum()),
        "islanding_outages_skipped": int(island.sum()),
        "full_ac_runs": n_full, "full_ac_runtime_s": round(t_full, 1), "non_converged": len(failed),
        "overloads": len(truth), "contingency_hours_with_overload": len(truth_pairs),
        "distinct_critical_outages": int(full.outage.nunique()) if len(full) else 0,
        "distinct_overloaded_branches": int(full.monitored.nunique()) if len(full) else 0,
        "worst_loading_pct_of_lte": round(100 * full.ac_loading.max(), 1) if len(full) else 0.0,
        "overloads_on_C10_reactive": int((full.monitored_uid == "C10").sum()) if len(full) else 0,
    }
    return summary, pd.DataFrame(screens), full


def main():
    rts = load_rts()
    ts = load_timeseries(rts)
    y = dict(np.load(CACHE / "phase1.npz"))
    project = json.loads((RES / "project.json").read_text())
    hours = critical_hours(ts, y, dc_limit=rts.branch.rate_a.values / rts.mva_factor)
    # add the sunniest hours (when the project produces most)
    shape = ts.pv_shape_area[project["bus"] // 100]
    hours = np.unique(np.r_[hours, np.argsort(-shape)[:12]])
    out = {}
    for name, proj in (("baseline", None), ("with_project", {"bus": project["bus"], "mw": project["mw"]})):
        s, scr, full = run_case(rts, ts, y, hours, proj)
        out[name] = {"summary": s, "screening": scr.to_dict("records")}
        full.to_csv(RES / f"phase3_overloads_{name}.csv", index=False)
        print(name, json.dumps(s, indent=2))
        print(scr.to_string(index=False))
    out["project"] = {k: v for k, v in project.items() if k != "candidates"}
    b, p = out["baseline"]["summary"], out["with_project"]["summary"]
    ob = pd.read_csv(RES / "phase3_overloads_baseline.csv")
    op = pd.read_csv(RES / "phase3_overloads_with_project.csv")
    kb = set(zip(ob.hour, ob.outage_uid, ob.monitored_uid))
    kp = set(zip(op.hour, op.outage_uid, op.monitored_uid))
    new = pd.DataFrame(sorted(kp - kb), columns=["hour", "outage", "monitored"])
    out["delta"] = {"overloads_baseline": b["overloads"], "overloads_with_project": p["overloads"],
                    "new_overload_cases": int(len(new)), "relieved_cases": int(len(kb - kp)),
                    "top_new_outage_monitored": {f"{o} -> {m}": int(n) for (o, m), n in
                                                 new.groupby(["outage", "monitored"]).size()
                                                 .sort_values(ascending=False).head(6).items()},
                    "worst_loading_change_pct_points": round(p["worst_loading_pct_of_lte"] - b["worst_loading_pct_of_lte"], 1)}
    (RES / "phase3_summary.json").write_text(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
