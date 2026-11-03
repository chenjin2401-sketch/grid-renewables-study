"""Phase 5: sensitivities and tipping points (full-year DC-OPF reruns, run in parallel).

  A. load growth 0-20 %         -> cost, LMP, congestion, project firm / flexible size, N-1 (DC)
  B. existing utility PV x1-x2  -> renewable curtailment, LMP, project firm size
  C. battery cost               -> break-even $/kWh against the cheapest line upgrade (from Phase 4)
"""
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

RES, FIG, CACHE = ROOT / "results", ROOT / "figures", ROOT / "cache"


def scenario(args):
    name, load_scale, pv_scale = args
    sys.path.insert(0, str(ROOT / "src"))
    import warnings
    warnings.filterwarnings("ignore")
    from gridstudy.contingency import critical_hours, dc_screen, project_dispatch
    from gridstudy.dcopf import DCOPF
    from gridstudy.hosting import firm_from_headroom, opf_headroom
    from gridstudy.mitigation import curtailment, size_for_curtailment
    from gridstudy.rts import load_rts, load_timeseries
    rts = load_rts()
    ts = load_timeseries(rts, load_scale=load_scale)
    pv = (rts.gen.type == "PV").values
    ts.avail[:, pv] *= pv_scale
    opf = DCOPF(rts)
    nh = len(ts.index)
    pg = np.zeros((nh, len(rts.gen)))
    flow = np.zeros((nh, len(rts.branch)))
    lmp = np.zeros((nh, len(rts.bus)))
    cost = np.zeros(nh)
    shed = np.zeros(nh)
    for h in range(nh):
        r = opf.solve(ts.pd_bus[h] - ts.rtpv_bus[h], avail=ts.avail[h])
        pg[h], flow[h], lmp[h], cost[h], shed[h] = r.pg, r.flow, r.lmp, r.cost, r.shed
    proj = json.loads((RES / "project.json").read_text())
    shape = ts.pv_shape_area[proj["bus"] // 100]
    sunny = np.where(shape >= 0.05)[0]
    head = np.full(nh, 3000.0)
    head[sunny] = opf_headroom(opf, ts, proj["bus"], sunny)
    firm = firm_from_headroom(head[sunny], shape[sunny])
    var = rts.gen.type.isin(["PV", "WIND", "HYDRO"]).values
    rate = rts.branch.rate_a.values
    y = {"pg": pg, "flow": flow}
    hours = critical_hours(ts, y, dc_limit=opf.rating)
    n1_base, _ = dc_screen(rts, flow[hours])
    _, flow_p, _ = project_dispatch(rts, ts, y, hours, {"bus": proj["bus"], "mw": proj["mw"]}, opf=opf)
    n1_proj, _ = dc_screen(rts, flow_p)
    return {
        "scenario": name, "load_scale": load_scale, "pv_scale": pv_scale,
        "annual_cost_musd": cost.sum() / 1e6, "avg_lmp": lmp.mean(),
        "lmp_spread_hours": int((lmp.max(1) - lmp.min(1) > 1).sum()),
        "congested_hours": int(np.isclose(np.abs(flow), opf.rating, atol=1e-3).any(1).sum()),
        "renewable_curtailment_pct": 100 * (1 - pg[:, var].sum() / np.nan_to_num(ts.avail[:, var]).sum()),
        "unserved_mwh": shed.sum(),
        "project_firm_mw_dc_network": firm,
        "project_size_at_3pct_mw": size_for_curtailment(shape, head, 3),
        "project_curtailment_pct": 100 * curtailment(proj["mw"], shape, head)[0],
        "n1_dc_overloads_baseline": len(n1_base), "n1_dc_overloads_with_project": len(n1_proj),
    }


def main():
    cases = [("load +0%", 1.0, 1.0), ("load +5%", 1.05, 1.0), ("load +10%", 1.10, 1.0),
             ("load +15%", 1.15, 1.0), ("load +20%", 1.20, 1.0),
             ("PV x1.5", 1.0, 1.5), ("PV x2", 1.0, 2.0), ("PV x2.5", 1.0, 2.5)]
    with ProcessPoolExecutor(max_workers=2) as ex:
        rows = list(ex.map(scenario, cases))
    df = pd.DataFrame(rows)
    df.round(3).to_csv(RES / "phase5_sensitivities.csv", index=False)
    print(df.round(2).to_string(index=False))

    # C. battery break-even vs. cheapest line upgrade
    bat = pd.read_csv(RES / "phase4_battery_sweep.csv")
    up = pd.read_csv(RES / "phase4_line_upgrades.csv").dropna(subset=["usd_per_avoided_mwh"])
    out = {"scenarios": df.round(3).to_dict("records")}
    if len(up) and (bat.avoided_mwh > 0).any():
        best = up.loc[up.usd_per_avoided_mwh.idxmin()]
        sys.path.insert(0, str(ROOT / "src"))
        from gridstudy.mitigation import annualized
        fom, rate, years = 0.04, 0.07, 20
        be = []
        for _, r in bat[(bat["mode"] == "capture") & (bat.avoided_mwh > 0)].iterrows():
            # $/kWh at which the battery's $/MWh-avoided equals the best line upgrade
            per_kwh = best.usd_per_avoided_mwh * r.avoided_mwh / annualized(r.battery_mwh * 1000, rate, years, fom)
            be.append({"battery_mw": int(r.battery_mw), "battery_mwh": int(r.battery_mwh),
                       "avoided_mwh": round(r.avoided_mwh), "breakeven_usd_per_kwh": round(per_kwh, 1)})
        out["breakeven_vs"] = best.to_dict()
        out["battery_breakeven"] = be
        print(pd.DataFrame(be).to_string(index=False))
    (RES / "phase5_summary.json").write_text(json.dumps(out, indent=2, default=float))
    figures(df)


def figures(df):
    load = df[df.pv_scale == 1.0]
    pvs = df[df.load_scale == 1.0]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.3))
    g = 100 * (load.load_scale - 1)
    axes[0].plot(g, load.project_firm_mw_dc_network, "o-", label="firm (no curtailment)")
    axes[0].plot(g, load.project_size_at_3pct_mw, "s-", label="flexible (<=3 % curtailment)")
    axes[0].set_xlabel("load growth (%)")
    axes[0].set_ylabel("project size the grid can take (MW)")
    axes[0].set_title("Project bus: DC network headroom vs. load growth")
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    axes[1].plot(g, load.n1_dc_overloads_baseline, "o-", label="baseline")
    axes[1].plot(g, load.n1_dc_overloads_with_project, "s-", label="with project")
    axes[1].set_xlabel("load growth (%)")
    axes[1].set_ylabel("N-1 overloads on critical hours (DC screen)")
    axes[1].set_title("Contingency overloads vs. load growth")
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    axes[2].plot(pvs.pv_scale, pvs.renewable_curtailment_pct, "o-", label="system renewable curtailment (%)")
    ax2 = axes[2].twinx()
    ax2.plot(pvs.pv_scale, pvs.avg_lmp, "s--", color="tab:orange", label="average LMP ($/MWh)")
    axes[2].set_xlabel("existing utility PV fleet (x today)")
    axes[2].set_ylabel("curtailment (%)")
    ax2.set_ylabel("average LMP ($/MWh)")
    axes[2].set_title("More solar: curtailment and prices")
    axes[2].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG / "phase5_sensitivities.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
