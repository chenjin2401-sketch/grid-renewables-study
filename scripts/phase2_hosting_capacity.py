"""Phase 2: solar hosting capacity at every bus of RTS-GMLC (N-0, firm, unity power factor).

1. DC screen (OPF-based): hourly network headroom on about 700 critical sunny hours per bus
   (sunniest, lowest net load, most congested, random), firm capacity = min headroom / shape.
2. AC verification on each bus's 24 tightest/sunniest hours with the redispatched fleet,
   bisecting down if AC shows a new thermal or voltage violation.
Buses run in parallel (4 processes). Needs cache/phase1.npz.
Outputs results/phase2_hosting_capacity.csv, results/phase2_summary.json, figures/phase2_*.png.
"""
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

CACHE, RES, FIG = ROOT / "cache", ROOT / "results", ROOT / "figures"
_W = {}


def _init():
    import warnings
    warnings.filterwarnings("ignore")
    from gridstudy.dcopf import DCOPF
    from gridstudy.rts import build_ac, load_rts, load_timeseries
    rts = load_rts()
    with np.load(CACHE / "phase1.npz") as z:                # only the arrays the study needs
        y = {k: z[k] for k in ("flow", "conv", "s_ac", "vm")}
    _W.update(rts=rts, ts=load_timeseries(rts), y=y, opf=DCOPF(rts), ac=build_ac(rts))


def study_bus(b, q_of_v=None):
    """Hosting capacity at bus b; q_of_v = Volt-VAR curve for the new plant (None = unity pf)."""
    import pandapower as pp
    from gridstudy.hosting import ac_critical_hours, ac_hosting_capacity, candidate_hours, firm_from_headroom, \
        opf_headroom
    rts, ts, y, opf, ac = (_W[k] for k in ("rts", "ts", "y", "opf", "ac"))
    area = b // 100
    shape = ts.pv_shape_area[area]
    hours = candidate_hours(ts, y, area, dc_limit=opf.rating)
    head = opf_headroom(opf, ts, b, hours)
    dc_cap = firm_from_headroom(head, shape[hours])
    worst = hours[np.argmin(np.where(shape[hours] >= 0.05, head / np.maximum(shape[hours], 1e-9), np.inf))]
    sg = pp.create_sgen(ac.net, b, p_mw=0.0, q_mvar=0.0, name="new solar")
    ac_hours = ac_critical_hours(hours, head, shape[hours])
    # add light-load sunny hours and hours with the highest baseline voltage near the POI
    # (Phase 4 showed these set the AC limit at voltage-limited buses)
    from gridstudy.hosting import study_area
    sunny = hours[shape[hours] >= 0.3]
    net = (ts.pd_bus.sum(1) - ts.rtpv_bus.sum(1))[sunny]
    vloc = np.nan_to_num(y["vm"][sunny][:, study_area(rts, b)]).max(1)
    ac_hours = np.unique(np.r_[ac_hours, sunny[np.argsort(net)[:12]], sunny[np.argsort(-vloc)[:12]]])
    ac_cap, why, where = ac_hosting_capacity(rts, opf, ac, ts, y, b, dc_cap, ac_hours, sg, q_of_v=q_of_v)
    ac.net.sgen.drop(sg, inplace=True)
    return {"bus": int(b), "area": area, "kv": float(rts.bus.base_kv[rts.bus.id == b].iloc[0]),
            "hours_screened": int(len(hours)), "dc_firm_mw": round(dc_cap, 1),
            "dc_limit_hour": str(ts.index[worst]), "ac_firm_mw": round(ac_cap, 1), "ac_limit": why,
            "ac_limit_hour": str(ts.index[where]) if where is not None else ""}


def main():
    sys.path.insert(0, str(ROOT / "src"))
    from gridstudy.rts import load_rts
    rts = load_rts()
    t0 = time.perf_counter()
    partial = CACHE / "phase2_partial.csv"                  # resume after an interruption
    done = pd.read_csv(partial) if partial.exists() else pd.DataFrame()
    rows = done.to_dict("records")
    todo = [int(b) for b in rts.bus.id if done.empty or b not in set(done.bus)]
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 2   # 2 processes fit in ~2 GB of free RAM
    with ProcessPoolExecutor(max_workers=workers, initializer=_init) as ex:
        for r in ex.map(study_bus, todo):
            rows.append(r)
            pd.DataFrame(rows).to_csv(partial, index=False)
            print(f"bus {r['bus']}: DC {r['dc_firm_mw']:7.1f} MW, AC {r['ac_firm_mw']:7.1f} MW  {r['ac_limit']}",
                  flush=True)
    df = pd.DataFrame(rows).sort_values("bus").reset_index(drop=True)
    df.to_csv(RES / "phase2_hosting_capacity.csv", index=False)
    hc = df.ac_firm_mw.values
    pos = hc[hc > 0]
    best = df.loc[df.ac_firm_mw.idxmax()]
    worst = df.loc[df.ac_firm_mw.idxmin()]
    summary = {
        "buses_studied": int(len(df)), "hours_screened_per_bus": int(df.hours_screened.median()),
        "runtime_s": round(time.perf_counter() - t0, 1),
        "ac_firm_hc_mw": {"min": float(hc.min()), "p25": float(np.percentile(hc, 25)), "median": float(np.median(hc)),
                          "p75": float(np.percentile(hc, 75)), "max": float(hc.max())},
        "total_firm_hc_gw_if_each_bus_alone": round(float(hc.sum()) / 1000, 2),
        "best_bus": int(best.bus), "best_bus_hc_mw": float(best.ac_firm_mw),
        "worst_bus": int(worst.bus), "worst_bus_hc_mw": float(worst.ac_firm_mw),
        "best_to_worst_ratio": round(float(hc.max() / pos.min()), 1) if len(pos) else None,
        "best_to_median_ratio": round(float(hc.max() / np.median(hc)), 2),
        "buses_zero_hc": int((hc == 0).sum()),
        "share_limited_by_dc_thermal_pct": round(100 * df.ac_limit.str.startswith("thermal (DC").mean(), 1),
        "share_ac_thermal_pct": round(100 * df.ac_limit.str.match(r"thermal [A-Z]").mean(), 1),
        "share_ac_voltage_pct": round(100 * df.ac_limit.str.startswith("voltage").mean(), 1),
        "buses_where_ac_below_dc": int((df.ac_firm_mw < df.dc_firm_mw - 5).sum()),
        "median_ac_reduction_where_lower_pct": round(float(np.median(
            (100 * (1 - df.ac_firm_mw / df.dc_firm_mw))[df.ac_firm_mw < df.dc_firm_mw - 5])), 1)
        if (df.ac_firm_mw < df.dc_firm_mw - 5).any() else 0.0,
        "median_hc_by_kv": {str(int(k)): float(g.ac_firm_mw.median()) for k, g in df.groupby("kv")},
        "median_hc_by_area": {str(a): float(g.ac_firm_mw.median()) for a, g in df.groupby("area")},
    }
    (RES / "phase2_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    plot(rts, df)


def plot(rts, df):
    from gridstudy.rts import DATA
    geo = pd.read_csv(DATA / "SourceData" / "bus.csv").set_index("Bus ID")
    fig, ax = plt.subplots(figsize=(10, 8))
    for _, br in rts.branch.iterrows():
        a, b = geo.loc[br.f], geo.loc[br.t]
        ax.plot([a.lng, b.lng], [a.lat, b.lat], color="0.75", lw=0.8, zorder=1)
    d = df.set_index("bus")
    sc = ax.scatter(geo.loc[d.index, "lng"], geo.loc[d.index, "lat"], c=d.ac_firm_mw, cmap="viridis",
                    s=30 + d.ac_firm_mw / 10, edgecolor="k", lw=0.5, zorder=2)
    for b, r in d.iterrows():
        if r.ac_firm_mw >= d.ac_firm_mw.quantile(0.9) or r.ac_firm_mw <= d.ac_firm_mw.quantile(0.1):
            ax.annotate(f"{b}: {r.ac_firm_mw:.0f}", (geo.loc[b, "lng"], geo.loc[b, "lat"]), fontsize=7,
                        xytext=(4, 4), textcoords="offset points")
    plt.colorbar(sc, ax=ax, label="Firm solar hosting capacity (MW, AC-verified)")
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.set_title("RTS-GMLC solar hosting capacity by bus (N-0, no curtailment, unity pf)")
    fig.tight_layout()
    fig.savefig(FIG / "phase2_hosting_capacity_map.png", dpi=150)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(13, 4.2))
    s = df.sort_values("ac_firm_mw", ascending=False)
    col = s.ac_limit.map(lambda x: "tab:red" if x.startswith("voltage") else
                         ("tab:orange" if x.startswith("thermal ") and "(DC" not in x else "tab:blue"))
    ax.bar(s.bus.astype(str), s.ac_firm_mw, color=col)
    ax.bar(s.bus.astype(str), (s.dc_firm_mw - s.ac_firm_mw).clip(lower=0), bottom=s.ac_firm_mw, color="0.85")
    ax.set_ylabel("MW")
    ax.set_title("Firm hosting capacity by bus. Blue: DC thermal limit holds in AC; orange: lower AC thermal "
                 "limit; red: voltage limit; grey: DC estimate removed by the AC check", fontsize=9)
    ax.tick_params(axis="x", labelrotation=90, labelsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "phase2_hosting_capacity_bars.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
