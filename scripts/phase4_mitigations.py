"""Phase 4: what to do when a solar project is larger than the firm hosting capacity.

Project: results/project.json (bus, size MW; chosen by select_project.py from the Phase 2 map).
Hourly network headroom at the project bus for the whole year:
  - DC (OPF-based) for all 8,784 hours,
  - AC-verified (same criteria as Phase 2) for every sunny hour, by bisection below the DC value.
Cases compared:
  0. firm interconnection (size = firm hosting capacity)
  1. flexible interconnection: full size, curtailed when the grid is full
  2. + co-located 4-hour battery (size sweep), charged only from the plant
  3. + reconductor / rebuild of the branch that limits the most curtailed energy
  4. inverter Volt-VAR (IEEE 1547 Cat B) instead of unity power factor
Economics: annualized cost (7 %, 20 years) per MWh of curtailment avoided; battery cost swept.
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

# --- cost assumptions (sources in docs/RESULTS.md) -------------------------------------------
RECONDUCTOR_USD_PER_MILE = {138: 0.37e6, 230: 0.45e6}     # MISO MTEP22 cost guide via Illinois IPA (2024);
REBUILD_USD_PER_MILE = {138: 1.7e6, 230: 2.2e6}           # 230 kV interpolated between the 138 and 345 kV values
RECONDUCTOR_UPRATE = 1.35                                  # reconductoring is enough up to 135 % loading
REBUILD_UPRATE = 2.0
BATTERY_USD_PER_KWH = [150, 250, 350, 450]                 # installed 4-hour system (sensitivity)
BATTERY_FOM = 0.04                                         # NREL ATB: fixed O&M = 4 % of capex per year
BATTERY_RTE = 0.85                                         # NREL ATB round-trip efficiency
RATE, YEARS = 0.07, 20
WORKERS = 2
SWEEP_FACTOR = 1.5          # hourly AC headroom is computed exactly for plants up to 1.5 x the project size
_W = {}


def _init(with_ac):
    import warnings
    warnings.filterwarnings("ignore")
    from gridstudy.dcopf import DCOPF
    from gridstudy.rts import build_ac, load_rts, load_timeseries
    rts = load_rts()
    _W.update(rts=rts, ts=load_timeseries(rts), opf=DCOPF(rts))
    if with_ac:
        import pandapower as pp
        with np.load(CACHE / "phase1.npz") as z:
            _W["y"] = {k: z[k] for k in ("conv", "s_ac", "vm")}
        ac = build_ac(rts)
        _W["ac"] = ac
        _W["sg"] = {}


def _dc_chunk(args):
    hours, bus, rating = args
    rts, ts, opf = _W["rts"], _W["ts"], _W["opf"]
    pos = {b: i for i, b in enumerate(rts.bus.id)}
    out = [opf.max_injection(ts.pd_bus[h] - ts.rtpv_bus[h], ts.avail[h], pos[bus], rating=rating,
                             return_binding=True) for h in hours]
    return np.array([o[0] for o in out]), np.array([o[1] for o in out])


def _ac_chunk(args):
    hours, dc_head, bus, size, need_mw, volt_var, rating = args
    import pandapower as pp
    from gridstudy.hosting import INV_OVERSIZE, ACChecker, volt_var_ieee1547
    ac = _W["ac"]
    if bus not in _W["sg"]:
        _W["sg"][bus] = pp.create_sgen(ac.net, bus, p_mw=0.0, q_mvar=0.0, name="project solar")
    chk = ACChecker(_W["rts"], _W["opf"], ac, _W["ts"], _W["y"], bus, _W["sg"][bus],
                    q_of_v=volt_var_ieee1547 if volt_var else None, s_rated=INV_OVERSIZE * size, rate=rating)
    shape = _W["ts"].pv_shape_area[bus // 100]
    # headroom is exact for every plant up to need_mw (the largest size evaluated)
    res = [chk.hourly_headroom(h, d, need=need_mw * shape[h], try_dc_first=not volt_var)
           for h, d in zip(hours, dc_head)]
    return np.array([r[0] for r in res]), [r[1] for r in res]


def dc_headroom(ex, hours, bus, rating=None):
    parts = list(ex.map(_dc_chunk, [(c, bus, rating) for c in np.array_split(hours, 4 * WORKERS)]))
    return np.concatenate([p[0] for p in parts]), np.concatenate([p[1] for p in parts])


def ac_headroom(ex, hours, dc_head, bus, size, need_mw=None, volt_var=False, rating=None):
    if len(hours) == 0:
        return np.array([]), []
    need_mw = size if need_mw is None else need_mw
    n = max(1, min(8 * WORKERS, len(hours)))
    chunks = [(h, d, bus, size, need_mw, volt_var, rating) for h, d in zip(np.array_split(hours, n),
                                                                           np.array_split(dc_head, n))]
    parts = list(ex.map(_ac_chunk, chunks))
    return np.concatenate([p[0] for p in parts]), sum((p[1] for p in parts), [])


def main():
    from gridstudy.mitigation import annualized, battery_dispatch, curtailment, size_for_curtailment
    from gridstudy.rts import load_rts, load_timeseries
    t0 = time.perf_counter()
    rts = load_rts()
    ts = load_timeseries(rts)
    with np.load(CACHE / "phase1.npz") as z:
        lmp_all = z["lmp"]
    proj = json.loads((RES / "project.json").read_text())
    bus, size = proj["bus"], proj["mw"]
    shape = ts.pv_shape_area[bus // 100]
    pos = {b: i for i, b in enumerate(rts.bus.id)}
    lmp = lmp_all[:, pos[bus]]
    nh = len(ts.index)
    sunny = np.where(shape >= 0.01)[0]
    out = {"project": {k: v for k, v in proj.items() if k != "candidates"}}

    cache_file = CACHE / "phase4_headroom.npz"
    reuse = "--reuse" in sys.argv and cache_file.exists()
    sweep_max = SWEEP_FACTOR * size                 # headroom is exact for plants up to this size
    if reuse:
        with np.load(cache_file) as z:
            reuse = "reasons" in z and int(z["bus"]) == bus and float(z["size"]) == size and \
                float(z.get("need_mw", 0)) == sweep_max
    if reuse:                                   # hourly headroom of a previous run (same project)
        with np.load(cache_file) as z:
            head_dc, head, binding, why = z["head_dc"], z["head"], z["binding"], list(z["reasons"])
        print("reusing cached hourly headroom", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=WORKERS, initializer=_init, initargs=(False,)) as ex:
            head_dc, binding = dc_headroom(ex, np.arange(nh), bus)
        print(f"DC headroom done {time.perf_counter() - t0:.0f} s", flush=True)
    with ProcessPoolExecutor(max_workers=WORKERS, initializer=_init, initargs=(True,)) as ex:
        if not reuse:
            h_ac, why = ac_headroom(ex, sunny, head_dc[sunny], bus, size, need_mw=sweep_max)
            head = head_dc.copy()
            head[sunny] = h_ac
            print(f"AC headroom done {time.perf_counter() - t0:.0f} s", flush=True)
            np.savez_compressed(cache_file, head_dc=head_dc, head=head, binding=binding,
                                reasons=np.array(why), bus=bus, size=size, need_mw=sweep_max)
        reasons = pd.Series(why, index=sunny)

        firm_dc = float(np.min(head_dc[shape >= 0.05] / shape[shape >= 0.05]))
        firm = float(np.min(head[shape >= 0.05] / shape[shape >= 0.05]))
        out["firm_mw_full_year"] = {"dc": round(firm_dc, 1), "ac": round(firm, 1),
                                    "phase2_ac_on_critical_hours": proj["firm_mw"]}

        # 0/1. firm vs flexible --------------------------------------------------------------
        cur, pot, deliv_flex = curtailment(size, shape, head)
        flex = {f"{p}%": round(size_for_curtailment(shape, head, p, hi=sweep_max), 1) for p in (1, 3, 5, 10)}
        out["size_curve_exact_up_to_mw"] = sweep_max
        firm_energy = (firm * shape).sum()
        out["flexible"] = {"size_mw": size, "curtailment_pct": round(100 * cur, 2), "potential_mwh": round(pot),
                           "delivered_mwh": round(deliv_flex), "curtailed_mwh": round(pot - deliv_flex),
                           "firm_plant_delivered_mwh": round(firm_energy),
                           "extra_energy_vs_firm_pct": round(100 * (deliv_flex / firm_energy - 1), 1),
                           "size_at_curtailment_cap_mw": flex,
                           "size_gain_vs_firm_pct": {k: round(100 * (v / firm - 1), 1) for k, v in flex.items()}}
        short = sunny[size * shape[sunny] > head[sunny] + 1e-6]
        lost = pd.Series(size * shape[short] - head[short], index=short)
        lim = reasons.loc[short].replace({"DC": "DC binding"})
        out["curtailed_mwh_by_limit"] = {k: round(float(v)) for k, v in
                                         lost.groupby(lim).sum().sort_values(ascending=False).head(6).items()}

        # 4. Volt-VAR: recompute every sunny hour (Volt-VAR can also make passing hours fail) ---
        vv_hours = sunny
        head_vv = head.copy()
        h_vv, why_vv = ac_headroom(ex, vv_hours, head_dc[vv_hours], bus, size, volt_var=True)
        head_vv[vv_hours] = h_vv
        np.savez_compressed(CACHE / "phase4_headroom_voltvar.npz", head=head_vv, reasons=np.array(why_vv),
                            hours=vv_hours)
        vv_worse = int(np.sum(np.minimum(h_vv, size * shape[vv_hours]) <
                              np.minimum(head[vv_hours], size * shape[vv_hours]) - 1e-6))
        vv_better = int(np.sum(np.minimum(h_vv, size * shape[vv_hours]) >
                               np.minimum(head[vv_hours], size * shape[vv_hours]) + 1e-6))
        c_vv, _, d_vv = curtailment(size, shape, head_vv)
        out["volt_var"] = {"hours_recomputed": int(len(vv_hours)), "hours_better": vv_better,
                           "hours_worse": vv_worse, "curtailment_pct": round(100 * c_vv, 2),
                           "delivered_mwh": round(d_vv), "extra_mwh_vs_unity": round(d_vv - deliv_flex),
                           "firm_mw": round(float(np.min(head_vv[shape >= 0.05] / shape[shape >= 0.05])), 1)}
        print(f"Volt-VAR done {time.perf_counter() - t0:.0f} s", flush=True)

        # 3. line upgrades: branches that limit the most curtailed energy -------------------
        lim_branches = [k.split()[1] for k in out["curtailed_mwh_by_limit"] if k.startswith("thermal ")]
        dc_lim = pd.Series(binding[short].sum(0), index=rts.branch.uid.values)
        lim_branches += [u for u in dc_lim[dc_lim > 0].sort_values(ascending=False).index if u not in lim_branches]
        upgrades = []
        prev = RES / "phase4_line_upgrades.csv"
        if "--reuse-upgrades" in sys.argv and prev.exists():      # unity-pf headroom is deterministic
            upgrades = pd.read_csv(prev).to_dict("records")       # same headroom -> same upgrade results
            lim_branches = []
        for uid in lim_branches[:2]:
            k = int(np.where(rts.branch.uid.values == uid)[0][0])
            br = rts.branch.iloc[k]
            kv = int(rts.bus.base_kv[rts.bus.id == br.f].iloc[0])
            redo = np.unique(np.r_[short, np.where(binding[:, k])[0]])
            for kind, factor, cost_mi in (("reconductor", RECONDUCTOR_UPRATE, RECONDUCTOR_USD_PER_MILE),
                                          ("rebuild", REBUILD_UPRATE, REBUILD_USD_PER_MILE)):
                r2 = rts.branch.rate_a.values.copy()
                r2[k] *= factor
                hd2, _ = dc_headroom(ex, redo, bus, rating=r2)
                sun_redo = redo[shape[redo] >= 0.01]
                ha2, _ = ac_headroom(ex, sun_redo, hd2[shape[redo] >= 0.01], bus, size, rating=r2)
                head2 = head.copy()
                head2[redo] = hd2
                head2[sun_redo] = ha2
                c2, _, d2 = curtailment(size, shape, head2)
                capex = cost_mi[kv] * max(br.length_mi, 1.0)
                avoided = d2 - deliv_flex
                upgrades.append({"branch": uid, "from": int(br.f), "to": int(br.t), "kv": kv,
                                 "miles": float(br.length_mi), "upgrade": kind, "new_rating_mva": round(r2[k]),
                                 "capex_musd": round(capex / 1e6, 2), "curtailment_pct": round(100 * c2, 2),
                                 "avoided_mwh": round(avoided),
                                 "usd_per_avoided_mwh": round(annualized(capex, RATE, YEARS) / avoided, 2)
                                 if avoided > 1 else None,
                                 "new_firm_mw": round(float(np.min(head2[shape >= 0.05] / shape[shape >= 0.05])), 1)})
                print(f"upgrade {uid} {kind}: curtailment {100 * c2:.2f} %", flush=True)
        up = pd.DataFrame(upgrades)
        up.to_csv(RES / "phase4_line_upgrades.csv", index=False)
        out["line_upgrades"] = up.to_dict("records")

    # 2. battery sweep ---------------------------------------------------------------------------
    pot_t = size * shape
    rows = []
    for mode in ("capture", "arbitrage"):           # capture: charge only would-be-curtailed solar
        for p_mw in (0, 25, 50, 100, 150, 200, 300):
            x = battery_dispatch(pot_t, head, lmp, p_mw, 4 * p_mw, rte=BATTERY_RTE, capture_only=mode == "capture")
            rows.append({"mode": mode, "battery_mw": p_mw, "battery_mwh": 4 * p_mw, "delivered_mwh": x.sum(),
                         "curtailment_pct": 100 * (1 - x.sum() / pot_t.sum()),
                         "revenue_usd": float((x * lmp).sum())})
    bat = pd.DataFrame(rows)
    base = bat[(bat.battery_mw == 0)].iloc[0]
    bat["avoided_mwh"] = bat.delivered_mwh - base.delivered_mwh
    bat["extra_revenue_usd"] = bat.revenue_usd - base.revenue_usd
    for c in BATTERY_USD_PER_KWH:
        ann = annualized(bat.battery_mwh * 1000 * c, RATE, YEARS, BATTERY_FOM)
        bat[f"usd_per_avoided_mwh_at_{c}"] = np.where(bat.avoided_mwh > 1, ann / bat.avoided_mwh.clip(lower=1), np.nan)
        bat[f"net_value_usd_per_yr_at_{c}"] = bat.extra_revenue_usd - ann
    # installed $/kWh at which the extra market revenue pays for the battery
    per_kwh_ann = annualized(1000.0, RATE, YEARS, BATTERY_FOM)          # annual cost of $1,000 per MWh
    bat["breakeven_usd_per_kwh"] = np.where(bat.battery_mwh > 0,
                                            bat.extra_revenue_usd / (bat.battery_mwh * per_kwh_ann), np.nan)
    bat.round(2).to_csv(RES / "phase4_battery_sweep.csv", index=False)
    out["battery"] = bat.round(2).to_dict("records")
    out["runtime_s"] = round(time.perf_counter() - t0, 1)
    out["assumptions"] = {"reconductor_usd_per_mile": RECONDUCTOR_USD_PER_MILE, "rebuild_usd_per_mile": REBUILD_USD_PER_MILE,
                          "reconductor_uprate": RECONDUCTOR_UPRATE, "rebuild_uprate": REBUILD_UPRATE,
                          "battery_usd_per_kwh": BATTERY_USD_PER_KWH, "battery_fom_frac": BATTERY_FOM,
                          "battery_rte": BATTERY_RTE, "discount_rate": RATE, "years": YEARS}
    (RES / "phase4_summary.json").write_text(json.dumps(out, indent=2, default=float))
    figures(bat, up, size, shape, head, head_vv, flex, firm)
    print(json.dumps(out, indent=2, default=float))


def figures(bat, up, size, shape, head, head_vv, flex, firm):
    from gridstudy.mitigation import curtailment
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    sizes = np.linspace(0, SWEEP_FACTOR * size, 80)          # headroom is exact only up to here
    sizes_vv = sizes[sizes <= size]                           # Volt-VAR headroom exact up to the project size
    axes[0].plot(sizes, [100 * curtailment(s, shape, head)[0] for s in sizes], label="unity pf")
    axes[0].plot(sizes_vv, [100 * curtailment(s, shape, head_vv)[0] for s in sizes_vv], "--", label="Volt-VAR")
    axes[0].axvline(firm, color="k", ls="--", lw=1, label=f"firm {firm:.0f} MW")
    axes[0].axvline(size, color="tab:red", ls=":", label=f"project {size:.0f} MW")
    for k, v in flex.items():
        axes[0].plot(v, float(k.rstrip("%")), "o", color="tab:green")
        axes[0].annotate(f"{k}: {v:.0f} MW", (v, float(k.rstrip("%"))), xytext=(4, -10),
                         textcoords="offset points", fontsize=8)
    axes[0].set_xlabel("plant size (MW)")
    axes[0].set_ylabel("annual curtailment (%)")
    axes[0].set_title("Flexible interconnection: size vs. curtailment")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)
    cap = bat[bat["mode"] == "capture"]
    arb = bat[bat["mode"] == "arbitrage"]
    axes[1].plot(cap.battery_mwh, cap.curtailment_pct, "o-", label="curtailment (capture mode)")
    axes[1].set_xlabel("battery energy (MWh, 4-hour)")
    axes[1].set_ylabel("annual curtailment (%)")
    ax2 = axes[1].twinx()
    ax2.plot(arb.battery_mwh, arb.extra_revenue_usd / 1e6, "s--", color="tab:orange", label="extra revenue (arbitrage)")
    ax2.set_ylabel("extra market revenue (M$/yr)")
    axes[1].set_title(f"Co-located battery, {size:.0f} MW plant")
    axes[1].grid(alpha=0.3)
    for c in (150, 350):
        axes[2].plot(cap.battery_mwh[1:], cap[f"usd_per_avoided_mwh_at_{c}"][1:], "o-", label=f"battery ${c}/kWh")
    for _, r in up.iterrows():
        if pd.notna(r.usd_per_avoided_mwh):
            axes[2].axhline(r.usd_per_avoided_mwh, ls="--", lw=1, label=f"{r.upgrade} {r.branch} (${r.capex_musd}M)",
                            color="tab:red" if r.upgrade == "reconductor" else "tab:purple")
    axes[2].set_yscale("log")
    axes[2].set_xlabel("battery energy (MWh)")
    axes[2].set_ylabel("annualized $ per MWh of curtailment avoided")
    axes[2].set_title("Battery vs. line upgrade")
    axes[2].legend(fontsize=7)
    axes[2].grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(FIG / "phase4_mitigations.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
