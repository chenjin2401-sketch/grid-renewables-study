"""Phase 1: RTS-GMLC model validation and full-year hourly baseline (DC-OPF + AC power flow).

Outputs
  cache/phase1.npz                    hourly arrays (git-ignored, reused by later phases)
  results/phase1_validation.csv       model vs. MATPOWER reference
  results/phase1_summary.json         headline numbers
  results/phase1_branch_stats.csv     congestion / loading statistics per branch
  figures/phase1_*.png
"""
import json
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gridstudy.ac_hourly import V_LIMITS, set_hour, solve_hour   # noqa: E402
from gridstudy.dcopf import DCOPF                                  # noqa: E402
from gridstudy.rts import build_ac, load_rts, load_timeseries     # noqa: E402
from gridstudy.validate_rts import base_case_ac, reference_ac_pf, reference_dc_opf  # noqa: E402

CACHE, RES, FIG = ROOT / "cache", ROOT / "results", ROOT / "figures"


def validate(rts):
    vm_ref, va_ref, pf_ref, loss_ref = reference_ac_pf()
    vm, va, pf, loss = base_case_ac(rts, build_ac)
    ids = list(rts.bus.id)
    rows = [
        {"check": "AC PF |V| (pu)", "max_abs_diff": np.abs(vm - [vm_ref[b] for b in ids]).max(),
         "note": "reference printed to 3 decimals"},
        {"check": "AC PF angle (deg)", "max_abs_diff": np.abs(va - [va_ref[b] for b in ids]).max(),
         "note": "reference printed to 3 decimals"},
        {"check": "AC PF branch P from (MW)", "max_abs_diff": np.abs(pf - pf_ref).max(),
         "note": "reference printed to 2 decimals"},
        {"check": "AC PF total losses (MW)", "max_abs_diff": abs(loss - loss_ref), "note": f"reference {loss_ref}"},
    ]
    obj_ref, lam_ref = reference_dc_opf()
    r = DCOPF(rts, mode="matpower").solve(rts.bus.pd.values)
    rows += [{"check": "DC-OPF objective ($/h)", "max_abs_diff": abs(r.cost - obj_ref), "note": f"reference {obj_ref}"},
             {"check": "DC-OPF LMP ($/MWh)", "max_abs_diff": np.abs(r.lmp - lam_ref).max(),
              "note": f"reference {lam_ref} at every bus"}]
    return pd.DataFrame(rows)


def calibrate_mva_factors(rts, ts, n_hours=600, cap=1.10):
    """AC MVA / DC MW ratio of each branch on the most stressed hours of an uncorrected DC-OPF run.

    The DC-OPF limits MW, the rating is in MVA: a branch dispatched to exactly its MW limit is
    over its MVA rating once reactive flow and losses are included. Each DC limit is therefore
    reduced by the 99th-percentile ratio observed when the branch is heavily loaded (>= 70 %),
    capped at 1.10; branches needing more (reactive-driven) are left as baseline findings.
    """
    opf = DCOPF(rts, mva_factor=np.ones(len(rts.branch)))
    nh = ts.pd_bus.shape[0]
    flows = np.zeros((nh, len(rts.branch)))
    for h in range(nh):
        flows[h] = opf.solve(ts.pd_bus[h] - ts.rtpv_bus[h], avail=ts.avail[h]).flow
    rate = rts.branch.rate_a.values
    hours = np.argsort(-(np.abs(flows) / rate).max(1))[:n_hours]
    y = run_year(rts, ts, hours=hours, mva_factor=np.ones(len(rate)))
    p = np.abs(y["flow"])
    heavy = (p > 0.7 * rate) & y["conv"][:, None]
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(heavy, y["s_ac"] / np.maximum(p, 1.0), np.nan)
    raw = np.array([np.nanpercentile(ratio[:, k], 99) if heavy[:, k].sum() >= 10 else 1.0
                    for k in range(len(rate))])
    raw = np.nan_to_num(raw, nan=1.0)
    f = np.clip(raw, 1.0, cap)
    out = pd.DataFrame({"uid": rts.branch.uid, "mva_factor": np.round(f, 4), "observed_ratio": np.round(raw, 4),
                        "heavy_hours": heavy.sum(0), "capped": raw > cap})
    out.to_csv(RES / "phase1_mva_factors.csv", index=False)
    return f, hours, y


def run_year(rts, ts, hours=None, mva_factor=None):
    hours = np.arange(ts.pd_bus.shape[0]) if hours is None else np.asarray(hours)
    nh, nb, nl, ng = len(hours), len(rts.bus), len(rts.branch), len(rts.gen)
    opf = DCOPF(rts, mva_factor=mva_factor)
    ac = build_ac(rts)
    out = {k: np.zeros(s) for k, s in [("pg", (nh, ng)), ("flow", (nh, nl)), ("lmp", (nh, nb)), ("cost", nh),
                                       ("shed", nh), ("vm", (nh, nb)), ("s_ac", (nh, nl)), ("p_ac", (nh, nl)),
                                       ("loss_ac", nh)]}
    out["conv"] = np.zeros(nh, bool)
    t0 = time.perf_counter()
    for j, h in enumerate(hours):
        r = opf.solve(ts.pd_bus[h] - ts.rtpv_bus[h], avail=ts.avail[h])
        out["pg"][j], out["flow"][j], out["lmp"][j], out["cost"][j], out["shed"][j] = \
            r.pg, r.flow, r.lmp, r.cost, r.shed
        set_hour(ac, rts, ts.pd_bus[h], ts.qd_bus[h], ts.rtpv_bus[h], r.pg)
        a = solve_hour(ac, rts)
        out["conv"][j], out["vm"][j], out["s_ac"][j], out["p_ac"][j], out["loss_ac"][j] = \
            a.converged, a.vm, a.s_max, a.p_from, a.loss
        if j % 1000 == 0:
            print(f"hour {j}/{nh}  {time.perf_counter() - t0:.0f} s", flush=True)
    out["runtime_s"] = time.perf_counter() - t0
    return out


def figures(rts, ts, y):
    FIG.mkdir(exist_ok=True)
    idx = ts.index
    # 1. energy mix by week
    fuel = rts.gen.type.replace({"CT": "Gas/oil CT", "CC": "Gas CC", "STEAM": "Coal/oil steam"}).values
    groups = ["NUCLEAR", "Coal/oil steam", "Gas CC", "Gas/oil CT", "HYDRO", "WIND", "PV"]
    weekly = pd.DataFrame({g: y["pg"][:, fuel == g].sum(1) for g in groups}, index=idx)
    weekly["Rooftop PV"] = ts.rtpv_bus.sum(1)
    weekly = weekly.resample("W").mean()
    fig, ax = plt.subplots(figsize=(12, 4.8))
    ax.stackplot(weekly.index, weekly.T.values, labels=weekly.columns, alpha=0.9)
    ax.plot(idx, ts.pd_bus.sum(1), color="k", lw=0.3, alpha=0.4, label="Hourly load")
    ax.set_ylabel("MW (weekly average)")
    ax.set_title("RTS-GMLC 2020 baseline dispatch (hourly DC-OPF), weekly average by resource")
    ax.legend(loc="upper left", ncol=4, fontsize=8)
    ax.margins(x=0)
    fig.tight_layout()
    fig.savefig(FIG / "phase1_energy_mix.png", dpi=150)
    plt.close(fig)
    # 2. LMP duration + voltage / loading distributions
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    lmp_sys = np.sort(y["lmp"].mean(1))[::-1]
    axes[0].plot(np.arange(len(lmp_sys)) / len(lmp_sys) * 100, lmp_sys)
    axes[0].set_xlabel("% of hours")
    axes[0].set_ylabel("Average LMP ($/MWh)")
    axes[0].set_title("LMP duration curve")
    axes[0].grid(alpha=0.3)
    vm = y["vm"][y["conv"]]
    axes[1].hist(vm.min(1), bins=60, alpha=0.7, label="lowest bus")
    axes[1].hist(vm.max(1), bins=60, alpha=0.7, label="highest bus")
    for v in V_LIMITS:
        axes[1].axvline(v, color="k", ls="--")
    axes[1].set_xlabel("|V| (pu)")
    axes[1].set_ylabel("hours")
    axes[1].set_title("Hourly voltage extremes (AC)")
    axes[1].legend()
    load = y["s_ac"][y["conv"]] / rts.branch.rate_a.values
    axes[2].plot(np.sort(load.max(1))[::-1] * 100)
    axes[2].axhline(100, color="r", ls="--")
    axes[2].set_xlabel("hours (sorted)")
    axes[2].set_ylabel("most loaded branch (% of normal rating)")
    axes[2].set_title("Branch loading duration (AC)")
    axes[2].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIG / "phase1_lmp_voltage_loading.png", dpi=150)
    plt.close(fig)


def main():
    for d in (CACHE, RES, FIG):
        d.mkdir(exist_ok=True)
    rts = load_rts(use_calibration=False)          # recalibrated below
    val = validate(rts)
    val.to_csv(RES / "phase1_validation.csv", index=False)
    print(val.to_string(index=False))
    ts = load_timeseries(rts)
    t_cal = time.perf_counter()
    f, cal_hours, y_cal = calibrate_mva_factors(rts, ts)
    t_cal = time.perf_counter() - t_cal
    rts.mva_factor = f
    y = run_year(rts, ts, mva_factor=f)
    np.savez_compressed(CACHE / "phase1.npz", **{k: v for k, v in y.items() if k != "runtime_s"})

    conv = y["conv"]
    rate = rts.branch.rate_a.values
    # effect of the MVA calibration on the same stressed hours
    over_before = ((y_cal["s_ac"] > rate) & y_cal["conv"][:, None]).any(1).sum()
    over_after = ((y["s_ac"][cal_hours] > rate) & conv[cal_hours][:, None]).any(1).sum()
    var = rts.gen.type.isin(["PV", "WIND", "HYDRO"]).values
    avail = np.nan_to_num(ts.avail[:, var]).sum()
    used = y["pg"][:, var].sum()
    over_ac = (y["s_ac"] > rate * 1.0) & conv[:, None]
    binding_dc = np.isclose(np.abs(y["flow"]), rate / f, atol=1e-3)
    vhi = (y["vm"] > V_LIMITS[1]) & conv[:, None]
    vlo = (y["vm"] < V_LIMITS[0]) & conv[:, None]
    lmp_mean = y["lmp"].mean(1)
    summary = {
        "hours": int(len(conv)), "ac_converged_hours": int(conv.sum()),
        "runtime_s_total": round(y["runtime_s"], 1),
        "annual_load_twh": round(ts.pd_bus.sum() / 1e6, 3),
        "peak_load_mw": round(ts.pd_bus.sum(1).max(), 1),
        "annual_cost_musd": round(y["cost"].sum() / 1e6, 2),
        "avg_lmp_usd_mwh": round(lmp_mean.mean(), 2),
        "lmp_p95_usd_mwh": round(np.percentile(lmp_mean, 95), 2),
        "hours_with_congestion_dc": int(binding_dc.any(1).sum()),
        "hours_with_lmp_spread_gt_1usd": int((y["lmp"].max(1) - y["lmp"].min(1) > 1).sum()),
        "renewable_curtailment_pct": round(100 * (1 - used / avail), 2),
        "renewable_share_of_load_pct": round(100 * (used + ts.rtpv_bus.sum()) / ts.pd_bus.sum(), 1),
        "unserved_energy_mwh": round(y["shed"].sum(), 1),
        "ac_loss_pct_of_load": round(100 * y["loss_ac"][conv].sum() / ts.pd_bus[conv].sum(), 2),
        "ac_hours_any_branch_over_100pct": int(over_ac.any(1).sum()),
        "ac_branches_ever_over_100pct": int(over_ac.any(0).sum()),
        "ac_max_branch_loading_pct": round(100 * np.nanmax(y["s_ac"][conv] / rate), 1),
        "ac_hours_any_v_above_1p05": int(vhi.any(1).sum()),
        "ac_hours_any_v_below_0p95": int(vlo.any(1).sum()),
        "ac_v_range_pu": [round(np.nanmin(y["vm"][conv]), 4), round(np.nanmax(y["vm"][conv]), 4)],
        "max_dc_vs_ac_branch_mw": round(np.nanmax(np.abs(y["flow"][conv] - y["p_ac"][conv])), 1),
        "mva_calibration": {"stressed_hours": int(len(cal_hours)), "runtime_s": round(t_cal, 1),
                            "branches_derated": int((f > 1.0001).sum()), "max_factor": round(float(f.max()), 3),
                            "ac_overload_hours_before": int(over_before), "ac_overload_hours_after": int(over_after)},
    }
    (RES / "phase1_summary.json").write_text(json.dumps(summary, indent=2))
    stats = pd.DataFrame({
        "uid": rts.branch.uid, "from": rts.branch.f, "to": rts.branch.t, "rate_a": rate,
        "dc_binding_hours": binding_dc.sum(0), "ac_over_100_hours": over_ac.sum(0),
        "ac_max_loading_pct": np.round(100 * np.nanmax(y["s_ac"][conv] / rate, axis=0), 1)})
    stats.sort_values("ac_max_loading_pct", ascending=False).to_csv(RES / "phase1_branch_stats.csv", index=False)
    figures(rts, ts, y)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
