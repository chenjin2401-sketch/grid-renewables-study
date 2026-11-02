"""Phase 0: reproduce the 5-bus 24-hour study in Python and validate it against MATPOWER.

Writes results/phase0_*.csv, figures/phase0_*.png and prints the validation table.
Needs data/private/ (see data/README.md).
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

from gridstudy.networks import branch_labels, five_bus          # noqa: E402
from gridstudy.scenarios import economic_dispatch, fixed_schedule, load_profile  # noqa: E402
from gridstudy.validate import compare, load_reference          # noqa: E402

PRIVATE = ROOT / "data" / "private"
RES = ROOT / "results"
FIG = ROOT / "figures"


def save_hourly(res, gens, labels, tag):
    hours = np.arange(res.vm.shape[0])
    g = pd.DataFrame({"hour": hours})
    for k, n in enumerate(gens.name):
        g[f"P{n}_MW"] = res.pg[:, k]
    for k, n in enumerate(gens.name):
        g[f"Q{n}_MVAr"] = res.qg[:, k]
    g["loss_MW"] = res.loss
    g["cost_usd"] = res.cost
    g.to_csv(RES / f"phase0_{tag}_generation.csv", index=False)
    v = pd.DataFrame(np.c_[hours, res.vm, res.va],
                     columns=["hour"] + [f"V{b}_pu" for b in range(1, 6)] + [f"ang{b}_deg" for b in range(1, 6)])
    v.to_csv(RES / f"phase0_{tag}_voltages.csv", index=False)
    b = pd.DataFrame(np.c_[hours, res.p_from, res.p_to, res.s_max],
                     columns=["hour"] + [f"Pfrom_{l}" for l in labels] + [f"Pto_{l}" for l in labels]
                     + [f"Smax_{l}" for l in labels])
    b.to_csv(RES / f"phase0_{tag}_branch_flows.csv", index=False)


def plot_overlay(results, refs, labels):
    """Python (lines) vs. MATPOWER (markers) for voltages and branch flows, both scenarios."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    hours = np.arange(24)
    titles = {"A": "Scenario A - fixed schedule", "B": "Scenario B - economic dispatch"}
    for col, tag in enumerate(["A", "B"]):
        res, ref = results[tag], refs.get(tag)
        ax = axes[0, col]
        for b in range(5):
            line, = ax.plot(hours, res.vm[:, b], lw=1.6, label=f"Bus {b + 1}")
            if ref is not None:
                ax.plot(hours, ref["vm"][:, b], "o", ms=4, mfc="none", color=line.get_color())
        ax.axhline(0.95, ls="--", color="k", lw=1)
        ax.set_title(titles[tag])
        ax.set_ylabel("|V| (pu)")
        ax.grid(alpha=0.3)
        ax = axes[1, col]
        for k, lab in enumerate(labels):
            line, = ax.plot(hours, res.p_from[:, k], lw=1.6, label=f"Line {lab}")
            if ref is not None:
                ax.plot(hours, ref["p_from"][:, k], "o", ms=4, mfc="none", color=line.get_color())
        ax.set_ylabel("P from-end (MW)")
        ax.set_xlabel("Hour of the day")
        ax.grid(alpha=0.3)
    axes[0, 1].legend(loc="center left", bbox_to_anchor=(1.01, 0.5))
    axes[1, 1].legend(loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.suptitle("Phase 0 validation: pandapower (lines) vs. MATPOWER (circles)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIG / "phase0_validation_overlay.png", dpi=160)
    plt.close(fig)


def main():
    if not (PRIVATE / "load_profile_5bus.csv").exists():
        sys.exit("data/private/load_profile_5bus.csv missing - run scripts/import_course_reference.py first")
    RES.mkdir(exist_ok=True)
    FIG.mkdir(exist_ok=True)
    p_load, q_load = load_profile(PRIVATE / "load_profile_5bus.csv")
    labels = branch_labels(five_bus()[0])

    t0 = time.perf_counter()
    res_a, gens, _ = fixed_schedule(p_load, q_load)
    t1 = time.perf_counter()
    res_b, _, info_b = economic_dispatch(p_load, q_load)
    t2 = time.perf_counter()
    results = {"A": res_a, "B": res_b}
    save_hourly(res_a, gens, labels, "A_fixed")
    save_hourly(res_b, gens, labels, "B_dispatch")

    ref_dir = PRIVATE / "matpower_reference"
    refs, tables = {}, {}
    for tag, mtag in (("A", "task1"), ("B", "task2")):
        if (ref_dir / f"{mtag}_voltages.csv").exists():
            refs[tag] = load_reference(ref_dir, mtag)
            tables[tag] = compare(results[tag], refs[tag])
            tables[tag].to_csv(RES / f"phase0_validation_{tag}.csv", index=False)

    plot_overlay(results, refs, labels)

    summary = {
        "converged_hours": {k: int(r.converged.sum()) for k, r in results.items()},
        "runtime_s": {"A": round(t1 - t0, 3), "B": round(t2 - t1, 3)},
        "lp_vs_merit_order_max_mw": info_b["lp_vs_merit_max_mw"],
        "daily_cost_usd": {k: float(r.cost.sum()) for k, r in results.items()},
        "daily_loss_mwh": {k: float(r.loss.sum()) for k, r in results.items()},
        "min_voltage_pu": {k: float(r.vm.min()) for k, r in results.items()},
        "voltage_violation_bus_hours": {k: int((r.vm < 0.95).sum()) for k, r in results.items()},
        "max_branch_loading_mva": {k: float(r.s_max.max()) for k, r in results.items()},
    }
    (RES / "phase0_summary.json").write_text(json.dumps(summary, indent=2))

    pd.set_option("display.width", 140)
    print(json.dumps(summary, indent=2))
    for tag, t in tables.items():
        print(f"\nValidation vs. MATPOWER - scenario {tag}")
        print(t.to_string(index=False, float_format=lambda x: f"{x:.3e}"))


if __name__ == "__main__":
    main()
