"""Compare hourly results with MATPOWER reference CSVs."""
import numpy as np
import pandas as pd

LINES = ["1_2", "1_4", "1_5", "2_3", "3_4", "4_5"]


def load_reference(ref_dir, tag):
    """Read <tag>_generation/_voltages/_branch_flows.csv written by the MATLAB study."""
    g = pd.read_csv(ref_dir / f"{tag}_generation.csv")
    v = pd.read_csv(ref_dir / f"{tag}_voltages.csv")
    b = pd.read_csv(ref_dir / f"{tag}_branch_flows.csv")
    return {
        "vm": v[[f"V{i}_pu" for i in range(1, 6)]].to_numpy(),
        "va": v[[f"ang{i}_deg" for i in range(1, 6)]].to_numpy(),
        "p_from": b[[f"Pfrom_L{l}" for l in LINES]].to_numpy(),
        "p_to": -b[[f"Parrive_L{l}" for l in LINES]].to_numpy(),   # Parrive = -Pt
        "s_max": b[[f"Smax_L{l}" for l in LINES]].to_numpy(),
        "pg": g[[f"PG{i}_MW" for i in range(1, 6)]].to_numpy(),
        "qg": g[[f"QG{i}_MVAr" for i in range(1, 6)]].to_numpy(),
        "loss": g["loss_MW"].to_numpy(),
        "cost": g["cost_usd"].to_numpy(),
    }


QUANTITIES = [("vm", "Voltage magnitude", "pu"), ("va", "Voltage angle", "deg"),
              ("p_from", "Branch P (from end)", "MW"), ("p_to", "Branch P (to end)", "MW"),
              ("s_max", "Branch loading", "MVA"), ("pg", "Generator P", "MW"),
              ("qg", "Generator Q", "MVAr"), ("loss", "Total losses", "MW"),
              ("cost", "Hourly cost", "$")]


def compare(res, ref):
    """Max absolute and max relative deviation per quantity over all hours."""
    rows = []
    for key, name, unit in QUANTITIES:
        a = np.asarray(getattr(res, key), dtype=float)
        b = np.asarray(ref[key], dtype=float)
        diff = np.abs(a - b)
        scale = np.maximum(np.abs(b), 1e-6)
        big = np.abs(b) > 1.0 if unit not in ("pu", "deg") else np.ones_like(b, dtype=bool)
        rel = (diff / scale)[big].max() * 100 if big.any() else 0.0
        rows.append({"quantity": name, "unit": unit, "max_abs_diff": diff.max(),
                     "max_rel_diff_pct": rel, "reference_range": f"{b.min():.4g} .. {b.max():.4g}"})
    return pd.DataFrame(rows)
