"""Phase 0 scenarios on the 5-bus system: fixed schedule vs. economic dispatch."""
import numpy as np
import pandas as pd

from .dispatch import economic_dispatch_lp, merit_order
from .networks import five_bus
from .timeseries import run_hourly

PF = 0.95
Q_PER_P = np.sqrt(1 - PF ** 2) / PF          # Q = P * tan(acos(pf)) = 0.3287 P
SLACK_CAP_MW = 300.0                         # slack (hydro) capped in dispatch to hold reserve


def load_profile(path):
    """24 x 5 real-power load profile (MW) -> (P, Q) arrays, Q at pf 0.95."""
    df = pd.read_csv(path)
    p = df[[f"bus{b}" for b in range(1, 6)]].to_numpy(dtype=float)
    return p, p * Q_PER_P


def fixed_schedule(p_load, q_load):
    """Scenario A: non-slack units hold their specified MW with Qg = 0; the slack balances."""
    net, gens = five_bus()
    nh = p_load.shape[0]
    pg = np.tile(np.asarray(gens.p_init_mw, dtype=float), (nh, 1))
    qg = np.zeros_like(pg)
    res = run_hourly(net, gens, p_load, q_load, pg, qg)
    return res, gens, {"pg_sched": pg, "qg_sched": qg}


def economic_dispatch(p_load, q_load):
    """Scenario B: hourly lossless LP dispatch (slack capped at 300 MW), Q at pf 0.95,
    then AC power flow with the slack unit covering the losses."""
    net, gens = five_bus()
    price = np.asarray(gens.price, dtype=float)
    p_max = np.asarray(gens.p_max_mw, dtype=float)
    p_max[gens.slack] = SLACK_CAP_MW
    p_min = np.zeros_like(p_max)
    nh = p_load.shape[0]
    pg = np.zeros((nh, len(price)))
    lam = np.zeros(nh)
    max_diff = 0.0
    for h in range(nh):
        demand = p_load[h].sum()
        pg[h], lam[h] = economic_dispatch_lp(demand, price, p_min, p_max)
        max_diff = max(max_diff, np.abs(pg[h] - merit_order(demand, price, p_min, p_max)).max())
    qg = pg * Q_PER_P
    res = run_hourly(net, gens, p_load, q_load, pg, qg)
    return res, gens, {"pg_sched": pg, "qg_sched": qg, "lambda": lam,
                       "lp_vs_merit_max_mw": max_diff}
