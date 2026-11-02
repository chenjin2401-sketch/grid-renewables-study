"""Hour-by-hour AC power flow."""
from dataclasses import dataclass

import numpy as np
import pandapower as pp


@dataclass
class HourlyResults:
    vm: np.ndarray        # (hours, buses) voltage magnitude, pu
    va: np.ndarray        # (hours, buses) voltage angle, deg
    p_from: np.ndarray    # (hours, lines) MW injected at the from end
    p_to: np.ndarray      # (hours, lines) MW injected at the to end
    q_from: np.ndarray
    q_to: np.ndarray
    s_max: np.ndarray     # (hours, lines) MVA loading, larger of the two ends
    pg: np.ndarray        # (hours, gens) MW (slack from the power flow)
    qg: np.ndarray        # (hours, gens) MVAr
    loss: np.ndarray      # (hours,) total series I^2R loss, MW
    cost: np.ndarray      # (hours,) $ = sum(price * pg)
    converged: np.ndarray


def run_hourly(net, gens, pd, qd, pg, qg):
    """Solve one AC power flow per row of the schedules.

    pd, qd : (hours, buses) load in MW / MVAr
    pg, qg : (hours, gens) scheduled generator output; the slack entries are ignored
    """
    nh, nb, nl, ng = pd.shape[0], len(net.bus), len(net.line), len(gens.name)
    out = {k: np.zeros((nh, n)) for k, n in
           [("vm", nb), ("va", nb), ("p_from", nl), ("p_to", nl), ("q_from", nl),
            ("q_to", nl), ("pg", ng), ("qg", ng)]}
    loss = np.zeros(nh)
    conv = np.zeros(nh, dtype=bool)
    price = np.asarray(gens.price, dtype=float)

    for h in range(nh):
        net.load["p_mw"] = pd[h]
        net.load["q_mvar"] = qd[h]
        for k, (kind, idx) in enumerate(gens.element):
            if kind == "sgen":
                net.sgen.at[idx, "p_mw"] = pg[h, k]
                net.sgen.at[idx, "q_mvar"] = qg[h, k]
        pp.runpp(net, algorithm="nr", calculate_voltage_angles=True, init="flat",
                 tolerance_mva=1e-9, max_iteration=30, numba=False)
        conv[h] = net.converged
        out["vm"][h] = net.res_bus.vm_pu.values
        out["va"][h] = net.res_bus.va_degree.values
        out["p_from"][h] = net.res_line.p_from_mw.values
        out["p_to"][h] = net.res_line.p_to_mw.values
        out["q_from"][h] = net.res_line.q_from_mvar.values
        out["q_to"][h] = net.res_line.q_to_mvar.values
        loss[h] = net.res_line.pl_mw.sum()
        for k, (kind, idx) in enumerate(gens.element):
            res = net.res_sgen if kind == "sgen" else net.res_ext_grid
            out["pg"][h, k] = res.at[idx, "p_mw"]
            out["qg"][h, k] = res.at[idx, "q_mvar"]

    s_from = np.hypot(out["p_from"], out["q_from"])
    s_to = np.hypot(out["p_to"], out["q_to"])
    return HourlyResults(**out, s_max=np.maximum(s_from, s_to), loss=loss,
                         cost=out["pg"] @ price, converged=conv)
