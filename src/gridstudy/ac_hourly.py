"""AC power flow of one hour on the RTS-GMLC model, given a DC-OPF dispatch."""
from dataclasses import dataclass

import numpy as np
import pandapower as pp

from .rts import branch_flows

V_LIMITS = (0.95, 1.05)


@dataclass
class ACHour:
    converged: bool
    vm: np.ndarray            # (buses,)
    p_from: np.ndarray        # (branches,)
    s_max: np.ndarray         # (branches,) MVA, larger end
    loss: float               # MW
    q_new: float = 0.0        # MVAr from the new solar plant (Volt-VAR studies)


def set_hour(ac, rts, pd_bus, qd_bus, rtpv_bus, pg, extra_sgen=None):
    """Write one hour's load and dispatch into the pandapower model.

    pg: (gens,) MW from the DC-OPF. Conventional units at a bus are aggregated into one
    voltage-controlled plant. Plants keep regulating voltage even when the dispatch gives
    them 0 MW (units held online / condenser mode for voltage support; without this the
    lightly loaded 138/230 kV network reaches 1.2 pu from line charging). Slack is
    distributed over the plants in proportion to their output.
    """
    net = ac.net
    load_bus = net.load.bus.values
    pos = {bid: i for i, bid in enumerate(rts.bus.id)}
    idx = np.array([pos[b] for b in load_bus])
    net.load["p_mw"] = pd_bus[idx]
    net.load["q_mvar"] = qd_bus[idx]
    for b, si in ac.rtpv_sgen.items():
        net.sgen.at[si, "p_mw"] = rtpv_bus[pos[b]]
    for gi, si in ac.var_sgen.items():
        net.sgen.at[si, "p_mw"] = pg[gi] if np.isfinite(pg[gi]) else 0.0
    for b, gidx in ac.ctrl_gen.items():
        rows = ac.gens_at_bus[b]
        p = float(np.nansum(pg[rows]))
        net.gen.at[gidx, "p_mw"] = p
        net.gen.at[gidx, "slack_weight"] = max(p, 0.0)
    eg = net.ext_grid.index[0]
    net.ext_grid.at[eg, "slack_weight"] = max(float(np.nansum(pg[ac.slack_gens])), 1.0)
    if extra_sgen is not None:
        for si, (p, q) in extra_sgen.items():
            net.sgen.at[si, "p_mw"] = p
            net.sgen.at[si, "q_mvar"] = q


def solve_hour(ac, rts, enforce_q_lims=True):
    net = ac.net
    try:
        pp.runpp(net, algorithm="nr", calculate_voltage_angles=True, init="auto", distributed_slack=True,
                 enforce_q_lims=enforce_q_lims, max_iteration=30, tolerance_mva=1e-6)
    except pp.LoadflowNotConverged:
        nb, nl = len(rts.bus), len(rts.branch)
        return ACHour(False, np.full(nb, np.nan), np.full(nl, np.nan), np.full(nl, np.nan), np.nan)
    pf, qf, pt, qt = branch_flows(ac)
    vm = net.res_bus.loc[list(rts.bus.id), "vm_pu"].values
    loss = net.res_line.pl_mw.sum() + net.res_trafo.pl_mw.sum()
    return ACHour(True, vm, pf, np.maximum(np.hypot(pf, qf), np.hypot(pt, qt)), loss)
