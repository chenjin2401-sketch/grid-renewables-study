"""N-1 contingency analysis: DC screening with LODFs, AC verification, full-AC benchmark.

Criterion: after any single branch outage, no branch above its long-term emergency (LTE)
rating. Outages that island part of the network are reported separately and skipped.
"""
import time

import numpy as np
import pandas as pd

from .ac_hourly import set_hour, solve_hour
from .rts import lodf


def critical_hours(ts, y, n_each=12, dc_limit=None):
    """Peak-load, lowest-net-load (sunny/light load) and most-stressed (max DC loading) hours."""
    load = ts.pd_bus.sum(1)
    net = load - ts.rtpv_bus.sum(1) - np.nansum(ts.avail, 1)
    flow = np.abs(y["flow"])
    stress = (flow / dc_limit).max(1) if dc_limit is not None else flow.max(1)   # % of DC limit
    picks = np.r_[np.argsort(-load)[:n_each], np.argsort(net)[:n_each], np.argsort(-stress)[:n_each]]
    return np.unique(picks)


def project_dispatch(rts, ts, y, hours, project, opf=None):
    """Hourly dispatch and DC flows with a new solar project.

    The existing fleet is redispatched by the DC-OPF around the plant's output (same
    convention as the hosting-capacity study). If the network cannot take the full output in an
    hour, the plant is curtailed to its headroom first. Returns pg, flows, plant MW.
    """
    pg = y["pg"][hours].copy()
    flow = y["flow"][hours].copy()
    p_new = np.zeros(len(hours))
    if project is not None:
        from .dcopf import DCOPF
        from .hosting import redispatch_with_plant
        opf = opf or DCOPF(rts)
        pos = {b: i for i, b in enumerate(rts.bus.id)}
        shape = ts.pv_shape_area[project["bus"] // 100][hours]
        for j, h in enumerate(hours):
            p = project["mw"] * shape[j]
            r = redispatch_with_plant(opf, ts, h, project["bus"], p)
            if r.status != 0:                                  # curtail to the network headroom
                p = opf.max_injection(ts.pd_bus[h] - ts.rtpv_bus[h], ts.avail[h], pos[project["bus"]]) * 0.999
                r = redispatch_with_plant(opf, ts, h, project["bus"], p)
            pg[j], flow[j], p_new[j] = r.pg, r.flow, p
    return pg, flow, p_new


def dc_screen(rts, flow, threshold=1.0, q_base=None):
    """Post-contingency DC flows for every (hour, outage, monitored branch).

    q_base: optional (hours, branches) reactive flow (MVAr) from the intact-system AC power flow.
    Without it: loading = |P_post| / (LTE / MVA factor) (pure DC screen).
    With it ("hybrid" screen): loading = sqrt(P_post^2 + Q_base^2) / LTE, which also sees
    branches that are loaded by reactive flow (invisible to a DC model).
    Returns DataFrame of (hour_pos, outage, monitored, dc_loading) above threshold, and the
    islanding mask.
    """
    L, island = lodf(rts.ptdf, rts.branch, rts.bus)
    lte = rts.branch.rate_lte.values
    if q_base is None and rts.mva_factor is not None:       # same MW-for-MVA margin as the DC-OPF
        lte = lte / rts.mva_factor
    out = []
    for j, f in enumerate(flow):
        post = f[:, None] + L * f[None, :]                 # (monitored, outage)
        if q_base is None:
            load = np.abs(post) / lte[:, None]
        else:
            load = np.hypot(post, q_base[j][:, None]) / lte[:, None]
        load[:, island] = 0
        np.fill_diagonal(load, 0)
        m, k = np.where(load > threshold)
        out += [(j, kk, mm, load[mm, kk]) for mm, kk in zip(m, k)]
    return pd.DataFrame(out, columns=["hour_pos", "outage", "monitored", "dc_loading"]), island


def ac_base_case(rts, ac, ts, hours, pg, p_new, sgen=None):
    """Intact-system AC power flow per hour: (converged mask, max |Q| per branch end in MVAr)."""
    from .rts import branch_flows
    conv = np.zeros(len(hours), bool)
    q = np.zeros((len(hours), len(rts.branch)))
    for j, h in enumerate(hours):
        extra = {sgen: (p_new[j], 0.0)} if sgen is not None else None
        set_hour(ac, rts, ts.pd_bus[h], ts.qd_bus[h], ts.rtpv_bus[h], pg[j], extra_sgen=extra)
        a = solve_hour(ac, rts)
        conv[j] = a.converged
        if a.converged:
            _, qf, _, qt = branch_flows(ac)
            q[j] = np.maximum(np.abs(qf), np.abs(qt))
    return conv, q


def _set_outage(ac, k, state):
    et, ei = ac.branch_elem[k]
    getattr(ac.net, et).at[ei, "in_service"] = state


def ac_contingencies(rts, ac, ts, y, hours, pg, p_new, pairs=None, sgen=None):
    """AC N-1. pairs: iterable of (hour_pos, outage) to run; None = every non-islanding outage.

    Returns DataFrame (hour_pos, outage, monitored, ac_loading) of LTE overloads, the number of
    AC power flows solved, the non-converged cases and the runtime.
    """
    _, island = lodf(rts.ptdf, rts.branch, rts.bus)
    lte = rts.branch.rate_lte.values
    todo = {}
    if pairs is None:
        for j in range(len(hours)):
            todo[j] = [k for k in range(len(rts.branch)) if not island[k]]
    else:
        for j, k in pairs:
            todo.setdefault(j, []).append(k)
    res, n_pf, failed = [], 0, []
    t0 = time.perf_counter()
    for j, ks in todo.items():
        h = hours[j]
        extra = {sgen: (p_new[j], 0.0)} if sgen is not None else None
        set_hour(ac, rts, ts.pd_bus[h], ts.qd_bus[h], ts.rtpv_bus[h], pg[j], extra_sgen=extra)
        for k in sorted(set(ks)):
            _set_outage(ac, k, False)
            a = solve_hour(ac, rts)
            _set_outage(ac, k, True)
            n_pf += 1
            if not a.converged:
                failed.append((j, k))
                continue
            load = a.s_max / lte
            load[k] = 0
            for m in np.where(load > 1.0)[0]:
                res.append((j, k, m, load[m]))
    df = pd.DataFrame(res, columns=["hour_pos", "outage", "monitored", "ac_loading"])
    return df, n_pf, failed, time.perf_counter() - t0
