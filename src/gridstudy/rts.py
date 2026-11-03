"""RTS-GMLC (73-bus) system: network, generator data, DC sensitivities and hourly time series.

Data: https://github.com/GridMod/RTS-GMLC (downloaded by scripts/download_rts_gmlc.py).

Modeling choices (documented in docs/phase1_results.md):
- Network from the MATPOWER case. Transformer taps are placed on the MATPOWER
  "from" side (pandapower's converter puts them on the HV side, which shifts
  voltages by up to 0.037 pu; with the fix the AC power flow matches MATPOWER).
- Ratings: continuous (normal) and long-term emergency (N-1) from branch.csv.
- Thermal units: no unit commitment; each unit can run from 0 to Pmax with the
  convex hull of its piecewise-linear cost curve through (0, 0).
- PV, wind and hydro: zero-cost, curtailable up to their hourly available MW.
  Rooftop PV (RTPV): fixed injection (behind-the-meter, netted from load).
- Not modeled: CSP plant, the existing 50 MW storage unit, the HVDC link (0 MW).
"""
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pandapower as pp
from pandapower.converter.pypower import from_ppc

from .matpower_io import read_case

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "rts_gmlc"
CASE = DATA / "FormattedData" / "MATPOWER" / "RTS_GMLC.m"
SLACK_BUS = 113
THERMAL = {"CT", "CC", "STEAM", "NUCLEAR"}
VARIABLE = {"PV", "WIND", "HYDRO"}
# units that regulate their bus voltage in the AC model (hydro has Q capability in RTS-GMLC;
# utility PV and wind have Qmin = Qmax = 0 in the data set)
VOLTAGE_CONTROL = THERMAL | {"HYDRO", "SYNC_COND"}


@dataclass
class RTS:
    ppc: dict
    bus: pd.DataFrame       # index 0..72: id, area, pd, qd, base_kv, vmax, vmin
    gen: pd.DataFrame       # index 0..157: uid, bus, type, fuel, pmin, pmax, qmin, qmax, vg, status
    branch: pd.DataFrame    # index 0..119: uid, f, t, r, x, b, ratio, rate_a, rate_lte, rate_ste, length_mi
    ptdf: np.ndarray        # (branches, buses) DC power transfer distribution factors, slack = bus 113
    cost_segments: list     # per gen: list of (width_mw, slope_$per_mwh) of the convexified cost from 0
    mva_factor: np.ndarray = None   # per branch AC MVA / DC MW margin used by the DC-OPF (Phase 1 calibration)


# ----------------------------------------------------------------------------- loading
def load_rts(case_path=CASE, use_calibration=True):
    ppc = read_case(case_path)
    b = ppc["bus"]
    bus = pd.DataFrame({"id": b[:, 0].astype(int), "type": b[:, 1].astype(int), "pd": b[:, 2], "qd": b[:, 3],
                        "gs": b[:, 4], "bs": b[:, 5], "area": b[:, 6].astype(int), "base_kv": b[:, 9],
                        "vmax": b[:, 11], "vmin": b[:, 12]})
    g = ppc["gen"]
    names = ppc["gen_name"]
    gen = pd.DataFrame({"uid": [n[0] for n in names], "type": [n[1] for n in names], "fuel": [n[2] for n in names],
                        "bus": g[:, 0].astype(int), "pg0": g[:, 1], "qmax": g[:, 3], "qmin": g[:, 4], "vg": g[:, 5],
                        "status": g[:, 7].astype(int), "pmax": g[:, 8], "pmin": g[:, 9]})
    br = ppc["branch"]
    csv = pd.read_csv(DATA / "SourceData" / "branch.csv")
    assert (csv["From Bus"].values == br[:, 0]).all() and (csv["To Bus"].values == br[:, 1]).all()
    branch = pd.DataFrame({"uid": csv["UID"], "f": br[:, 0].astype(int), "t": br[:, 1].astype(int),
                           "r": br[:, 2], "x": br[:, 3], "b": br[:, 4], "ratio": br[:, 8],
                           "rate_a": csv["Cont Rating"].astype(float), "rate_lte": csv["LTE Rating"].astype(float),
                           "rate_ste": csv["STE Rating"].astype(float), "length_mi": csv["Length"].astype(float)})
    factors = ROOT / "results" / "phase1_mva_factors.csv"
    mva = None
    if use_calibration and factors.exists():
        f = pd.read_csv(factors)
        assert (f.uid.values == branch.uid.values).all()
        mva = f.mva_factor.values
    return RTS(ppc=ppc, bus=bus, gen=gen, branch=branch, ptdf=make_ptdf(bus, branch),
               cost_segments=[_convex_segments(row) for row in ppc["gencost"]], mva_factor=mva)


def _convex_segments(gc):
    """Convex hull of the PWL cost points and the origin -> [(width, slope), ...] from 0 MW."""
    n = int(gc[3])
    pts = [(0.0, 0.0)] + [(gc[4 + 2 * k], gc[5 + 2 * k]) for k in range(n)]
    pts = sorted(set(pts))
    hull = []
    for p in pts:                                   # lower hull (monotone chain)
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            if (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1) <= 0:
                hull.pop()
            else:
                break
        hull.append(p)
    return [(x2 - x1, (y2 - y1) / (x2 - x1)) for (x1, y1), (x2, y2) in zip(hull[:-1], hull[1:]) if x2 > x1]


def make_ptdf(bus, branch):
    """DC PTDF matrix (MATPOWER makeBdc convention: b = 1 / (x * tap)), slack = SLACK_BUS."""
    nb, nl = len(bus), len(branch)
    pos = {bid: i for i, bid in enumerate(bus.id)}
    tap = np.where(branch.ratio.values == 0, 1.0, branch.ratio.values)
    bsr = 1.0 / (branch.x.values * tap)
    cf = np.zeros((nl, nb))
    fi = branch.f.map(pos).values
    ti = branch.t.map(pos).values
    cf[np.arange(nl), fi] = 1
    cf[np.arange(nl), ti] = -1
    bf = bsr[:, None] * cf
    bbus = cf.T @ bf
    s = pos[SLACK_BUS]
    keep = [i for i in range(nb) if i != s]
    ptdf = np.zeros((nl, nb))
    ptdf[:, keep] = bf[:, keep] @ np.linalg.inv(bbus[np.ix_(keep, keep)])
    return ptdf


def lodf(ptdf, branch, bus):
    """Line outage distribution factors: LODF[l, k] = change of flow on l per MW on k when k trips."""
    pos = {bid: i for i, bid in enumerate(bus.id)}
    fi = branch.f.map(pos).values
    ti = branch.t.map(pos).values
    ptdf_kk = ptdf[:, fi] - ptdf[:, ti]           # (l, k): flow on l for 1 MW injected f_k -> t_k
    denom = 1 - np.diag(ptdf_kk)
    islanding = np.abs(denom) < 1e-6
    out = ptdf_kk / np.where(islanding, np.nan, denom)[None, :]
    np.fill_diagonal(out, -1.0)
    return out, islanding


# ----------------------------------------------------------------------------- time series
def _read_ts(rel):
    df = pd.read_csv(DATA / "timeseries_data_files" / rel)
    idx = pd.to_datetime(dict(year=df.Year, month=df.Month, day=df.Day)) + pd.to_timedelta(df.Period - 1, unit="h")
    return df.drop(columns=["Year", "Month", "Day", "Period"]).set_index(idx)


@dataclass
class TimeSeries:
    index: pd.DatetimeIndex
    pd_bus: np.ndarray      # (hours, buses) MW load
    qd_bus: np.ndarray      # (hours, buses) MVAr load
    rtpv_bus: np.ndarray    # (hours, buses) MW rooftop PV (fixed injection)
    avail: np.ndarray       # (hours, gens) MW available for PV / WIND / HYDRO units (NaN for others)
    pv_shape_area: dict     # area -> (hours,) utility PV output per MW of installed capacity


def load_timeseries(rts, load_scale=1.0):
    load = _read_ts("Load/DAY_AHEAD_regional_Load.csv")
    nh, nb = len(load), len(rts.bus)
    pd_bus = np.zeros((nh, nb))
    qd_bus = np.zeros((nh, nb))
    for area in (1, 2, 3):
        m = (rts.bus.area == area).values
        share = rts.bus.pd.values[m] / rts.bus.pd.values[m].sum()
        scale = load[str(area)].values[:, None] * share[None, :] * load_scale
        pd_bus[:, m] = scale
        with np.errstate(invalid="ignore", divide="ignore"):
            pf_ratio = np.where(rts.bus.pd.values[m] > 0, rts.bus.qd.values[m] / rts.bus.pd.values[m], 0)
        qd_bus[:, m] = scale * pf_ratio[None, :]

    avail = np.full((nh, len(rts.gen)), np.nan)
    series = pd.concat([_read_ts("PV/DAY_AHEAD_pv.csv"), _read_ts("WIND/DAY_AHEAD_wind.csv"),
                        _read_ts("Hydro/DAY_AHEAD_hydro.csv")], axis=1)
    for i, g in rts.gen.iterrows():
        if g.type in VARIABLE:
            avail[:, i] = np.minimum(series[g.uid].values, g.pmax)

    rtpv = _read_ts("RTPV/DAY_AHEAD_rtpv.csv")
    pos = {bid: i for i, bid in enumerate(rts.bus.id)}
    rtpv_bus = np.zeros((nh, nb))
    for uid in rtpv.columns:
        rtpv_bus[:, pos[int(uid.split("_")[0])]] += rtpv[uid].values

    pv_shape_area = {}
    for area in (1, 2, 3):
        idx = [i for i, g in rts.gen.iterrows() if g.type == "PV" and g.bus // 100 == area]
        cap = rts.gen.pmax.values[idx].sum()
        pv_shape_area[area] = np.nansum(avail[:, idx], axis=1) / cap
    return TimeSeries(index=load.index, pd_bus=pd_bus, qd_bus=qd_bus, rtpv_bus=rtpv_bus, avail=avail,
                      pv_shape_area=pv_shape_area)


# ----------------------------------------------------------------------------- AC network
@dataclass
class ACNet:
    net: object
    bus_pos: dict           # bus id -> pandapower bus index
    branch_elem: list       # per branch: ("line"|"trafo", index)
    ctrl_gen: dict          # bus id -> pandapower gen index (voltage-controlled plant at that bus)
    gens_at_bus: dict       # bus id -> list of conventional gen row indices
    var_sgen: dict          # gen row index -> sgen index (PV / WIND / HYDRO)
    rtpv_sgen: dict         # bus id -> sgen index
    slack_gens: list        # gen rows aggregated at the slack bus


def build_ac(rts):
    """pandapower model for hourly AC power flow with dispatch from the DC-OPF."""
    p = rts.ppc
    net = from_ppc({k: p[k] for k in ("version", "baseMVA", "bus", "gen", "branch")}, f_hz=60,
                   validate_conversion=False)
    lk = net._from_ppc_lookups
    # transformer taps on the MATPOWER from-bus side
    branch_elem = list(zip(lk["branch"].element_type, lk["branch"].element.astype(int)))
    for k, (et, ei) in enumerate(branch_elem):
        if et == "trafo":
            fbus = rts.branch.f[k]
            net.trafo.at[ei, "tap_side"] = "lv" if net.trafo.at[ei, "lv_bus"] == fbus else "hv"
    bus_pos = {bid: bid for bid in rts.bus.id}          # from_ppc keeps MATPOWER bus numbers as index
    # replace converted generators by: one controlled gen per plant bus + sgens for variable/RTPV
    net.gen.drop(net.gen.index, inplace=True)
    net.sgen.drop(net.sgen.index, inplace=True)
    net.ext_grid.drop(net.ext_grid.index, inplace=True)
    gens_at_bus, ctrl_gen, var_sgen = {}, {}, {}
    for i, g in rts.gen.iterrows():
        if g.type in VOLTAGE_CONTROL:
            gens_at_bus.setdefault(g.bus, []).append(i)
    for b, rows in gens_at_bus.items():
        vg = rts.gen.vg[rows[0]]
        if b == SLACK_BUS:
            pp.create_ext_grid(net, b, vm_pu=vg, va_degree=0.0, name=f"plant {b}")
            continue
        ctrl_gen[b] = pp.create_gen(net, b, p_mw=0.0, vm_pu=vg, name=f"plant {b}",
                                    min_q_mvar=rts.gen.qmin[rows].sum(), max_q_mvar=rts.gen.qmax[rows].sum())
    for i, g in rts.gen.iterrows():
        if g.type in VARIABLE and g.type not in VOLTAGE_CONTROL:          # PV and wind: P only
            var_sgen[i] = pp.create_sgen(net, g.bus, p_mw=0.0, q_mvar=0.0, name=g.uid)
    rtpv_sgen = {b: pp.create_sgen(net, b, p_mw=0.0, name=f"RTPV {b}") for b in rts.bus.id}
    return ACNet(net=net, bus_pos=bus_pos, branch_elem=branch_elem, ctrl_gen=ctrl_gen, gens_at_bus=gens_at_bus,
                 var_sgen=var_sgen, rtpv_sgen=rtpv_sgen, slack_gens=gens_at_bus.get(SLACK_BUS, []))


def branch_flows(ac):
    """Return arrays (p_from, q_from, p_to, q_to) in branch order from the last power flow."""
    n = ac.net
    out = np.zeros((4, len(ac.branch_elem)))
    for k, (et, ei) in enumerate(ac.branch_elem):
        if et == "line":
            r = n.res_line.loc[ei]
            out[:, k] = r.p_from_mw, r.q_from_mvar, r.p_to_mw, r.q_to_mvar
        else:
            r = n.res_trafo.loc[ei]
            hv_is_from = n.trafo.at[ei, "tap_side"] == "hv"       # tap side was set to the from bus
            a = (r.p_hv_mw, r.q_hv_mvar, r.p_lv_mw, r.q_lv_mvar)
            out[:, k] = a if hv_is_from else (a[2], a[3], a[0], a[1])
    return out
