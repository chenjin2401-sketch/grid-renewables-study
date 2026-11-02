"""Test networks built in pandapower.

five_bus() is the Phase 0 validation case: MATPOWER case5 (PJM 5-bus system,
F. Li and R. Bo, 2010) with the modifications used in the author's ECE 451
coursework: bus 4 is the slack, all other buses are PQ, every branch is rated
400 MVA, and the generators have linear costs.
"""
from dataclasses import dataclass, field

import numpy as np
import pandapower as pp

BASE_MVA = 100.0
F_HZ = 60.0


@dataclass
class GenInfo:
    """Generator data kept next to the pandapower net (index = generator number - 1)."""
    name: list
    bus: list            # 1-based bus numbers, as in MATPOWER
    p_init_mw: list      # specified (initial) dispatch
    p_max_mw: list
    q_max_mvar: list
    price: list          # $/MWh
    slack: int           # 0-based index of the slack generator
    element: list = field(default_factory=list)   # ("sgen", idx) or ("ext_grid", idx)


# MATPOWER-style tables -------------------------------------------------------
#        bus  type  Pd   Qd     baseKV
_BUS = [(1, 1, 0, 0, 230),
        (2, 1, 300, 98.61, 230),
        (3, 1, 300, 98.61, 230),
        (4, 3, 400, 131.47, 230),
        (5, 1, 0, 0, 230)]
#        fbus tbus  r        x       b (total, pu)  rateA (MVA)
_BRANCH = [(1, 2, 0.00281, 0.0281, 0.00712, 400),
           (1, 4, 0.00304, 0.0304, 0.00658, 400),
           (1, 5, 0.00064, 0.0064, 0.03126, 400),
           (2, 3, 0.00108, 0.0108, 0.01852, 400),
           (3, 4, 0.00297, 0.0297, 0.00674, 400),
           (4, 5, 0.00297, 0.0297, 0.00674, 400)]
#        name, bus, Pg0, Pmax, Qmax, $/MWh
_GEN = [("G1", 1, 20, 100, 30, 10),
        ("G2", 1, 100, 500, 127.5, 20),
        ("G3", 3, 50, 100, 390, 30),
        ("G4", 4, 230, 1000, 150, 5),
        ("G5", 5, 200, 300, 450, 80)]


def five_bus():
    """Return (net, gens) for the modified 5-bus system.

    Branch data are converted from per unit on 100 MVA / 230 kV to ohms and nF
    for a 1 km line, so the pandapower pi model equals the MATPOWER branch model:
        Z_base = kV^2 / MVA,  R = r * Z_base,  X = x * Z_base,
        B_total = b / Z_base  ->  C = B_total / (2*pi*f).
    """
    net = pp.create_empty_network(sn_mva=BASE_MVA, f_hz=F_HZ)
    bus_idx = {}
    for (b, btype, pd, qd, kv) in _BUS:
        bus_idx[b] = pp.create_bus(net, vn_kv=kv, name=f"Bus {b}", index=b - 1)
        pp.create_load(net, bus_idx[b], p_mw=pd, q_mvar=qd, name=f"Load {b}")

    for (f, t, r, x, b, rate) in _BRANCH:
        kv = _BUS[f - 1][4]
        z_base = kv ** 2 / BASE_MVA
        c_nf = b / z_base / (2 * np.pi * F_HZ) * 1e9
        pp.create_line_from_parameters(
            net, bus_idx[f], bus_idx[t], length_km=1.0,
            r_ohm_per_km=r * z_base, x_ohm_per_km=x * z_base, c_nf_per_km=c_nf,
            max_i_ka=rate / (np.sqrt(3) * kv), name=f"Line {f}-{t}")

    slack_bus = next(b for (b, btype, *_rest) in _BUS if btype == 3)
    gens = GenInfo(name=[g[0] for g in _GEN], bus=[g[1] for g in _GEN],
                   p_init_mw=[g[2] for g in _GEN], p_max_mw=[g[3] for g in _GEN],
                   q_max_mvar=[g[4] for g in _GEN], price=[g[5] for g in _GEN],
                   slack=next(i for i, g in enumerate(_GEN) if g[1] == slack_bus))
    for i, (name, b, pg, pmax, qmax, price) in enumerate(_GEN):
        if i == gens.slack:
            idx = pp.create_ext_grid(net, bus_idx[b], vm_pu=1.0, va_degree=0.0, name=name)
            gens.element.append(("ext_grid", idx))
        else:
            # generator at a PQ bus: fixed P and Q injection (pandapower "static generator")
            idx = pp.create_sgen(net, bus_idx[b], p_mw=pg, q_mvar=0.0, name=name)
            gens.element.append(("sgen", idx))
    return net, gens


def branch_labels(net):
    return [f"{net.bus.at[f, 'name'].split()[-1]}-{net.bus.at[t, 'name'].split()[-1]}"
            for f, t in zip(net.line.from_bus, net.line.to_bus)]
