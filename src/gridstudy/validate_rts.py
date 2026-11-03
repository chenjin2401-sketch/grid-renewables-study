"""Validate the RTS-GMLC model against the reference MATPOWER output shipped with RTS-GMLC."""
import numpy as np
import pandapower as pp
from pandapower.converter.pypower import from_ppc

from .rts import CASE, DATA, branch_flows

REF = DATA / "FormattedData" / "MATPOWER" / "MATPOWER-out.txt"


def _section(lines, title):
    start = next(i for i, l in enumerate(lines) if "MATPOWER Version" in l and title in l)
    end = next((i for i in range(start + 1, len(lines)) if "MATPOWER Version" in lines[i]), len(lines))
    return lines[start:end]


def _table(sec, header_start, n):
    i0 = next(i for i, l in enumerate(sec) if l.startswith(header_start))
    return [l.split() for l in sec[i0 + 3:i0 + 3 + n]]


def reference_ac_pf(nb=73, nl=120):
    lines = REF.read_text().splitlines()
    sec = _section(lines, "AC Power Flow")
    bus = _table(sec, " Bus      Voltage", nb)
    br = _table(sec, "Brnch   From", nl)
    vm = {int(r[0]): float(r[1]) for r in bus}
    va = {int(r[0]): float(r[2].rstrip("*")) for r in bus}
    p_from = np.array([float(r[3]) for r in br])
    loss = float(next(l for l in sec if "Losses (I^2 * Z)" in l).split("Z)")[1].split()[0])
    return vm, va, p_from, loss


def reference_dc_opf():
    lines = REF.read_text().splitlines()
    sec = _section(lines, "DC Optimal Power Flow")
    obj = float(next(l for l in sec if "Objective Function Value" in l).split("=")[1].split()[0])
    lam = float(next(l for l in sec if l.startswith("Lambda P")).split()[2])
    return obj, lam


def base_case_ac(rts, ac_builder):
    """AC power flow of the MATPOWER base case (original dispatch) with the tap fix."""
    ac = ac_builder(rts)
    p = rts.ppc
    net = from_ppc({k: p[k] for k in ("version", "baseMVA", "bus", "gen", "branch")}, f_hz=60,
                   validate_conversion=False)
    for k, (et, ei) in enumerate(ac.branch_elem):
        if et == "trafo":
            net.trafo.at[ei, "tap_side"] = "lv" if net.trafo.at[ei, "lv_bus"] == rts.branch.f[k] else "hv"
    pp.runpp(net, calculate_voltage_angles=True, init="flat", numba=False)
    ac.net = net
    pf, _, _, _ = branch_flows(ac)
    vm = net.res_bus.loc[list(rts.bus.id), "vm_pu"].values
    va = net.res_bus.loc[list(rts.bus.id), "va_degree"].values
    loss = net.res_line.pl_mw.sum() + net.res_trafo.pl_mw.sum()
    return vm, va, pf, loss


__all__ = ["reference_ac_pf", "reference_dc_opf", "base_case_ac", "CASE"]
