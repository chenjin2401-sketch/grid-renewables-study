"""Read a MATPOWER .m case file (numeric matrices and name cell arrays) without MATLAB."""
import re

import numpy as np


def _matrix(text, name):
    m = re.search(r"mpc\." + name + r"\s*=\s*\[(.*?)\];", text, re.S)
    if not m:
        return None
    rows = []
    for line in m.group(1).splitlines():
        line = line.split("%")[0].strip().rstrip(";").strip()
        if line:
            rows.append([float(v) for v in line.replace(";", " ").split()])
    return np.array(rows)


def _cells(text, name):
    m = re.search(r"mpc\." + name + r"\s*=\s*\{(.*?)\};", text, re.S)
    if not m:
        return None
    out = []
    for line in m.group(1).splitlines():
        vals = re.findall(r"'([^']*)'", line)
        if vals:
            out.append(vals if len(vals) > 1 else vals[0])
    return out


def read_case(path):
    """Return a ppc-style dict: baseMVA, bus, gen, branch, gencost, dcline, bus_name, gen_name."""
    text = open(path, encoding="utf-8").read()
    base = float(re.search(r"mpc\.baseMVA\s*=\s*([\d.]+)", text).group(1))
    ppc = {"version": "2", "baseMVA": base}
    for name in ("bus", "gen", "branch", "gencost", "dcline", "areas"):
        mat = _matrix(text, name)
        if mat is not None:
            ppc[name] = mat
    for name in ("bus_name", "gen_name"):
        cells = _cells(text, name)
        if cells is not None:
            ppc[name] = cells
    return ppc
