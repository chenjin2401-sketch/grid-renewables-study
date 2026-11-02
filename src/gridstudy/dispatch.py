"""Economic dispatch."""
import numpy as np
import pyomo.environ as pyo


def economic_dispatch_lp(demand_mw, price, p_min, p_max, solver="appsi_highs"):
    """Lossless single-period economic dispatch as a linear program.

        min  sum_i price_i * P_i
        s.t. sum_i P_i = demand,   p_min_i <= P_i <= p_max_i

    Returns (P, lambda): the optimal outputs and the system incremental cost
    (dual of the power-balance constraint, $/MWh).
    """
    n = len(price)
    m = pyo.ConcreteModel()
    m.I = pyo.RangeSet(0, n - 1)
    m.P = pyo.Var(m.I, bounds=lambda m, i: (p_min[i], p_max[i]))
    m.cost = pyo.Objective(expr=sum(price[i] * m.P[i] for i in m.I), sense=pyo.minimize)
    m.balance = pyo.Constraint(expr=sum(m.P[i] for i in m.I) == demand_mw)
    m.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)
    res = pyo.SolverFactory(solver).solve(m)
    if res.solver.termination_condition != pyo.TerminationCondition.optimal:
        raise RuntimeError(f"dispatch not optimal: {res.solver.termination_condition}")
    p = np.array([pyo.value(m.P[i]) for i in m.I])
    return p, abs(m.dual[m.balance])


def merit_order(demand_mw, price, p_min, p_max):
    """Closed-form merit-order dispatch for linear costs (used to check the LP)."""
    p = np.array(p_min, dtype=float)
    rest = demand_mw - p.sum()
    for i in np.argsort(price):
        add = min(p_max[i] - p[i], rest)
        p[i] += add
        rest -= add
    if rest > 1e-9:
        raise ValueError("demand exceeds total capacity")
    return p
