"""Mitigations for a solar project above the firm hosting capacity.

- flexible interconnection: accept curtailment when the grid has no room (hourly headroom)
- co-located battery: store solar that would be curtailed and export it later (yearly LP)
- line upgrade: raise the rating of the limiting branch and recompute the headroom
- economics: annualized cost per MWh of curtailment avoided
"""
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import diags, eye, hstack, vstack, csr_matrix


def curtailment(size_mw, shape, head):
    """Fraction of potential solar energy curtailed when export is capped by the hourly headroom."""
    pot = size_mw * shape
    deliv = np.minimum(pot, np.maximum(head, 0))
    return 1 - deliv.sum() / pot.sum(), pot.sum(), deliv.sum()


def size_for_curtailment(shape, head, cap_pct, hi=3000.0, tol=1.0):
    """Largest plant size (MW) whose annual curtailment stays at or below cap_pct."""
    lo = 0.0
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if curtailment(mid, shape, head)[0] * 100 <= cap_pct:
            lo = mid
        else:
            hi = mid
    return lo


def battery_dispatch(potential, head, price, p_mw, e_mwh, rte=0.85, capture_only=False):
    """Yearly LP for a battery charged only from the co-located solar plant.

    max  sum_t (price_t + 1e-3) * x_t          x_t = exported MW (solar used + discharge - charge)
    s.t. x_t = s_t - c_t + d_t,  0 <= s_t <= potential_t,  c_t <= s_t
         0 <= x_t <= head_t,  0 <= c_t, d_t <= p_mw
         e_{t+1} = e_t + eta c_t - d_t / eta,  0 <= e_t <= e_mwh,  e_0 = e_T = 0
    capture_only: charge only solar that the grid cannot take that hour
        (c_t <= max(potential_t - head_t, 0)) -> isolates curtailment recovery from arbitrage.
    The small 1e-3 $/MWh bonus breaks ties toward delivering more energy.
    Returns exported MW per hour.
    """
    n = len(potential)
    head = np.maximum(head, 0)
    if p_mw <= 0 or e_mwh <= 0:
        return np.minimum(potential, head)
    eta = np.sqrt(rte)
    # variables: s (n), c (n), d (n), e (n+1)
    I = eye(n, format="csr")
    Z = csr_matrix((n, n))
    zcol = csr_matrix((n, n + 1))
    x_row = hstack([I, -I, I, zcol])                         # export = s - c + d
    # energy balance: e_{t+1} - e_t - eta c_t + d_t / eta = 0
    shift = hstack([csr_matrix((n, 1)), eye(n)]) - hstack([eye(n), csr_matrix((n, 1))])
    a_eq = hstack([Z, -eta * I, I / eta, shift]).tocsr()
    b_eq = np.zeros(n)
    a_ub = vstack([x_row, -x_row, hstack([-I, I, Z, zcol])]).tocsr()   # x <= head, x >= 0, c <= s
    b_ub = np.r_[head, np.zeros(n), np.zeros(n)]
    w = price + 1e-3
    c = -np.r_[w, -w, w, np.zeros(n + 1)]
    c_cap = np.minimum(p_mw, np.maximum(potential - head, 0)) if capture_only else np.full(n, p_mw)
    bounds = ([(0, p) for p in potential] + [(0, c) for c in c_cap] + [(0, p_mw)] * n
              + [(0, 0)] + [(0, e_mwh)] * (n - 1) + [(0, 0)])
    res = linprog(c, A_ub=a_ub, b_ub=b_ub, A_eq=a_eq, b_eq=b_eq, bounds=bounds, method="highs")
    if res.status != 0:
        raise RuntimeError(res.message)
    s, ch, d = res.x[:n], res.x[n:2 * n], res.x[2 * n:3 * n]
    return s - ch + d


def crf(rate, years):
    """Capital recovery factor."""
    return rate * (1 + rate) ** years / ((1 + rate) ** years - 1)


def annualized(capex, rate=0.07, years=20, fom_frac=0.0):
    return capex * crf(rate, years) + capex * fom_frac
