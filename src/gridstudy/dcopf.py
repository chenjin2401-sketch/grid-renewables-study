"""DC optimal power flow (economic dispatch with transmission limits) as an LP, PTDF form.

    min   sum_g sum_s slope_gs * x_gs  +  VOLL * sum_b shed_b
    s.t.  sum_g P_g + sum_b shed_b = sum_b d_b                       (power balance)
          -R_l <= PTDF_l . (Cg P + shed - d) <= R_l                   (branch limits)
          P_g = Pmin_g + sum_s x_gs,  0 <= x_gs <= width_gs           (piecewise-linear cost)

d_b is net load (load - rooftop PV). LMP_b = lambda + PTDF[:, b] . (mu_up - mu_dn).
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, hstack, vstack

from .rts import THERMAL, VARIABLE

VOLL = 10_000.0


@dataclass
class DCOPFResult:
    pg: np.ndarray          # MW per generator
    flow: np.ndarray        # MW per branch
    lmp: np.ndarray         # $/MWh per bus
    cost: float             # $/h (including fixed cost at Pmin in MATPOWER mode)
    shed: float             # MW of unserved load
    binding: np.ndarray     # bool per branch: at its limit
    status: int


class DCOPF:
    """Reusable LP structure; only bounds and right-hand sides change from hour to hour."""

    def __init__(self, rts, rating="rate_a", mode="relaxed", mva_factor=None):
        """mode='relaxed': thermal units 0..Pmax with convexified costs (time series).
        mode='matpower': committed units only, Pmin enforced, original PWL costs (validation).
        mva_factor: per-branch AC MVA / DC MW ratio; DC limits become rating / mva_factor so the
        dispatch leaves room for reactive flow (default: rts.mva_factor if set, else 1)."""
        self.rts, self.mode = rts, mode
        gen, nb = rts.gen, len(rts.bus)
        pos = {bid: i for i, bid in enumerate(rts.bus.id)}
        self.bus_of_gen = gen.bus.map(pos).values
        cols, slopes, widths, self.fixed_mw, self.fixed_cost = [], [], [], np.zeros(len(gen)), 0.0
        for i, g in gen.iterrows():
            if mode == "matpower":
                if g.status == 0:
                    continue
                gc = rts.ppc["gencost"][i]
                n = int(gc[3])
                x = gc[4:4 + 2 * n:2]
                y = gc[5:5 + 2 * n:2]
                self.fixed_mw[i] = x[0]
                self.fixed_cost += y[0]
                segs = [(x[k + 1] - x[k], (y[k + 1] - y[k]) / (x[k + 1] - x[k])) for k in range(n - 1)
                        if x[k + 1] > x[k]]
            elif g.type in THERMAL or g.type in VARIABLE:
                segs = rts.cost_segments[i] if g.type in THERMAL else [(g.pmax, 0.0)]
            else:
                continue
            room = g.pmax - self.fixed_mw[i]                # never exceed Pmax
            for w, s in segs:
                w = min(w, max(room, 0.0))
                room -= w
                if w <= 0:
                    continue
                cols.append(i)
                widths.append(w)
                slopes.append(s)
        self.seg_gen = np.array(cols)
        self.seg_width = np.array(widths, dtype=float)
        nseg = len(cols)
        # variables: [segments | shed per bus]
        self.c = np.r_[np.array(slopes), np.full(nb, VOLL)]
        inj_bus = np.r_[self.bus_of_gen[self.seg_gen], np.arange(nb)]           # bus of each variable
        ptdf = rts.ptdf
        a_flow = ptdf[:, inj_bus]                                               # (lines, vars)
        self.a_ub = csr_matrix(np.vstack([a_flow, -a_flow]))
        self.a_eq = csr_matrix(np.ones((1, nseg + nb)))
        if mva_factor is None:
            mva_factor = getattr(rts, "mva_factor", None)
        self.mva_factor = np.ones(len(rts.branch)) if mva_factor is None or mode == "matpower" \
            else np.asarray(mva_factor, dtype=float)
        self.rating = rts.branch[rating].values.astype(float) / self.mva_factor
        self.nseg, self.nb = nseg, nb
        self.var_seg_idx = {}                                                   # gen -> its segment columns
        for j, gi in enumerate(self.seg_gen):
            self.var_seg_idx.setdefault(gi, []).append(j)

    def max_injection(self, d_bus, avail, bus_idx, rating=None, return_binding=False):
        """Largest new zero-cost injection at bus_idx the network can take this hour, with every
        existing unit (thermal and renewable) redispatched within its limits and all branches
        within rating (OPF-based hosting capacity). Returns MW (0 if the base case is infeasible);
        with return_binding=True also the boolean mask of branches at their limit."""
        ptdf = self.rts.ptdf
        fixed_inj = np.zeros(self.nb)
        np.add.at(fixed_inj, self.bus_of_gen, self.fixed_mw)
        net_d = d_bus - fixed_inj
        r = self.rating if rating is None else np.asarray(rating, float) / self.mva_factor
        base_flow = ptdf @ (-net_d)
        a_ub = hstack([self.a_ub, csr_matrix(np.r_[ptdf[:, bus_idx], -ptdf[:, bus_idx]][:, None])]).tocsr()
        a_eq = hstack([self.a_eq, csr_matrix(np.ones((1, 1)))]).tocsr()
        ub = np.r_[self._seg_ub(avail), np.zeros(self.nb), 1e5]            # no load shedding
        c = np.zeros(self.nseg + self.nb + 1)
        c[-1] = -1.0
        res = linprog(c, A_ub=a_ub, b_ub=np.r_[r - base_flow, r + base_flow], A_eq=a_eq,
                      b_eq=np.array([net_d.sum()]), bounds=np.c_[np.zeros_like(ub), ub], method="highs")
        mw = res.x[-1] if res.status == 0 else 0.0
        if not return_binding:
            return mw
        if res.status != 0:
            return mw, np.zeros(len(r), bool)
        nl = len(r)
        binding = (np.abs(res.ineqlin.marginals[:nl]) > 1e-7) | (np.abs(res.ineqlin.marginals[nl:]) > 1e-7)
        return mw, binding

    def _seg_ub(self, avail):
        ub = self.seg_width.copy()
        if avail is not None:
            for gi, cols in self.var_seg_idx.items():
                if not np.isnan(avail[gi]):
                    cap = avail[gi]
                    for j in cols:
                        ub[j] = min(self.seg_width[j], max(cap, 0))
                        cap -= ub[j]
        return ub

    def solve(self, d_bus, avail=None, extra_inj=None, rating=None):
        """d_bus: (buses,) net load MW. avail: (gens,) MW cap for variable units (NaN = not variable).
        extra_inj: (buses,) fixed MW injections (e.g. a new solar plant). rating: override (branches,)."""
        rts = self.rts
        ptdf = rts.ptdf
        fixed_inj = np.zeros(self.nb)
        np.add.at(fixed_inj, self.bus_of_gen, self.fixed_mw)
        if extra_inj is not None:
            fixed_inj = fixed_inj + extra_inj
        net_d = d_bus - fixed_inj
        r = self.rating if rating is None else np.asarray(rating, float) / self.mva_factor
        base_flow = ptdf @ (-net_d)
        b_ub = np.r_[r - base_flow, r + base_flow]
        b_eq = np.array([net_d.sum()])
        ub = np.r_[self._seg_ub(avail), np.maximum(d_bus, 0)]
        bounds = np.c_[np.zeros_like(ub), ub]
        res = linprog(self.c, A_ub=self.a_ub, b_ub=b_ub, A_eq=self.a_eq, b_eq=b_eq, bounds=bounds,
                      method="highs")
        if res.status != 0:
            return DCOPFResult(np.full(len(rts.gen), np.nan), np.full(len(r), np.nan),
                               np.full(self.nb, np.nan), np.nan, np.nan, np.zeros(len(r), bool), res.status)
        x = res.x
        pg = self.fixed_mw.copy()
        np.add.at(pg, self.seg_gen, x[:self.nseg])
        shed = x[self.nseg:]
        inj = -net_d + shed
        np.add.at(inj, self.bus_of_gen, pg - self.fixed_mw)
        flow = ptdf @ inj
        mu = res.ineqlin.marginals
        nl = len(r)
        lmp = res.eqlin.marginals[0] + ptdf.T @ (mu[:nl] - mu[nl:])
        binding = np.abs(np.abs(flow) - r) < 1e-4
        return DCOPFResult(pg=pg, flow=flow, lmp=lmp, cost=res.fun + self.fixed_cost, shed=shed.sum(),
                           binding=binding, status=0)
