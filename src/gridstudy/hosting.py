"""Solar hosting capacity: how many MW of new solar a bus can take without network violations.

Definition (OPF-based hosting capacity):
- In each hour, the network headroom at bus k is the largest new zero-cost injection at k
  for which a DC optimal power flow can still keep every branch within its normal rating,
  redispatching existing units (thermal down, renewables curtailed only if needed).
  This separates network limits from the market: a bus is limited only when no redispatch
  can make room.
- New plant output = capacity x the hourly utility-PV shape of its area, unity power factor.
- Firm hosting capacity = largest capacity whose output never exceeds the headroom
  (no curtailment in any hour), evaluated on a set of critical hours per bus.
- AC verification: at the firm capacity, the redispatched AC power flow must not create a new
  thermal overload attributable to the plant (branch distribution factor >= 5 %, loading
  increase >= 5 % of the plant MW; an existing overload must also grow by >= 1 % of its rating)
  or, in the study area (buses within STUDY_HOPS branches of the POI), a voltage outside
  0.95-1.05 pu (or more than 0.005 pu beyond the reference value where that already violates).
  The reference is the same redispatched fleet with the plant's MW taken by the system slack
  (superposition), so only the plant's own injection is judged, not the market's response.
"""
import numpy as np

from .ac_hourly import V_LIMITS, set_hour, solve_hour
from .rts import THERMAL, VARIABLE

V_TOL = 0.005
MIN_SHAPE = 0.05         # ignore hours with almost no solar when converting MW to capacity
DFAX_CUTOFF = 0.05       # attribution threshold for AC thermal overloads
WORSEN_FLOOR = 0.01      # an already-overloaded branch must load up by >= 1 % of its rating to count
AC_DC_ITER = 3           # AC-to-DC limit corrections per check (MW vs. MVA on DC-binding branches)
INV_OVERSIZE = 1.1       # inverter kVA = 1.1 x plant MW, so 0.44 pu reactive power is available at full output
STUDY_HOPS = 2           # voltage criteria are monitored within this many branches of the POI


def study_area(rts, bus_id, hops=STUDY_HOPS):
    """Boolean mask (buses) of the buses within `hops` branches of bus_id."""
    adj = {}
    for f, t in zip(rts.branch.f, rts.branch.t):
        adj.setdefault(f, set()).add(t)
        adj.setdefault(t, set()).add(f)
    seen, frontier = {bus_id}, {bus_id}
    for _ in range(hops):
        frontier = {n for b in frontier for n in adj.get(b, ())} - seen
        seen |= frontier
    return rts.bus.id.isin(seen).values


def candidate_hours(ts, y, area, n_sunny=300, n_low=300, n_rand=200, seed=0, dc_limit=None):
    """Critical sunny hours for the DC screen: sunniest, lowest net load, most congested, random."""
    shape = ts.pv_shape_area[area]
    sunny = np.where(shape >= MIN_SHAPE)[0]
    net = (ts.pd_bus.sum(1) - ts.rtpv_bus.sum(1))[sunny]
    rng = np.random.default_rng(seed + area)
    load = np.abs(y["flow"][sunny]) / (dc_limit if dc_limit is not None else 1.0)     # % of DC limit
    congested = sunny[np.argsort(-load.max(1))[:n_rand]]
    picks = np.r_[sunny[np.argsort(-shape[sunny])[:n_sunny]], sunny[np.argsort(net)[:n_low]],
                  congested, rng.choice(sunny, size=min(n_rand, len(sunny)), replace=False)]
    return np.unique(picks)


def opf_headroom(opf, ts, bus_id, hours, rating=None):
    """Network headroom (MW of new injection) at bus_id for each hour in `hours`."""
    pos = {b: i for i, b in enumerate(opf.rts.bus.id)}
    return np.array([opf.max_injection(ts.pd_bus[h] - ts.rtpv_bus[h], ts.avail[h], pos[bus_id], rating=rating)
                     for h in hours])


def firm_from_headroom(head, shape_h):
    ok = shape_h >= MIN_SHAPE
    return float(np.min(head[ok] / shape_h[ok])) if ok.any() else np.inf


def redispatch_with_plant(opf, ts, h, bus_id, p_new, rating=None):
    """DC-OPF dispatch of the existing fleet with the new plant injecting p_new MW at bus_id."""
    pos = {b: i for i, b in enumerate(opf.rts.bus.id)}
    extra = np.zeros(len(opf.rts.bus))
    extra[pos[bus_id]] = p_new
    return opf.solve(ts.pd_bus[h] - ts.rtpv_bus[h], avail=ts.avail[h], extra_inj=extra, rating=rating)


def static_dfax(rts, bus_id):
    """Distribution factors of an injection at bus_id sunk by the thermal fleet pro-rata to Pmax.

    Used, as in PJM/MISO studies, to decide which branches a plant can be held responsible for
    (|DFAX| >= 5 %); hour-to-hour OPF redispatch elsewhere in the system is not attributed to it.
    """
    pos = {b: i for i, b in enumerate(rts.bus.id)}
    thermal = rts.gen.type.isin(THERMAL).values
    w = np.where(thermal, rts.gen.pmax.values, 0.0)
    w = w / w.sum()
    gpos = rts.gen.bus.map(pos).values
    return rts.ptdf[:, pos[bus_id]] - rts.ptdf[:, gpos] @ w


class ACChecker:
    """AC criteria for a plant at one bus (attribution sets computed once)."""

    def __init__(self, rts, opf, ac, ts, y, bus_id, sgen_idx, q_of_v=None, s_rated=None, rate=None):
        self.rts, self.opf, self.ac, self.ts, self.y = rts, opf, ac, ts, y
        self.bus, self.sgen, self.q_of_v, self.s_rated = bus_id, sgen_idx, q_of_v, s_rated
        self.custom = rate is not None                            # e.g. an upgraded branch rating
        self.rate = rts.branch.rate_a.values if rate is None else np.asarray(rate, float)
        self.bi = {b: i for i, b in enumerate(rts.bus.id)}[bus_id]
        self.attributable = np.abs(static_dfax(rts, bus_id)) >= DFAX_CUTOFF
        self.local = study_area(rts, bus_id)

    def _ac_with_plant(self, h, r, p_new, s_rated):
        """AC power flow with the plant (and its Volt-VAR control, if any). Returns ACHour."""
        rts, ac, ts = self.rts, self.ac, self.ts
        set_hour(ac, rts, ts.pd_bus[h], ts.qd_bus[h], ts.rtpv_bus[h], r.pg, extra_sgen={self.sgen: (p_new, 0.0)})
        a = solve_hour(ac, rts)
        if self.q_of_v is not None and a.converged:              # damped fixed-point Volt-VAR iteration
            unity, q, settled = a, 0.0, None
            for _ in range(15):
                q = 0.5 * q + 0.5 * self.q_of_v(a.vm[self.bi], p_new, s_rated)
                ac.net.sgen.at[self.sgen, "q_mvar"] = q
                a = solve_hour(ac, rts)
                if not a.converged:
                    break
                if abs(self.q_of_v(a.vm[self.bi], p_new, s_rated) - q) < 0.5:
                    settled = a
                    break
            # an inverter that cannot settle on the curve falls back to unity power factor
            a = settled if settled is not None else unity
            ac.net.sgen.at[self.sgen, "q_mvar"] = 0.0
        return a

    def check(self, h, p_new, s_rated=None):
        """(ok, reason) for the plant injecting p_new MW in hour h.

        AC feedback to the dispatch: if an attributable branch is over its MVA rating in AC while
        the DC-OPF holds it exactly at its DC (MW) limit, that DC limit is tightened by the observed
        MVA/MW ratio and the fleet is redispatched (up to AC_DC_ITER times), as an operator's
        security-constrained dispatch would do. The plant fails only if redispatch cannot clear it.
        """
        rts, ac, ts, y, rate = self.rts, self.ac, self.ts, self.y, self.rate
        if not y["conv"][h]:
            return True, ""                                       # baseline AC failed this hour
        s_rated = s_rated or self.s_rated
        dc_rate = rate.copy()
        for it in range(AC_DC_ITER + 1):
            r = redispatch_with_plant(self.opf, ts, h, self.bus, p_new,
                                      rating=dc_rate if (self.custom or it > 0) else None)
            if r.status != 0:
                return False, "DC infeasible" if it == 0 else f"thermal {why_t} (after redispatch)"
            a = self._ac_with_plant(h, r, p_new, s_rated)
            if not a.converged:
                return False, "AC not converged"
            # Reference: the same redispatched fleet with the plant's MW taken by the system slack.
            # Comparing against it isolates the plant's own effect from the market redispatch.
            s_plant, v_plant = a.s_max, a.vm
            set_hour(ac, rts, ts.pd_bus[h], ts.qd_bus[h], ts.rtpv_bus[h], r.pg, extra_sgen={self.sgen: (0.0, 0.0)})
            b = solve_hour(ac, rts)
            base_s, base_v = (b.s_max, b.vm) if b.converged else (y["s_ac"][h], y["vm"][h])
            inc = s_plant - base_s
            # new overload: branch was within rating and the plant pushes it over (>= 5 % of its MW)
            new_over = (base_s <= rate) & (s_plant > rate + 1e-6) & (inc >= DFAX_CUTOFF * p_new)
            # existing overload: only a material worsening counts (>= 5 % of plant MW and >= 1 % of rating)
            worse = (base_s > rate) & (inc >= np.maximum(DFAX_CUTOFF * p_new, WORSEN_FLOOR * rate))
            over = (new_over | worse) & self.attributable
            if not over.any():
                break
            why_t = rts.branch.uid[np.argmax(np.where(over, s_plant / rate, 0))]
            fixable = over & r.binding                           # at the DC limit: MW/MVA gap, redispatch
            if not fixable.any() or it == AC_DC_ITER:
                return False, f"thermal {why_t}"
            dc_rate[fixable] *= rate[fixable] / s_plant[fixable] * 0.995
        hi = np.maximum(V_LIMITS[1], base_v + V_TOL)
        lo = np.minimum(V_LIMITS[0], base_v - V_TOL)
        bad = ((v_plant > hi) | (v_plant < lo)) & self.local
        if bad.any():
            return False, f"voltage {rts.bus.id[np.argmax(bad)]}"
        return True, ""

    def hourly_headroom(self, h, dc_head, need=None, tol=5.0, try_dc_first=True):
        """Largest injection in hour h that passes the AC criteria (bisection below dc_head).

        With `need`, stop as soon as `need` is known to pass: the result is then `need` with
        reason "AC ok", a lower bound on the headroom that is exact for any plant whose output
        that hour is <= need. Callers must pass the output of the LARGEST plant size they will
        evaluate. Bisection runs only when `need` itself fails.
        """
        if try_dc_first:
            ok, why = self.check(h, dc_head)
            if ok:
                return dc_head, "DC"
        else:
            why = ""
        hi = dc_head
        if need is not None and need < dc_head:
            ok, reason = self.check(h, need)
            if ok:
                return need, "AC ok"
            hi, why = need, reason
        lo = 0.0
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            ok, reason = self.check(h, mid)
            if ok:
                lo = mid
            else:
                hi, why = mid, reason
        return lo, why


def ac_check(rts, opf, ac, ts, y, hours, bus_id, cap_mw, sgen_idx, q_of_v=None):
    """AC power flow with the plant at cap_mw on the given hours. Returns (ok, reason, hour)."""
    chk = ACChecker(rts, opf, ac, ts, y, bus_id, sgen_idx, q_of_v, s_rated=INV_OVERSIZE * cap_mw)
    shape = ts.pv_shape_area[bus_id // 100]
    for h in hours:
        ok, why = chk.check(h, cap_mw * shape[h])
        if not ok:
            return False, why, h
    return True, "", None


def ac_critical_hours(hours, head, shape_h, n=24):
    """The tightest DC hours for this bus plus the sunniest ones."""
    ok = shape_h >= MIN_SHAPE
    hh, rr, ss = hours[ok], (head / np.where(shape_h > 0, shape_h, np.nan))[ok], shape_h[ok]
    return np.unique(np.r_[hh[np.argsort(rr)[: n // 2]], hh[np.argsort(-ss)[: n // 2]]])


def ac_hosting_capacity(rts, opf, ac, ts, y, bus_id, dc_cap, ac_hours, sgen_idx, q_of_v=None, tol=5.0):
    """Bisection on capacity with AC checks, starting from the DC result."""
    hi = dc_cap
    ok, reason, h = ac_check(rts, opf, ac, ts, y, ac_hours, bus_id, hi, sgen_idx, q_of_v)
    if ok:
        return hi, "thermal (DC, N-0)", None
    lo, why, where = 0.0, reason, h
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        ok, reason, h = ac_check(rts, opf, ac, ts, y, ac_hours, bus_id, mid, sgen_idx, q_of_v)
        if ok:
            lo = mid
        else:
            hi, why, where = mid, reason, h
    return lo, why, where


# ---------------------------------------------------------------- planning "sink" dispatch (N-1)
def displacement_weights(rts, pg_t):
    """(hours, gens) share of displaced MW per unit: thermal pro-rata to output; in hours with no
    thermal output, the dispatched PV / wind / hydro pro-rata. Also returns thermal MW per hour."""
    thermal = rts.gen.type.isin(THERMAL).values
    variable = rts.gen.type.isin(VARIABLE).values
    p = np.clip(np.nan_to_num(pg_t), 0, None)
    w = np.where(thermal[None, :], p, 0.0)
    tot = w.sum(1, keepdims=True)
    w = np.where(tot > 1.0, w, np.where(variable[None, :], p, 0.0))
    s = w.sum(1, keepdims=True)
    return np.divide(w, s, out=np.zeros_like(w), where=s > 0), tot[:, 0]


def injection_sensitivity(rts, bus_id, weights):
    """(hours, branches) MW change on each branch per MW injected at bus_id (balanced by weights)."""
    pos = {b: i for i, b in enumerate(rts.bus.id)}
    gpos = rts.gen.bus.map(pos).values
    return rts.ptdf[:, pos[bus_id]][None, :] - weights @ rts.ptdf[:, gpos].T


def volt_var_ieee1547(v, p_mw, s_rated):
    """IEEE 1547-2018 Category B default Volt-VAR curve (Q in MVAr, + = injection).

    V1..V4 = 0.92, 0.98, 1.02, 1.08 pu; Q = +/-44 % of rated apparent power, limited so that
    P^2 + Q^2 <= S_rated^2.
    """
    qmax = 0.44 * s_rated
    q = float(np.interp(v, [0.0, 0.92, 0.98, 1.02, 1.08, 2.0], [qmax, qmax, 0.0, 0.0, -qmax, -qmax]))
    q_lim = np.sqrt(max(s_rated ** 2 - p_mw ** 2, 0.0))
    return float(np.clip(q, -q_lim, q_lim))
