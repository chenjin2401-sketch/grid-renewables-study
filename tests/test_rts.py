"""RTS-GMLC model checks (skipped when the data has not been downloaded)."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gridstudy.rts import CASE  # noqa: E402

pytestmark = pytest.mark.skipif(not CASE.exists(), reason="run scripts/download_rts_gmlc.py first")


@pytest.fixture(scope="module")
def rts():
    from gridstudy.rts import load_rts
    return load_rts()


def test_ac_matches_matpower_reference(rts):
    from gridstudy.rts import build_ac
    from gridstudy.validate_rts import base_case_ac, reference_ac_pf
    vm_ref, va_ref, pf_ref, loss_ref = reference_ac_pf()
    vm, va, pf, loss = base_case_ac(rts, build_ac)
    ids = list(rts.bus.id)
    assert np.abs(vm - [vm_ref[b] for b in ids]).max() < 6e-4        # reference has 3 decimals
    assert np.abs(va - [va_ref[b] for b in ids]).max() < 6e-4
    assert np.abs(pf - pf_ref).max() < 6e-3                           # 2 decimals
    assert abs(loss - loss_ref) < 6e-3


def test_dcopf_matches_matpower_reference(rts):
    from gridstudy.dcopf import DCOPF
    from gridstudy.validate_rts import reference_dc_opf
    obj, lam = reference_dc_opf()
    r = DCOPF(rts, mode="matpower").solve(rts.bus.pd.values)
    assert r.cost == pytest.approx(obj, abs=0.01)
    assert np.allclose(r.lmp, lam, atol=0.01)


def test_lmp_equals_finite_difference(rts):
    """LMP from the LP duals equals the cost change for 1 MW more load at a bus (congested hour)."""
    from gridstudy.dcopf import DCOPF
    opf = DCOPF(rts)
    rating = rts.branch.rate_a.values.copy()
    rating[rts.branch.uid.values == "A27"] = 150.0          # force congestion so LMPs differ
    d = rts.bus.pd.values.copy()
    base = opf.solve(d, rating=rating)
    assert base.lmp.max() - base.lmp.min() > 1.0
    for i in (0, 30, 60):
        d2 = d.copy()
        d2[i] += 1.0
        assert opf.solve(d2, rating=rating).cost - base.cost == pytest.approx(base.lmp[i], abs=0.05)


def test_ptdf_matches_dc_power_flow(rts):
    """PTDF flows satisfy Kirchhoff's current law at every non-slack bus."""
    rng = np.random.default_rng(0)
    inj = rng.normal(size=len(rts.bus))
    inj -= inj.mean()
    flow = rts.ptdf @ inj
    pos = {b: i for i, b in enumerate(rts.bus.id)}
    f = rts.branch.f.map(pos).values
    t = rts.branch.t.map(pos).values
    bal = np.zeros(len(rts.bus))
    np.add.at(bal, f, flow)
    np.add.at(bal, t, -flow)
    s = pos[113]
    mask = np.arange(len(rts.bus)) != s
    assert np.allclose(bal[mask], inj[mask], atol=1e-9)
