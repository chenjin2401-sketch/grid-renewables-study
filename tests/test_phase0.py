"""Phase 0 checks: dispatch logic always; MATPOWER validation when data/private/ exists."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gridstudy.dispatch import economic_dispatch_lp, merit_order   # noqa: E402
from gridstudy.networks import five_bus                            # noqa: E402

PRIVATE = ROOT / "data" / "private"
HAVE_REF = (PRIVATE / "load_profile_5bus.csv").exists() and (PRIVATE / "matpower_reference").exists()

PRICE = np.array([10, 20, 30, 5, 80.0])
PMAX = np.array([100, 500, 100, 300, 300.0])


@pytest.mark.parametrize("demand, lam", [(640, 20), (950, 30), (1237, 80)])
def test_lp_matches_merit_order(demand, lam):
    p, dual = economic_dispatch_lp(demand, PRICE, np.zeros(5), PMAX)
    assert np.allclose(p, merit_order(demand, PRICE, np.zeros(5), PMAX), atol=1e-6)
    assert p.sum() == pytest.approx(demand)
    assert dual == pytest.approx(lam)


def test_five_bus_base_case_converges():
    import pandapower as pp
    net, gens = five_bus()
    pp.runpp(net, calculate_voltage_angles=True, numba=False)
    assert net.converged
    # power balance: generation = load + losses
    gen = net.res_sgen.p_mw.sum() + net.res_ext_grid.p_mw.sum()
    assert gen == pytest.approx(net.res_load.p_mw.sum() + net.res_line.pl_mw.sum(), abs=1e-6)


@pytest.fixture(scope="module")
def phase0():
    if not HAVE_REF:
        pytest.skip("data/private/ not present (coursework reference data)")
    from gridstudy.scenarios import economic_dispatch, fixed_schedule, load_profile
    from gridstudy.validate import load_reference
    p, q = load_profile(PRIVATE / "load_profile_5bus.csv")
    ref_dir = PRIVATE / "matpower_reference"
    return {
        "A": (fixed_schedule(p, q)[0], load_reference(ref_dir, "task1")),
        "B": (economic_dispatch(p, q)[0], load_reference(ref_dir, "task2")),
    }


TOL = {"vm": 1e-5, "va": 1e-3, "p_from": 1e-3, "p_to": 1e-3, "pg": 1e-3, "qg": 1e-3,
       "loss": 1e-3, "cost": 0.05}


@pytest.mark.parametrize("scenario", ["A", "B"])
@pytest.mark.parametrize("key", list(TOL))
def test_matches_matpower(phase0, scenario, key):
    res, ref = phase0[scenario]
    assert res.converged.all()
    diff = np.abs(np.asarray(getattr(res, key)) - ref[key]).max()
    assert diff < TOL[key], f"{scenario}/{key}: max diff {diff:.3e}"
