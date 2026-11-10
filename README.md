# Grid Renewables Study

An open, reproducible **interconnection study for solar + battery storage** on the 73-bus RTS-GMLC
transmission test system with a full year of hourly data: how much solar each bus can take, what fails
under N-1 contingencies, and what is the cheapest way to connect a project that is larger than the grid
can take today. Every result is measured against an explicit baseline.

**Full results: [docs/RESULTS.md](docs/RESULTS.md)** · Plan: [docs/PROJECT_OUTLINE.md](docs/PROJECT_OUTLINE.md)

## Headline results

| Question | Answer |
|---|---|
| Is the model right? | Matches the MATPOWER reference output of RTS-GMLC to its print precision (4.8e-4 pu, $0.0014/h on $225,806/h); the 5-bus port matches to 2e-10 pu |
| How much solar can each bus take? | 17-1,966 MW firm (115x spread); median 38 MW at 138 kV vs. 574 MW at 230 kV; a DC-only screen overstates it 1-37x bus by bus (median 5x) |
| What limits it? | Mostly voltage rise at light load (75 % of buses), then thermal ratings |
| An 870 MW project at a bus with 245 MW firm capacity? | Flexible interconnection: 1.18 % curtailment for **3.5x the energy** of a firm-sized plant |
| Cheapest way to recover the curtailment? | Evaluated one at a time: **inverter Volt-VAR** (a settings change) recovers 9.2 GWh/yr; among capital options, reconductoring the limiting line costs **$135/MWh**, 3x less than the cheapest battery ($418/MWh at $150/kWh) |
| Does Volt-VAR help elsewhere? | At the 35 buses limited by their own voltage, median firm capacity rises **65 -> 468 MW** (32 of 35 gain more than 10 %, none worse) |
| Can N-1 screening be fast and complete? | Pure DC (LODF) screening finds only 25 % of AC overloads; adding base-case reactive flow finds **99 % at 3.8x** the speed of full AC |

## Phases

| Phase | Topic | Script |
|---|---|---|
| 0 | Python port of a 5-bus 24-h study, validated against MATPOWER ([notes](docs/phase0_results.md)) | `phase0_reproduce.py` |
| 1 | RTS-GMLC model validation + full-year DC-OPF / AC baseline, MVA calibration of DC limits | `phase1_rts_baseline.py` |
| 2 | Solar hosting-capacity map, all 73 buses (OPF-based, AC-verified) | `phase2_hosting_capacity.py` |
| 3 | Automated N-1 screening: DC, hybrid, full AC | `phase3_n1_screening.py` |
| 4 | Mitigations for an 870 MW project: flexible interconnection, Volt-VAR, line upgrade, battery | `phase4_mitigations.py` |
| 4b | Volt-VAR at voltage-limited buses | `phase4b_volt_var_map.py` |
| 5 | Sensitivities: load growth, solar penetration, battery cost | `phase5_sensitivities.py` |

## Quick start

```bash
python -m venv .venv
.venv\Scripts\activate            # Windows  (source .venv/bin/activate on Linux/macOS)
pip install -e .[dev] numba
python scripts/download_rts_gmlc.py     # about 3.6 MB from GridMod/RTS-GMLC
pytest
python scripts/run_all.py               # all phases, several hours on a laptop (2 worker processes)
```

## Repository layout

```
src/gridstudy/   RTS-GMLC loader, DC-OPF (PTDF LP, LMPs, max injection), hourly AC power flow,
                 hosting capacity with AC feedback, N-1 contingency analysis, storage / upgrade economics
scripts/         one script per phase; writes results/ and figures/
tests/           model validation, LMP finite-difference, PTDF, dispatch checks
docs/            outline and results
```

## Tools

pandapower (AC power flow) · SciPy/HiGHS and Pyomo (linear programs) · NumPy / pandas · matplotlib · pytest

## Notes

- RTS-GMLC data: https://github.com/GridMod/RTS-GMLC (downloaded, not redistributed).
- The 5-bus validation case is MATPOWER `case5` (F. Li and R. Bo, 2010) with documented modifications.
- Results were checked by independent numeric and code reviews; issues they found (an MW-vs-MVA
  artifact in the AC checks, a headroom lower bound used beyond its range) were fixed and all phases rerun.
