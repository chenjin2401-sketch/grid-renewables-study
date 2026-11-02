# Grid Renewables Study

An open, reproducible **interconnection study for solar + battery storage** on transmission test
systems: how much solar each bus can take, what fails under N-1 contingencies, and when a battery
beats a line upgrade. Every result is measured against an explicit baseline.

Full plan: [docs/PROJECT_OUTLINE.md](docs/PROJECT_OUTLINE.md)

## Status

| Phase | Topic | Status | Key result |
|---|---|---|---|
| 0 | Python port of a 5-bus 24-h study, validated against MATPOWER | ✅ done | Matches MATPOWER within 2e-10 pu / 1e-6 MW ([notes](docs/phase0_results.md)) |
| 1 | Scale-up to the 73-bus RTS-GMLC system with real time series | ⏳ next | |
| 2 | Solar hosting-capacity map | | |
| 3 | Automated N-1 screening | | |
| 4 | Mitigations: Volt-VAR, flexible interconnection, battery vs. line upgrade | | |
| 5 | Sensitivities and tipping points | | |
| 6 | Dashboard and write-up | | |

## Quick start

```bash
python -m venv .venv
.venv\Scripts\activate            # Windows  (source .venv/bin/activate on Linux/macOS)
pip install -e .[dev]
pytest                            # Phase 0 validation tests skip if data/private/ is absent
python scripts/phase0_reproduce.py
```

## Repository layout

```
src/gridstudy/   networks, hourly AC power flow, dispatch LP, validation helpers
scripts/         one script per phase; writes results/ and figures/
tests/           pytest checks
data/            see data/README.md (coursework-derived inputs stay local)
docs/            outline and per-phase results notes
```

## Tools

pandapower (AC power flow) · Pyomo + HiGHS (optimization) · pandas · matplotlib · pytest

## Notes

- The 5-bus validation case is MATPOWER `case5` (PJM 5-bus system, F. Li and R. Bo, 2010) with
  documented modifications.
- Implementation was AI-assisted. The study design, modeling choices and analysis are the
  author's.
