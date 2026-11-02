# Grid Renewables Study — Project Outline

**One-line pitch:** an open, reproducible *interconnection study* for solar + battery storage on a
realistic transmission test system: how much solar can each bus take, what breaks under
contingencies, and when a battery beats a line upgrade, with every result measured against a
clear baseline.

**Target roles:** utility transmission planning / operations, renewable and storage developers.

---

## 1. The question the project answers

> A developer wants to connect a solar plant (and possibly a battery) to the grid.
> Where can it connect, how large can it be, what violations does it cause (normal and N-1),
> and what is the cheapest way to fix them: inverter controls, a battery, curtailment, or a
> line upgrade?

This is the question a utility planning engineer answers in an interconnection study, and the
one a developer has to understand to get a project built.

## 2. System and data

| Item | Choice | Why |
|---|---|---|
| Validation system | 5-bus system (MATPOWER `case5` with documented modifications) | Small enough to check by hand; existing MATPOWER results to validate against |
| Main system | **RTS-GMLC** (73 buses, open source, GridMod/RTS-GMLC on GitHub) | Realistic grid with full-year hourly load, solar and wind data |
| Cost data | NREL Annual Technology Baseline (battery, solar); ISO transmission cost guides (line upgrades) | Public, citable sources |
| Tools | Python, pandapower (AC power flow), Pyomo + HiGHS (optimization), pandas, matplotlib, pytest, Streamlit | Open-source stack anyone can rerun |

## 3. Phases, deliverables and the numbers each phase produces

Every phase ends with (a) code + tests, (b) figures, (c) a short results note with the
quantified outcome, and (d) a learning checkpoint (questions you must be able to answer).

| Phase | Work | Baseline → comparison | Headline numbers |
|---|---|---|---|
| **0. Python port + validation** | Re-implement the 24-h 5-bus study (fixed schedule and economic dispatch) in pandapower/Pyomo | MATPOWER results | Max deviation in V (pu), branch MW, losses, cost |
| **1. Scale-up** | Load RTS-GMLC into pandapower; time-series AC power flow on representative days (summer peak, spring high-solar/low-load, winter) | Published RTS-GMLC base case (convergence, flows) | System size, hours simulated, runtime, baseline violation count |
| **2. Hosting capacity** | For each candidate bus, increase solar MW until the first thermal/voltage violation (worst hours, N-0) | No-mitigation hosting capacity | MW per bus (map); best vs. worst site ratio; binding constraint per bus |
| **3. N-1 screening** | Automated contingency analysis (every line/transformer outage) at critical hours, before/after solar | Full AC for every case vs. DC pre-screen + AC check | Contingencies screened, critical overloads found, runtime speed-up, missed-violation rate |
| **4. Mitigations** | (a) Inverter Volt-VAR (IEEE 1547-style), (b) flexible interconnection with capped curtailment, (c) battery sized by multi-period DC-OPF, (d) line upgrade | Phase 2/3 results without mitigation | Hosting-capacity gain (MW, %), overloads removed, curtailment (MWh, %), $ cost vs. line upgrade, payback |
| **5. Sensitivities** | Solar penetration 10–60 %, load growth 0–20 %, battery cost −30 %, seasons | Phase 4 base case | Tipping points ("violations start at X % penetration", "battery wins below $Y/kWh") |
| **6. Packaging** | README with results, Streamlit dashboard, short technical write-up | — | One-command reproduction; demoable app |

**Minimum showcase:** Phases 0–3 (planning story). **Full story:** through Phase 5.

## 4. Experiment matrix (baseline vs. intervention)

| Experiment | Baseline | Intervention | Metric |
|---|---|---|---|
| E1 Model accuracy | MATPOWER | pandapower port | max \|ΔV\|, max \|ΔP\|, Δloss, Δcost |
| E2 Siting | Worst / average bus | Best bus | Hosting capacity (MW), losses |
| E3 Reactive support | Inverters at unity pf | Volt-VAR | Hosting capacity gain (%), voltage violations |
| E4 Flexible interconnection | Firm (no curtailment) | ≤ 1 %, 3 %, 5 % annual curtailment | Extra MW connected per % curtailed |
| E5 Contingency risk | Pre-solar grid | With solar (± battery) | N-1 overload count, worst loading (%) |
| E6 Non-wires alternative | Line reconductoring | Battery (MW/MWh from optimization) | $ per MW of relief, payback, side benefits |
| E7 Battery sizing | No battery | 0 → N MWh sweep | Marginal value curve, knee point |
| E8 Screening speed | Full AC N-1 | DC pre-screen + AC verification | Runtime ×, missed violations |

Rules: the baseline is what industry would actually do (no straw men); every number is produced
by a script; report absolute and relative values ("+40 %, 50 → 70 MW").

## 5. Repository layout

```
grid-renewables-study/
├── src/gridstudy/        # library code (networks, power flow, dispatch, studies)
├── scripts/              # one script per phase/experiment, writes results/ and figures/
├── tests/                # pytest: validation tolerances, sanity checks
├── data/                 # small inputs; large external data downloaded by a script
├── results/              # CSV outputs (generated)
├── figures/              # PNG outputs (generated)
└── docs/                 # this outline, results notes, learning notes
```

## 6. Timeline (adjust around exams)

| When | Phase |
|---|---|
| Nov 2026 (wk 1) | Phase 0 |
| Nov 2026 (wk 2–4) | Phase 1 |
| Dec 2026 – winter break | Phases 2–3 |
| Jan 2027 | Phase 4 |
| Feb 2027 | Phases 5–6, publish |

## 7. Learning checkpoints (what you must be able to explain)

- **Phase 0:** Newton-Raphson power flow; PQ / PV / slack buses; why the slack picks up losses;
  merit-order dispatch and the incremental cost λ; how line charging affects reactive flow.
- **Phase 1:** per-unit system; transformer modeling; why convergence can fail and what to try.
- **Phase 2:** thermal vs. voltage limits; why hosting capacity depends on location and hour;
  why low load + high solar is the worst case.
- **Phase 3:** N-1 criterion (NERC TPL-001 idea); DC power flow assumptions and when they break;
  PTDF / LODF basics.
- **Phase 4:** Volt-VAR curves and reactive capability; multi-period DC-OPF with storage
  (state-of-charge, charge/discharge efficiency); non-wires alternatives.
- **Phase 5:** sensitivity vs. scenario analysis; how to present uncertainty.

## 8. Resume bullet templates (fill with real numbers when each phase is done)

- Built an open-source Python interconnection-study framework (pandapower, Pyomo) validated
  against MATPOWER to within **[x]** pu / **[y]** MW.
- Mapped solar hosting capacity across **[N]** buses of the 73-bus RTS-GMLC system
  (**[a–b]** MW); optimal siting gave **[k]×** the capacity of the worst site.
- Automated N-1 screening of **[M]** contingencies; DC pre-screening cut runtime **[r]×** with
  **[0]** missed violations.
- Showed a **[P MW / E MWh]** battery resolved **[K]** overloads at **[c]%** lower cost than
  reconductoring (**[n]**-year payback); Volt-VAR plus ≤3 % curtailment raised hosting capacity
  **[g]%**.

## 9. Ground rules

- Public repo contains only open data and original code: no course handouts, instructor code,
  or course reports.
- AI-assisted implementation is disclosed in the README; design decisions and analysis are the
  author's own.
