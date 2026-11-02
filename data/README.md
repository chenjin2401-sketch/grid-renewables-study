# Data

| Folder | Contents | In git? |
|---|---|---|
| `data/private/` | 24-hour load profile and MATPOWER reference results from the author's ECE 451 coursework, used only to validate the Phase 0 port | No (`.gitignore`) |
| `data/rts_gmlc/` | RTS-GMLC network and time series, downloaded by a script (Phase 1) | No (downloaded) |

Create `data/private/` from the coursework folder with:

```
python scripts/import_course_reference.py "<path to the MATLAB project folder>"
```

Phase 0 tests are skipped when `data/private/` is missing, so the public repo still runs.
