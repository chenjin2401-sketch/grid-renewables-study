"""Download the RTS-GMLC files used by this project into data/rts_gmlc/.

Source: https://github.com/GridMod/RTS-GMLC (branch master). Only the network
tables, the MATPOWER case with its reference output, and the hourly
day-ahead time series are fetched (about 3.6 MB).
"""
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "data" / "rts_gmlc"
BASE = "https://raw.githubusercontent.com/GridMod/RTS-GMLC/master/RTS_Data/"
FILES = [
    "SourceData/bus.csv",
    "SourceData/branch.csv",
    "SourceData/gen.csv",
    "SourceData/storage.csv",
    "SourceData/dc_branch.csv",
    "SourceData/timeseries_pointers.csv",
    "SourceData/README.md",
    "FormattedData/MATPOWER/RTS_GMLC.m",
    "FormattedData/MATPOWER/MATPOWER-out.txt",
    "timeseries_data_files/Load/DAY_AHEAD_regional_Load.csv",
    "timeseries_data_files/PV/DAY_AHEAD_pv.csv",
    "timeseries_data_files/RTPV/DAY_AHEAD_rtpv.csv",
    "timeseries_data_files/WIND/DAY_AHEAD_wind.csv",
    "timeseries_data_files/Hydro/DAY_AHEAD_hydro.csv",
]


def main():
    for rel in FILES:
        out = DEST / rel
        if out.exists():
            continue
        out.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(BASE + rel, timeout=60) as r:
            out.write_bytes(r.read())
        print(f"{out.stat().st_size:>9,d}  {rel}")
    print("done:", DEST)


if __name__ == "__main__":
    main()
