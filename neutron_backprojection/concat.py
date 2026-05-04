#!/usr/bin/env python3
"""
Stream files with pandas, renumber Event_ID across files, write TSV.
Usage: python concat_eventid_pandas_stream.py "test*.dat" out.tsv
"""
import sys, glob, re
import pandas as pd

def nk(s): return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]

if len(sys.argv) < 3:
    print("Usage: script.py '<glob>' out.tsv", file=sys.stderr); sys.exit(1)

pattern, out = sys.argv[1], sys.argv[2]
files = sorted(glob.glob(pattern), key=nk)
if not files: raise SystemExit("No files matched")

open(out, "w").close()  # clear output file

last_max = 0
first_file = True
chunksize = 100_000

for f in files:
    it = pd.read_csv(f, sep=None, engine="python", chunksize=chunksize)
    for ci, chunk in enumerate(it):
        if 'EventID' not in chunk.columns:
            raise SystemExit(f"'EventID' column not found in {f}")

        chunk['EventID'] = pd.to_numeric(chunk['EventID'], errors='raise').astype(int)

        if ci == 0:
            actual_first = int(chunk['EventID'].iloc[0])
            delta = 0 if first_file else (last_max + 1 - actual_first)

        chunk['EventID'] += delta
        last_max = max(last_max, int(chunk['EventID'].max()))

        chunk.to_csv(
            out,
            sep="\t",
            index=False,
            header=(first_file and ci == 0),
            mode='a'
        )

        print(f"wrote chunk from {f}")

    first_file = False

print(f"Done. Final EventID = {last_max}")