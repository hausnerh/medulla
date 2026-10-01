#!/usr/bin/env python3
"""
Write the per-event weight file list for run_systematics from a selection toml.

run_systematics reads GENIE/flux universe weights from `[input].weights`, which
it hands to a ROOT TChain. A TChain wildcard can't express the osc production's
nested */*/*.root layout across two productions (run2 + run4), but run_systematics
also accepts a .txt file listing one input per line. This builds that list from
the SAME globs the selection used for its MC sample, so the weights always come
from exactly the files that were selected.

Usage:
  make_weights_list.py SELECTION_TOML OUTPUT_TXT [--sample cvext]
"""
import argparse
import sys
from glob import glob

import toml

ap = argparse.ArgumentParser()
ap.add_argument("selection_toml")
ap.add_argument("output_txt")
ap.add_argument("--sample", default="cvext", help="sample whose files carry the weights")
args = ap.parse_args()

samples = [s for s in toml.load(args.selection_toml).get("sample", []) if s.get("name") == args.sample]
if not samples:
    sys.exit(f"error: no [[sample]] named '{args.sample}' in {args.selection_toml}")
globs = []
for s in samples:
    globs += s["path"] if isinstance(s["path"], list) else [s["path"]]
files = sorted({f for g in globs for f in glob(g)})
if not files:
    sys.exit(f"error: sample '{args.sample}' globs matched no files: {globs}")
with open(args.output_txt, "w") as fh:
    fh.write("\n".join(files) + "\n")
print(f"wrote {len(files)} weight files ({len(globs)} glob(s)) -> {args.output_txt}")
