#!/usr/bin/env python3
"""
Diagnose the systematic band on a low-stat sideband bin.

The grey 'Total syst.' band is built purely from MC universe weights
(spineplot/systematic.py). For a bin holding only 1-2 MC events, the FRACTIONAL
band equals the spread of those events' universe weights, so a single event with
a large multisim weight swing inflates the absolute band to ~1 candidate even
on a ~0.05 nominal. This script finds the culprit: it bins the CV tree on the
plotted variable, picks a bin, and for every systematic knob prints the events
whose universe weights swing the most (std across universes, and min/max).

Run inside the SL7 container with uproot on PATH (your plot venv), e.g.:
  python3 batch/diagnose_sys_band.py build/output_gOre_1g1p_sys.root \
      --tree selected_cc_Xg1p_nm1_delta_mass --var reco_delta_mass_P \
      --lo 1350 --hi 1500
"""
import argparse
import numpy as np
import uproot

ap = argparse.ArgumentParser()
ap.add_argument("root", help="systematics ROOT file (local path)")
ap.add_argument("--key", default="full", help="events/<key> directory (default: full)")
ap.add_argument("--tree", default="selected_cc_Xg1p_nm1_delta_mass",
                help="CV tree name under events/<key>")
ap.add_argument("--var", default="reco_delta_mass_P", help="plotted variable branch")
ap.add_argument("--lo", type=float, required=True, help="low edge of the window to inspect")
ap.add_argument("--hi", type=float, required=True, help="high edge of the window to inspect")
ap.add_argument("--systs", default="multisimTree,multisigmaTree,variationTree",
                help="comma-separated systematic tree suffixes")
ap.add_argument("--top", type=int, default=8, help="show this many worst knobs")
args = ap.parse_args()

f = uproot.open(args.root)
d = f[f"events/{args.key}"]

cv = d[args.tree].arrays([args.var], library="np")[args.var]
sel = (cv > args.lo) & (cv < args.hi)
idx = np.where(sel)[0]
print(f"CV tree events/{args.key}/{args.tree}: {len(cv)} rows; "
      f"{sel.sum()} in window ({args.lo}, {args.hi}) on {args.var}")
print(f"  nominal (unweighted) count in window = {sel.sum()}  "
      f"(row indices {list(idx)})")
if sel.sum() == 0:
    raise SystemExit("no events in window — widen --lo/--hi")

for suffix in args.systs.split(","):
    tname = f"{args.tree}_{suffix}"
    if tname not in [k.split(";")[0] for k in d.keys()]:
        print(f"\n[{suffix}] tree not present — skipping")
        continue
    wt = d[tname]
    knobs = [k for k in wt.keys() if k not in ("Run", "Subrun", "Evt")]
    print(f"\n[{suffix}] {tname}: {len(knobs)} knob branch(es)")
    rows = []
    for knob in knobs:
        arr = wt[knob].array(library="np")          # per-event vector of weights
        arr = np.stack(arr)[idx, :]                 # (n_window_events, n_universe_pts)
        # per-event spread across universes/sigma-points; report the worst event
        ev_std = arr.std(axis=1)
        j = int(np.argmax(ev_std))
        rows.append((knob, ev_std[j], arr[j].min(), arr[j].max(), arr.shape[1]))
    rows.sort(key=lambda r: r[1], reverse=True)
    print(f"  {'knob':<32} {'max ev std':>11} {'min w':>9} {'max w':>9} {'npts':>5}")
    for knob, s, lo, hi, npts in rows[: args.top]:
        print(f"  {knob:<32} {s:>11.3f} {lo:>9.3f} {hi:>9.3f} {npts:>5}")

print("\nInterpretation: a knob with 'max ev std' >> 1 (or min/max weights far "
      "from ~1) is a single event whose reweight is pathological; that one event "
      "drives the tail band. Multisim knobs use these directly; multisigma knobs "
      "are interpolated in spineplot, so their effective spread is comparable.")
