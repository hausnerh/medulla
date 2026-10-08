#!/usr/bin/env python3
"""Off-beam gate count for samples whose header livetime is empty (ReCAF2026).

Usage: offbeam_gates.py EXPOSURE_ROOT SELECTION_ROOT

EXPOSURE_ROOT is the merged output of selection/toml/gOre_osc_run2_data_exposure.toml
(a cut-free, one-row-per-event `events` tree per data sample). SELECTION_ROOT
is the merged selection output whose off-beam Livetime is to be filled.

medulla's data "Livetime" is a gate count (len(bnbinfo) + noffbeambnb + ...);
on-beam it is the BNB spill count. The ReCAF2026 off-beam files have neither,
so the off-beam gates are recovered as sum(gate_delta) (gates since the
previous trigger of the same type). Checks before a number is printed:
  * every data sample has the same number of events as the selection's
    per-event *_exposure tree (i.e. both runs saw the same files);
  * on-beam sum(gate_delta) agrees with the on-beam Livetime the selection
    took from the spill info (same unit as the off-beam number) to 5%.
Last line on success: "OFFBEAM_LIVETIME <value>". Exit 2 if a check fails.
"""
import sys

import numpy as np
import uproot

BRANCHES = ['gate_delta', 'nbnb', 'noffbeambnb', 'unfolded_nbnb', 'is_first_in_subrun']
ONBEAM_TOL = 0.05


def selection_info(f, s):
    d = f['events/' + s]
    names = {k.split(';')[0] for k in d.keys(recursive=False)}
    lt = d['Livetime'].values()[0] if 'Livetime' in names else None
    exp = sorted(n for n in names if n.endswith('_exposure'))
    n = d[exp[0]].num_entries if exp else None
    return lt, n


def main():
    if len(sys.argv) != 3:
        print(__doc__); return 1
    fe, fs = uproot.open(sys.argv[1]), uproot.open(sys.argv[2])
    ok, sums = True, {}
    for s in sorted({k.split(';')[0] for k in fe['events'].keys(recursive=False)}):
        t = fe[f'events/{s}/events']
        keys = set(t.keys())
        # medulla writes event-scope branches as event_<name>; accept either.
        names = {b: (f'event_{b}' if f'event_{b}' in keys else b) for b in BRANCHES}
        ids = [k for k in ('Run', 'Subrun', 'Evt') if k in keys]
        raw = t.arrays(list(names.values()) + ids, library='np')
        a = {b: raw[names[b]] for b in BRANCHES}
        n = len(a['gate_delta'])
        nuniq = len(set(zip(*(raw[k] for k in ids)))) if len(ids) == 3 else n
        gd = a['gate_delta']
        sums[s] = gd.sum()
        lt_sel, n_sel = selection_info(fs, s)
        print(f"{s:9s} events={n} (unique {nuniq}) sum(gate_delta)={gd.sum():.6g} "
              f"[min {gd.min():.0f} median {np.median(gd):.0f} max {gd.max():.0f}, <=0: {(gd <= 0).sum()}] "
              f"sum(nbnb)={a['nbnb'].sum():.6g} sum(unfolded_nbnb)={a['unfolded_nbnb'].sum():.6g} "
              f"sum(noffbeambnb)={a['noffbeambnb'].sum():.6g} subruns={int(a['is_first_in_subrun'].sum())}")
        print(f"{'':9s} selection: Livetime={lt_sel} per-event exposure rows={n_sel}")
        if nuniq != n:
            print(f"  WARN {s}: {n - nuniq} duplicate (Run, Subrun, Evt) rows — a file may be listed twice")
        if n_sel is not None and n_sel != n:
            print(f"  FAIL {s}: {n} events here vs {n_sel} in the selection — the two runs did not process the same files")
            ok = False
    if 'offbeam' not in sums or sums['offbeam'] <= 0:
        print("  FAIL: no off-beam gates (sum(gate_delta) <= 0)"); ok = False
    lt_on, _ = selection_info(fs, 'onbeam') if 'onbeam' in sums else (None, None)
    if lt_on:
        r = sums['onbeam'] / lt_on
        print(f"onbeam cross-check: sum(gate_delta) / selection Livetime = {r:.4f}")
        if abs(r - 1) > ONBEAM_TOL:
            print(f"  FAIL: on-beam gate_delta disagrees with the spill-info gate count by more than {ONBEAM_TOL:.0%} "
                  f"— gate_delta is not the same unit as Livetime here; do not use it for off-beam")
            ok = False
    else:
        print("  FAIL: no on-beam Livetime in the selection to cross-check against"); ok = False
    if not ok:
        return 2
    print(f"OFFBEAM_LIVETIME {sums['offbeam']:.10g}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
