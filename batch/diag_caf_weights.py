#!/usr/bin/env python3
"""How many reweight groups does each true neutrino carry in a set of flat CAFs?

run_systematics reads knob `index` i from rec.mc.nu.wgt[i] and aborts with
"WeightReader: Index out of range in 'get_nuniv()'" when a neutrino has
<= i groups. This tallies rec.mc.nu.wgt..length over every neutrino and
prints the neutrinos that carry fewer groups than the most common count.

Usage: diag_caf_weights.py 'GLOB_OR_FILE' [...] [--need 189]
"""
import argparse
import glob
import sys
from collections import Counter

import numpy as np
import uproot

BR = ['rec.hdr.run', 'rec.hdr.subrun', 'rec.hdr.evt',
      'rec.mc.nu.wgt..length', 'rec.mc.nu.pdg', 'rec.mc.nu.E', 'rec.mc.nu.iscc']


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paths', nargs='+')
    ap.add_argument('--need', type=int, default=189, help='groups the knob list needs (max index + 1)')
    ap.add_argument('--show', type=int, default=20, help='print at most this many short neutrinos')
    args = ap.parse_args()

    files = sorted({f for p in args.paths for f in (glob.glob(p) or [p])})
    counts, short, shown, nnu = Counter(), Counter(), 0, 0
    for fn in files:
        t = uproot.open(fn)['recTree']
        have = [b for b in BR if b in t.keys()]
        a = t.arrays(have, library='ak')
        for ev in range(len(a)):
            nw = a['rec.mc.nu.wgt..length'][ev]
            for i in range(len(nw)):
                nnu += 1
                n = int(nw[i])
                counts[n] += 1
                if n < args.need:
                    short[fn.split('/')[-1]] += 1
                    if shown < args.show:
                        shown += 1
                        extra = ' '.join(f"{b.split('.')[-1]}={a[b][ev][i]}" for b in have[4:])
                        print(f"  SHORT nwgt={n:3d} run={a['rec.hdr.run'][ev]} subrun={a['rec.hdr.subrun'][ev]} "
                              f"evt={a['rec.hdr.evt'][ev]} nu#{i} {extra}  [{fn.split('/')[-1]}]")
    print(f"\n{len(files)} files, {nnu} neutrinos. wgt-group count -> #neutrinos:")
    for n, c in sorted(counts.items()):
        flag = '   <-- fewer than --need' if n < args.need else ''
        print(f"  {n:4d}: {c}{flag}")
    if short:
        print(f"\n{sum(short.values())} short neutrinos in {len(short)} file(s); worst files:")
        for fn, c in short.most_common(10):
            print(f"  {c:6d}  {fn}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
