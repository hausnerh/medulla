#!/usr/bin/env python3
"""Report POT/Livetime histograms and _exposure-tree sums per sample directory.

Usage: diag_exposure.py FILE_OR_GLOB [...]
For each events/<sample> it prints whether the POT/Livetime histograms exist
(and their values, all cycles listed), plus the summed `pot`/`livetime` of the
first <tree>_exposure tree, so a missing histogram can be compared with what
the selection's exposure trees recorded.
"""
import glob
import sys

import uproot


def main():
    files = sorted({f for p in sys.argv[1:] for f in (glob.glob(p) or [p])})
    for fn in files:
        f = uproot.open(fn)
        if 'events' not in {k.split(';')[0] for k in f.keys(recursive=False)}:
            print(f'{fn}: no events/ directory'); continue
        for s in sorted({k.split(';')[0] for k in f['events'].keys(recursive=False)}):
            d = f['events/' + s]
            keys = d.keys(recursive=False)
            names = {k.split(';')[0] for k in keys}
            pots = [k for k in keys if k.split(';')[0] == 'POT']
            pot = d['POT'].values()[0] if 'POT' in names else 'MISSING'
            lt = d['Livetime'].values()[0] if 'Livetime' in names else 'MISSING'
            exp = sorted(n for n in names if n.endswith('_exposure'))
            es = ''
            if exp:
                a = d[exp[0]].arrays(['pot', 'livetime'], library='np')
                es = f"  {exp[0]}: n={len(a['pot'])} sum(pot)={a['pot'].sum():.4g} sum(livetime)={a['livetime'].sum():.4g}"
            print(f"{fn.split('/')[-1]:40s} {s:9s} POT={pot} ({len(pots)} cycles) Livetime={lt}{es}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
