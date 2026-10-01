#!/usr/bin/env python3
"""Truth-matching bookkeeping for a hadd'd medulla SELECTION output.

For every events/<sample>/<tree> it reports how the selected reco interactions
split by truth match, mirroring what run_systematics' add_weights does next:

  matched_nu  true_neutrino_id >= 0          -> events/full (one row per nu!)
  cosmic      true_neutrino_id <  0          -> <tree>_nonmatched (not plotted)
  no_match    true_neutrino_id is NaN        -> <tree>_nonmatched (not plotted)
              (empty SPINE match_ids: e.g. an overlay data cosmic, no truth)
  dup_rows    extra rows sharing (Run, Subrun, Evt, true_neutrino_id) with
              another selected row: one true neutrino reconstructed as >1
              selected reco interaction (SPINE reco->truth matching is
              many-to-one). add_weights keeps only the FIRST such row and
              silently drops the rest, while data counts every fragment.

Usage: diag_truth_match.py build/output_gOre_1g1p_signal.root [--trees PATTERN]
"""
import argparse
import re
import sys

import numpy as np
import uproot


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('selection_root')
    ap.add_argument('--trees', default=r'^selected_', help='regex on tree name (default: ^selected_)')
    args = ap.parse_args()

    f = uproot.open(args.selection_root)
    pat = re.compile(args.trees)
    hdr = f"{'sample/tree':55s} {'rows':>8s} {'matched_nu':>10s} {'cosmic':>7s} {'no_match':>8s} {'dup_rows':>8s} {'dup%':>6s}"
    print(hdr)
    print('-' * len(hdr))
    for sample in sorted(k.rstrip(';1').split('/')[-1] for k in f['events'].keys(recursive=False)):
        d = f[f'events/{sample}']
        for tname in sorted({k.split(';')[0] for k in d.keys(recursive=False)}):
            if not pat.search(tname):
                continue
            t = d[tname]
            if 'true_neutrino_id' not in t.keys():
                continue
            a = t.arrays(['Run', 'Subrun', 'Evt', 'true_neutrino_id'], library='np')
            nu = a['true_neutrino_id']
            n = len(nu)
            nan = np.isnan(nu)
            neg = ~nan & (nu < 0)
            ok = ~nan & (nu >= 0)
            keys = np.stack([a['Run'][ok], a['Subrun'][ok], a['Evt'][ok], nu[ok].astype(np.int64)], axis=1)
            n_unique = len(np.unique(keys, axis=0)) if len(keys) else 0
            dup = int(ok.sum()) - n_unique
            pct = 100.0 * dup / ok.sum() if ok.sum() else 0.0
            print(f"{sample + '/' + tname:55s} {n:8d} {int(ok.sum()):10d} {int(neg.sum()):7d} {int(nan.sum()):8d} {dup:8d} {pct:6.2f}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
