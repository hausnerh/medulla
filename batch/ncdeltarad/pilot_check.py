#!/usr/bin/env python3
"""
Validate the enhanced-radiative NC Δ sample from the gen-level sidecar CSV
written by dump_delta_decays.C. The orchestrator runs this as a gate after the
gen stage; it can also be run by hand on a pilot.

Required (non-zero exit if they fail):
  1. The patched table was used: the radiative fraction per W bin matches
     delta_br.py's enhanced prediction (and is far from nominal).
  2. The sidecar can be joined to the CAF / medulla output: (run, subrun,
     event) is unique, and every Δ has a two-body decay.

Informational (the weight is taken from the sidecar, not from CAF truth):
  3. CAF w (simb::MCNeutrino::W) vs the Δ mass GENIE decayed with. nugen fills
     MCNeutrino::W assuming a free nucleon at rest, so it differs from the
     Δ mass by Fermi motion and binding (pilot: median +31 MeV, up to +0.49 GeV).
  4. CAF-level radiative tag (status-1 photon > 30 MeV) vs the true channel.
  5. struck nucleon vs Δ charge.

    python batch/ncdeltarad/pilot_check.py delta_decays.csv [--br-gamma 0.663]
"""
import argparse
import csv
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import delta_br as D  # noqa: E402


def load(path):
    rows = list(csv.DictReader(open(path)))
    if not rows:
        raise SystemExit(f'{path}: no NC Δ decays')
    col = lambda k, f=float: np.array([f(r[k]) for r in rows])
    return dict(run=col('run', int), subrun=col('subrun', int), event=col('event', int),
                pdg=col('delta_pdg', int), w_delta=col('w_delta'), w_kine=col('w_kine'),
                hitnuc=col('hitnuc', int), nphot=col('n_photon_status1_gt30mev', int),
                nd=col('n_daughters', int), rad=np.array(['22' in r['channel'].split(':') for r in rows]),
                channel=[r['channel'] for r in rows])


def radiative_fraction_check(d, br_gamma):
    edges = np.array([1.07, 1.12, 1.18, 1.24, 1.30, 1.40, 1.60, 2.5])
    print('\n1. radiative fraction vs W (observed ± binomial, enhanced / nominal prediction)')
    chi2, ndf = 0.0, 0
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (d['w_delta'] >= lo) & (d['w_delta'] < hi)
        k = sel.sum()
        if k < 5:
            continue
        obs = d['rad'][sel].mean()
        pred = np.mean([D.p_gamma(p, w, br_gamma) for p, w in zip(d['pdg'][sel], d['w_delta'][sel])])
        nom = np.mean([D.p_gamma(p, w) for p, w in zip(d['pdg'][sel], d['w_delta'][sel])])
        err = np.sqrt(max(pred * (1 - pred), 1e-6) / k)
        chi2 += ((obs - pred) / err)**2
        ndf += 1
        print(f'   W [{lo:.2f},{hi:.2f}) n={k:5d}  obs {obs:.3f} ± {err:.3f}   enhanced {pred:.3f}   nominal {nom:.4f}')
    good = ndf > 0 and chi2 / ndf < 3
    print(f'   chi2/ndf = {chi2:.1f}/{ndf}  ->  {"OK" if good else "FAIL (table not applied, or a different BR)"}')
    return good


def join_check(d):
    keys = Counter(zip(d['run'], d['subrun'], d['event']))
    dup = sum(1 for v in keys.values() if v > 1)
    bad_nd = int((d['nd'] != 2).sum())
    good = dup == 0 and bad_nd == 0
    print(f'\n2. join key (run, subrun, event): {len(keys)} unique, {dup} duplicated; '
          f'Δ without a two-body decay: {bad_nd}  ->  {"OK" if good else "FAIL"}')
    return good


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('csv')
    ap.add_argument('--br-gamma', type=float, default=D.BR_GAMMA_DEFAULT)
    a = ap.parse_args()
    d = load(a.csv)
    print(f'{len(d["pdg"])} NC Δ decays ({(d["pdg"] == 2214).sum()} Δ+, {(d["pdg"] == 2114).sum()} Δ0)')

    ok = radiative_fraction_check(d, a.br_gamma)
    ok &= join_check(d)

    dw = d['w_kine'] - d['w_delta']
    print(f'\n3. [info] CAF w - Δ mass: median {np.median(dw):+.4f}, 1-99% '
          f'[{np.percentile(dw, 1):+.4f}, {np.percentile(dw, 99):+.4f}] GeV (weight uses the Δ mass)')
    tag = d['nphot'] >= 1
    mis_rad = np.mean(~tag[d['rad']]) if d['rad'].any() else 0
    mis_non = np.mean(tag[~d['rad']]) if (~d['rad']).any() else 0
    print(f'\n4. [info] photon tag: radiative untagged {mis_rad:.2%}, non-radiative tagged {mis_non:.2%}')
    tab = Counter(zip(d['hitnuc'], d['pdg']))
    print('\n5. [info] struck nucleon x Δ (count):  ' +
          ', '.join(f'hitnuc {h} -> {p}: {n}' for (h, p), n in sorted(tab.items())))

    w = D.weight(d['rad'], d['w_delta'], d['pdg'], a.br_gamma)
    print(f'\nweights: radiative median {np.median(w[d["rad"]]):.4f}, Nπ median {np.median(w[~d["rad"]]):.3f}, '
          f'max {w.max():.2f}; effective radiative statistics gain over nominal x{1 / np.mean(w[d["rad"]]):.0f}')
    print(f'\n{"PASS" if ok else "FAIL"}')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
