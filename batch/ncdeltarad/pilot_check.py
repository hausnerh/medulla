#!/usr/bin/env python3
"""
Validate a pilot of the enhanced-radiative NC Δ sample from the CSV written
by dump_delta_decays.C. Checks, in order:

  1. The patched table was used: the observed radiative fraction per W bin
     matches delta_br.py's enhanced prediction (and is far from nominal).
  2. W in CAF truth (simb::MCNeutrino::W, i.e. SRTrueInteraction::w) equals
     the Δ invariant mass GENIE decays with. The weight is a function of W,
     so if they differ the weight must be computed from the Δ mass instead.
  3. The CAF-level radiative tag (a status-1 photon above 30 MeV in an NC Δ
     event) agrees with the true decay channel.
  4. Struck nucleon -> Δ charge (proton -> Δ+, neutron -> Δ0).

Exit status is non-zero if any check fails.
    python batch/ncdeltarad/pilot_check.py delta_decays.csv [--br-gamma 0.663]
"""
import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import delta_br as D  # noqa: E402


def load(path):
    rows = list(csv.DictReader(open(path)))
    if not rows:
        raise SystemExit(f'{path}: no NC Δ decays')
    col = lambda k, f=float: np.array([f(r[k]) for r in rows])
    return dict(pdg=col('delta_pdg', int), w_delta=col('w_delta'), w_kine=col('w_kine'),
                hitnuc=col('hitnuc', int), nphot=col('n_photon_status1_gt30mev', int),
                rad=np.array(['22' in r['channel'].split(':') for r in rows]),
                channel=[r['channel'] for r in rows])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('csv')
    ap.add_argument('--br-gamma', type=float, default=D.BR_GAMMA_DEFAULT)
    ap.add_argument('--w-tol', type=float, default=1e-3, help='GeV, |W_kine - W_delta| tolerance')
    a = ap.parse_args()
    d = load(a.csv)
    n = len(d['pdg'])
    ok = True
    print(f'{n} NC Δ decays ({(d["pdg"] == 2214).sum()} Δ+, {(d["pdg"] == 2114).sum()} Δ0)')

    # 1. radiative fraction vs W against the enhanced and nominal predictions
    edges = np.array([1.07, 1.12, 1.18, 1.24, 1.30, 1.40, 1.60, 2.5])
    print('\n1. radiative fraction vs W (observed ± binomial, enhanced / nominal prediction)')
    chi2, ndf = 0.0, 0
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (d['w_delta'] >= lo) & (d['w_delta'] < hi)
        k = sel.sum()
        if k < 5:
            continue
        obs = d['rad'][sel].mean()
        pred = np.mean([D.p_gamma(p, w, a.br_gamma) for p, w in zip(d['pdg'][sel], d['w_delta'][sel])])
        nom = np.mean([D.p_gamma(p, w) for p, w in zip(d['pdg'][sel], d['w_delta'][sel])])
        err = np.sqrt(max(pred * (1 - pred), 1e-6) / k)
        chi2 += ((obs - pred) / err)**2
        ndf += 1
        print(f'   W [{lo:.2f},{hi:.2f}) n={k:5d}  obs {obs:.3f} ± {err:.3f}   enhanced {pred:.3f}   nominal {nom:.4f}')
    good1 = ndf > 0 and chi2 / ndf < 3
    print(f'   chi2/ndf = {chi2:.1f}/{ndf}  ->  {"OK" if good1 else "FAIL (table not applied, or a different BR)"}')
    ok &= good1

    # 2. CAF W vs the Δ mass GENIE decays with
    dw = d['w_kine'] - d['w_delta']
    frac = np.mean(np.abs(dw) < a.w_tol)
    good2 = frac > 0.99
    print(f'\n2. |W_kine - W_delta| < {a.w_tol} GeV for {frac:.1%} '
          f'(median {np.median(dw):+.4f}, 1-99% [{np.percentile(dw, 1):+.4f}, {np.percentile(dw, 99):+.4f}] GeV)  ->  '
          f'{"OK: weight from CAF w" if good2 else "MISMATCH: compute the weight from the Δ mass"}')
    ok &= good2

    # 3. CAF-level radiative tag vs truth
    tag = d['nphot'] >= 1
    mis_rad = np.mean(~tag[d['rad']]) if d['rad'].any() else 0
    mis_non = np.mean(tag[~d['rad']]) if (~d['rad']).any() else 0
    good3 = mis_rad < 0.01 and mis_non < 0.01
    print(f'\n3. photon tag: radiative untagged {mis_rad:.2%}, non-radiative tagged {mis_non:.2%}  ->  '
          f'{"OK" if good3 else "FAIL: revisit the 30 MeV photon tag"}')
    ok &= good3

    # 4. struck nucleon gives the Δ charge
    good4 = np.all(D.delta_pdg_from_hitnuc(d['hitnuc']) == d['pdg'])
    print(f'\n4. hitnuc -> Δ charge consistent for all decays: {"OK" if good4 else "FAIL"}')
    ok &= bool(good4)

    w = D.weight(d['rad'], d['w_delta'], d['pdg'], a.br_gamma)
    print(f'\nweights: radiative median {np.median(w[d["rad"]]):.4f}, Nπ median {np.median(w[~d["rad"]]):.3f}, '
          f'max {w.max():.2f}; effective radiative statistics gain over nominal x{1 / np.mean(w[d["rad"]]):.0f}')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
