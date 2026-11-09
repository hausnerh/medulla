#!/usr/bin/env python3
"""
Is the ncdeltarad sample what we wanted? Compares its medulla output, weighted
back to nominal by the gen sidecar, with the CV's NC resonant events, tree by
tree (stage1/2/3 and the true signal tree) and truth category by category.

For every category it prints
  cv       CV rows (all modes) and the NC RES part of them, which is what the
           ncdeltarad sample can replace;
  ncdr     raw ncdeltarad rows;
  pred     the ncdeltarad rows weighted and scaled to the CV POT, with their
           statistical error: closure means pred agrees with the CV NC RES
           count (pull = (pred - cv_res) / sqrt(err^2 + cv_res));
  gain     effective statistics relative to the CV: the nominal exposure the
           weighted rows are worth, POT_ncdr * sum(w) / sum(w^2), over POT_cv.
           It uses the weights, so it never credits the sample with more
           statistics than it has (gain 1 = as good as the CV).

Rows are reco interactions (true interactions in the true signal tree), as in
spineplot. The ncdeltarad sample's cosmic rows (category 9) come from the
overlay and are dropped: the CV and off-beam samples already model cosmics.
Rows join the sidecar on (Run, Subrun, Evt); NC RES events with no Δ+/Δ0
decay (heavier resonances) have no sidecar row and weight 1.

Usage
  sample_check.py --ncdr build/gOre_1g1p_signal_ncdeltarad_run2.root \\
      --cv build/output_gOre_1g1p_signal_osc_run2.root \\
      --weights <persistent_base>/gen_sidecar/ncdeltarad_weights.csv [--json out.json]
"""
import argparse
import csv
import json
import math
import sys

import numpy as np

TREES = ['selected_1g1p_stage1', 'selected_1g1p_stage2', 'selected_1g1p_stage3', 'signal_1g1p']
CATS = ['NC Δ→Nγ (1γ1p)', 'NC Δ→Nγ (other topo)', 'Other NC 1γXp post-FSI',
        'NC π0 (Δ res)', 'NC π0 (other)', 'NC π±', 'Other NC', 'CC νe',
        'Other CC', 'Cosmic']
COSMIC = 9
GENIE_RES = 1
BRANCHES = ['Run', 'Subrun', 'Evt', 'true_category', 'true_cc', 'true_interaction_mode',
            'true_baryon_res_code']


def load(path, key, trees=TREES):
    """{'pot': float, 'trees': {tree: {branch: array}}} for sample `key` of a medulla output."""
    import uproot
    f = uproot.open(path)
    samples = [k.split(';')[0] for k in f['events'].keys(recursive=False)]
    if key not in samples:
        raise SystemExit(f'{path}: no sample "{key}" (has {samples})')
    d = f['events'][key]
    pot = float(np.sum(d['POT'].values()))
    out = {}
    for t in trees:
        if t not in d:
            print(f'NOTE  {path}: sample {key} has no tree {t}')
            continue
        have = set(d[t].keys())
        miss = [b for b in BRANCHES if b not in have]
        if miss:
            raise SystemExit(f'{path}: {key}/{t} lacks {miss}')
        out[t] = {b: np.asarray(v, dtype=float) for b, v in d[t].arrays(BRANCHES, library='np').items()}
    return {'pot': pot, 'trees': out}


def read_weights(path):
    """{(run, subrun, event): (weight, radiative)} from the sidecar weights CSV."""
    side = {}
    rows = 0
    for r in csv.DictReader(open(path)):
        rows += 1
        side[(int(r['run']), int(r['subrun']), int(r['event']))] = (float(r['weight']), int(r['radiative']))
    if len(side) != rows:
        raise SystemExit(f'{path}: sidecar keys not unique ({rows} rows, {len(side)} keys)')
    return side


def nc_res(d):
    return (d['true_cc'] == 0) & (d['true_interaction_mode'] == GENIE_RES)


def join(d, side):
    """Per-row weight (1 without a sidecar row), joined mask and radiative flag (-1 if unjoined)."""
    keys = zip(d['Run'].astype(int).tolist(), d['Subrun'].astype(int).tolist(), d['Evt'].astype(int).tolist())
    got = [side.get(k) for k in keys]
    joined = np.array([g is not None for g in got], dtype=bool)
    w = np.array([g[0] if g else 1.0 for g in got])
    rad = np.array([g[1] if g else -1 for g in got], dtype=int)
    return w, joined, rad


def neutrino(d):
    """Rows matched to a true neutrino. Cosmic rows have category 9; reco
    interactions with no truth match at all have NaN truth (category included)."""
    c = d['true_category']
    return np.isfinite(c) & (c != COSMIC)


def consistency(d, joined, rad):
    """Problems with the ncdeltarad rows' truth and the sidecar join (list of (level, message))."""
    out = []
    nu = neutrino(d)
    n_bad = int(np.sum(nu & ~nc_res(d)))
    if n_bad:
        out.append(('FAIL', f'{n_bad} neutrino rows are not NC RES (the sample is NCRES only)'))
    delta = nu & (d['true_baryon_res_code'] == 0)
    n = int(np.sum(delta & ~joined))
    if n:
        out.append(('FAIL', f'{n} Δ(1232) rows have no sidecar row: their weight would wrongly be 1'))
    n = int(np.sum(nu & joined & (d['true_baryon_res_code'] != 0)))
    if n:
        out.append(('WARN', f'{n} rows join the sidecar but are not Δ(1232) (resnum != 0)'))
    sig = nu & np.isin(d['true_category'], [0, 1])
    n = int(np.sum(sig & (rad == 0)))
    if n:
        out.append(('INFO', f'{n}/{int(sig.sum())} NC Δ→Nγ-category rows come from a non-radiative Δ decay '
                            f'(FSI or de-excitation photon)'))
    return out


def summarize(cv, nd, w, pot_cv, pot_nd):
    """Per-category comparison rows for one tree (see the module docstring)."""
    scale = pot_cv / pot_nd
    res_cv = nc_res(cv)
    rows = []
    for c, name in enumerate(CATS):
        in_cv = cv['true_category'] == c
        m = (nd['true_category'] == c) & (c != COSMIC)
        sw, sw2 = float(np.sum(w[m])), float(np.sum(w[m] ** 2))
        n_res = int(np.sum(in_cv & res_cv))
        pred, err = sw * scale, math.sqrt(sw2) * scale
        den = math.sqrt(err ** 2 + n_res)
        rows.append({
            'category': c, 'name': name,
            'cv': int(np.sum(in_cv)), 'cv_res': n_res, 'ncdr': int(np.sum(m)),
            'sum_w': sw, 'pred': pred, 'err': err,
            'pull': (pred - n_res) / den if den > 0 else 0.0,
            'gain': (pot_nd * sw / sw2) / pot_cv if sw2 > 0 else 0.0,
        })
    return rows


def closure(rows):
    """Overall ratio pred / cv_res (with its error) and chi2/ndf over categories with any entries."""
    used = [r for r in rows if r['cv_res'] > 0 or r['ncdr'] > 0]
    p, n = sum(r['pred'] for r in used), sum(r['cv_res'] for r in used)
    var = sum(r['err'] ** 2 for r in used)
    ratio = p / n if n else float('nan')
    rerr = ratio * math.sqrt(var / p ** 2 + 1 / n) if n and p else float('nan')
    chi2 = sum(r['pull'] ** 2 for r in used)
    return ratio, rerr, chi2, len(used)


def fmt_table(tree, rows):
    out = [f'\n== {tree}',
           f'{"category":<26}{"cv":>8}{"cv_res":>8}{"ncdr":>8}{"pred":>10}{"±err":>8}{"pull":>7}{"gain":>8}']
    for r in rows:
        if not (r['cv'] or r['ncdr']):
            continue
        g = f'{r["gain"]:.1f}x' if r['ncdr'] else '-'
        out.append(f'{r["name"]:<26}{r["cv"]:>8}{r["cv_res"]:>8}{r["ncdr"]:>8}{r["pred"]:>10.1f}{r["err"]:>8.1f}'
                   f'{r["pull"]:>7.2f}{g:>8}')
    ratio, rerr, chi2, ndf = closure(rows)
    out.append(f'closure: pred/cv_res = {ratio:.3f} ± {rerr:.3f}; chi2 {chi2:.1f} / {ndf} categories')
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ncdr', required=True, help='medulla output of the ncdeltarad sample')
    ap.add_argument('--cv', required=True, help='medulla output with the CV sample (same selection config)')
    ap.add_argument('--weights', required=True, help='gen_sidecar/ncdeltarad_weights.csv')
    ap.add_argument('--ncdr-key', default='ncdeltarad')
    ap.add_argument('--cv-key', default='cv')
    ap.add_argument('--json', help='also write the tables here')
    a = ap.parse_args(argv)

    nd, cv = load(a.ncdr, a.ncdr_key), load(a.cv, a.cv_key)
    side = read_weights(a.weights)
    print(f'POT: ncdeltarad {nd["pot"]:.4g}, cv {cv["pot"]:.4g} (ratio {nd["pot"] / cv["pot"]:.3f}); '
          f'sidecar {len(side)} NC Δ events')
    nfail = 0
    report = {'pot_ncdr': nd['pot'], 'pot_cv': cv['pot'], 'trees': {}}
    for t in TREES:
        if t not in nd['trees'] or t not in cv['trees']:
            continue
        d = nd['trees'][t]
        w, joined, rad = join(d, side)
        nu = neutrino(d)
        n_nan = int(np.sum(~np.isfinite(d['true_category'])))
        n_cv_nan = int(np.sum(~np.isfinite(cv['trees'][t]['true_category'])))
        print(f'\n{t}: ncdeltarad {len(w)} rows ({int(np.sum(~nu)) - n_nan} cosmic + {n_nan} without any truth '
              f'match, dropped; the CV has {n_cv_nan} without a match, not counted either); '
              f'{int(np.sum(joined & nu))}/{int(nu.sum())} neutrino rows join the sidecar')
        for level, msg in consistency(d, joined, rad):
            nfail += level == 'FAIL'
            print(f'{level:<5} {msg}')
        rows = summarize(cv['trees'][t], d, w, cv['pot'], nd['pot'])
        print('\n'.join(fmt_table(t, rows)))
        report['trees'][t] = {'rows': rows, 'closure': dict(zip(('ratio', 'err', 'chi2', 'ndf'), closure(rows)))}
    if a.json:
        with open(a.json, 'w') as f:
            json.dump(report, f, indent=1)
    print(f'\n{nfail} FAIL')
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
