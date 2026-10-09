#!/usr/bin/env python3
"""Summarize the stage-1 cut-flow diagnostic (selection/toml/gOre_osc_run2_cutflow_diag.toml).

Usage: osc_cutflow_diag.py MERGED_ROOT

Compares the new production (cv_new) with the old production's nominal-POT cv
MC (cv_old), normalized per 1e20 POT, to find which stage1 cut (fiducial,
containment, flash [0,1.6] us, 1g1p topology, single proton) loses the new MC's
neutrinos, and where its reco flash times sit compared with data.
"""
import sys

import numpy as np
import uproot

FLASH_LO, FLASH_HI = 0.0, 1.6
EDGES = np.arange(-10.0, 10.01, 0.5)


def pot(f, s):
    try:
        return float(f[f'events/{s}/POT'].values()[0])
    except Exception:
        return float('nan')


def tree(f, s, name):
    key = f'events/{s}/{name}'
    try:
        return f[key].arrays(library='np')
    except Exception:
        return None


def is1(a):
    return np.nan_to_num(a.astype(float), nan=0.0) == 1


def in_window(t):
    t = t.astype(float)
    return np.isfinite(t) & (t >= FLASH_LO) & (t <= FLASH_HI)


def flash_hist(label, t):
    t = t.astype(float)
    t = t[np.isfinite(t)]
    if len(t) == 0:
        print(f'    {label}: no finite flash times'); return
    h, _ = np.histogram(t, bins=EDGES)
    under, over = (t < EDGES[0]).sum(), (t > EDGES[-1]).sum()
    pct = np.percentile(t, [5, 25, 50, 75, 95])
    print(f'    {label}: n={len(t)} p5/25/50/75/95 = ' + ' / '.join(f'{p:.3g}' for p in pct)
          + f'  in [{FLASH_LO},{FLASH_HI}]: {in_window(t).mean():.1%}')
    print(f'      0.5us bins from {EDGES[0]:g} to {EDGES[-1]:g} us (under {under}, over {over}):')
    print('      ' + ' '.join(f'{x}' for x in h))


def true_cutflow(f, s):
    a = tree(f, s, 'diag_true_nu')
    if a is None:
        return
    norm = 1e20 / pot(f, s)
    print(f'\n=== {s}: diag_true_nu (true neutrinos in the true FV), POT {pot(f, s):.4g}')
    if 'event_trigger_within_gate' in a:
        tg = a['event_trigger_within_gate'].astype(float)
        print(f'  trigger_within_gate: min {np.nanmin(tg):.4g} median {np.nanmedian(tg):.4g} max {np.nanmax(tg):.4g}')
    cat = a['true_gOre_category'].astype(float)
    for sel_name, sel in (('all', np.ones(len(cat), bool)), ('gOre_category==0', cat == 0)):
        matched = np.isfinite(a['reco_fiducial'].astype(float))
        steps = [
            ('true nu in FV', sel),
            ('reco-matched', matched),
            ('reco fiducial', is1(a['reco_fiducial'])),
            ('reco contained', is1(a['reco_containment'])),
            ('valid flash match', is1(a['reco_flash'])),
            (f'flash in [{FLASH_LO},{FLASH_HI}]', in_window(a['reco_flash_time'])),
            ('1g topology', is1(a['reco_gOre_topology'])),
            ('1 proton', np.nan_to_num(a['reco_n_protons'].astype(float), nan=-1) == 1),
        ]
        print(f'  -- {sel_name}: sequential (per 1e20 POT) | N-1 pass rate among reco-matched')
        cum = np.ones(len(cat), bool)
        base = sel & matched
        for name, m in steps:
            cum &= m
            n1 = (m & base).sum() / max(base.sum(), 1)
            print(f'    {name:24s} {cum.sum() * norm:12.1f}   {n1:7.1%}')
    m = np.isfinite(a['reco_fiducial'].astype(float))
    flash_hist('reco flash_time of reco-matched true nu', a['reco_flash_time'][m & is1(a['reco_flash'])])
    vf = is1(a['reco_flash'])
    print(f'    valid flash match among reco-matched: {vf[m].mean():.1%}')


def reco_trees(f, s):
    p = pot(f, s)
    norm = (1e20 / p) if np.isfinite(p) and p > 0 else float('nan')
    a = tree(f, s, 'diag_reco_noflash')
    if a is not None:
        n = len(a['reco_flash'])
        vf = is1(a['reco_flash'])
        win = in_window(a['reco_flash_time'])
        print(f'\n=== {s}: diag_reco_noflash (stage1 minus flash), n={n} ({n * norm:.1f} per 1e20 POT)'
              f'  valid flash {vf.mean() if n else 0:.1%}  in window {win.mean() if n else 0:.1%}')
        flash_hist('reco flash_time', a['reco_flash_time'][vf])
    a = tree(f, s, 'diag_reco_noprong')
    if a is not None:
        n = len(a['reco_flash_time'])
        print(f'=== {s}: diag_reco_noprong (FV + contained + in-time flash), n={n} ({n * norm:.1f} per 1e20 POT)')
        if n:
            topo = is1(a['reco_gOre_topology'])
            nsh = np.nan_to_num(a['reco_n_primary_showers'].astype(float), nan=-1).astype(int)
            npr = np.nan_to_num(a['reco_n_protons'].astype(float), nan=-1).astype(int)
            print(f'    1g topology pass {topo.mean():.1%};  +1 proton {(topo & (npr == 1)).mean():.1%}')
            for lab, v in (('n_primary_showers', nsh), ('n_protons', npr)):
                u, c = np.unique(np.clip(v, -1, 5), return_counts=True)
                print(f'    {lab:18s} ' + '  '.join(f'{k}:{x / n:.1%}' for k, x in zip(u, c)))


def main():
    if len(sys.argv) != 2:
        print(__doc__); return 1
    f = uproot.open(sys.argv[1])
    samples = sorted({k.split(';')[0] for k in f['events'].keys(recursive=False)})
    print('samples:', samples)
    for s in samples:
        if pot(f, s) > 0 and tree(f, s, 'diag_true_nu') is not None:
            true_cutflow(f, s)
    for s in samples:
        reco_trees(f, s)
    return 0


if __name__ == '__main__':
    sys.exit(main())
