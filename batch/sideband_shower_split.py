#!/usr/bin/env python3
"""
Shower-multiplicity split of a gOre sideband: can the pi0-rejection cuts tell
single-photon-like events from pi0-like ones where nothing is ever selected?

The non-fiducial NgNp sideband (outside_signal_box + >=1 shower + >=1 proton,
gOre_1g1p_sidebands_osc_run2.toml) is disjoint from the 1g1p selection by
construction, so it can be looked at freely. Inside it:
  signal-like   reco_n_primary_showers == 1   (looks like a lone photon)
  pi0-like      reco_n_primary_showers == 2   (looks like a pi0 decay)
Each cut of the pi0 rejection (and then e/gamma separation) is applied to both
classes, at the thresholds given (default: the 2D re-tune), and the pass
fractions compared, data vs MC. The cuts separate the classes if they keep a
much larger fraction of the signal-like class; the MC models that separation
if its pass fractions match the data's.

MC truth is grouped as true 1gamma (categories 0-2), true NC pi0 (3, 4), other
neutrino (5-8), cosmics (off-beam, scaled by livetime) and dirt (--dirt). The
sys ROOT's events/full holds only truth-matched neutrinos (run_systematics
routes overlay cosmics to *_nonmatched), so the off-beam sample is the only
cosmic source; dirt rows without a neutrino match are dropped likewise.

Usage
  sideband_shower_split.py build/output_gOre_1g1p_sidebands_osc_run2_sys.root --out build/sideband_split \\
      [--dirt build/gOre_1g1p_sidebands_osc_run2_dirt.root] [--gps 0.9985 --dm-halfwidth 95]
"""
import argparse
import os
import sys

import numpy as np

TREE = 'selected_nonfid_NgNp_stage1'
GPS = 'reco_leading_primary_gOre_primary_softmax'
PPS = 'reco_leading_primary_proton_primary_softmax'
DM = 'reco_delta_mass_P'
NSH = 'reco_n_primary_showers'
GROUPS = [('true 1$\\gamma$ (NC $\\Delta{\\to}N\\gamma$, other NC 1$\\gamma$)', '#E69900'),
          ('true NC $\\pi^0$', '#F2E640'), ('other $\\nu$', '#9980CC'),
          ('cosmic (off-beam)', '#808080'), ('dirt (out-of-FV $\\nu$)', '#8C564B')]
COSMIC, DIRT = 3, 4
CLASSES = [('signal-like (1 shower)', 1), ('pi0-like (2 showers)', 2)]


def group_of(cat):
    c = np.asarray(cat, dtype=float)
    g = np.full(len(c), 2)
    g[np.isin(c, [0, 1, 2])] = 0
    g[np.isin(c, [3, 4])] = 1
    g[~np.isfinite(c) | (c == 9)] = COSMIC
    return g


def cuts(thr):
    """(name, function of a dataframe -> bool array); cumulative in this order."""
    def gt(col, t):
        return lambda d: np.nan_to_num(d[col].to_numpy(dtype=float), nan=-np.inf) > t

    def lt(col, t):
        return lambda d: np.nan_to_num(d[col].to_numpy(dtype=float), nan=np.inf) < t
    dm = lambda d: np.nan_to_num(np.abs(d[DM].to_numpy(dtype=float) - 1232.0), nan=np.inf) < thr['dm']
    return [(f'gamma primary softmax > {thr["gps"]}', gt(GPS, thr['gps'])),
            (f'proton primary softmax > {thr["pps"]}', gt(PPS, thr['pps'])),
            (f'|M_Delta - 1232| < {thr["dm"]} MeV', dm),
            (f'start dE/dx > {thr["dedx"]}', gt('reco_leading_primary_gOre_start_dedx', thr['dedx'])),
            (f'axial spread > {thr["axial"]}', gt('reco_leading_primary_gOre_axial_spread', thr['axial'])),
            (f'directional spread < {thr["dirs"]}', lt('reco_leading_primary_gOre_directional_spread', thr['dirs'])),
            (f'photon softmax > {thr["phs"]}', gt('reco_leading_primary_gOre_photon_softmax', thr['phs']))]


def load(path, key, want):
    import uproot
    f = uproot.open(path)
    if 'events' not in f or key not in f['events']:
        raise SystemExit(f'{path}: no events/{key}')
    t = f[f'events/{key}/{TREE}']
    miss = [b for b in want if b not in t.keys()]
    if miss:
        raise SystemExit(f'{path}: events/{key}/{TREE} lacks {miss}')
    d = f[f'events/{key}']
    return (t.arrays(want, library='pd'), float(np.sum(d['POT'].values())), float(np.sum(d['Livetime'].values())))


def frac(w, m):
    """Weighted pass fraction and its error (binomial, effective entries)."""
    sw = w.sum()
    if sw <= 0:
        return float('nan'), float('nan')
    p = w[m].sum() / sw
    neff = sw ** 2 / (w ** 2).sum()
    return p, np.sqrt(max(p * (1 - p), 0) / neff)


def main(argv=None):
    import pandas as pd
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('sys', help='sidebands sys ROOT (events/{full,offbeam,onbeam})')
    ap.add_argument('--dirt', help='medulla output with the dirt MC run through the sidebands config (events/dirt)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--gps', type=float, default=0.9985, help='gamma primary softmax (2D re-tune; nominal 0.9991)')
    ap.add_argument('--dm-halfwidth', type=float, default=95.0, help='Delta mass halfwidth [MeV] (2D re-tune; nominal 59.7)')
    ap.add_argument('--pps', type=float, default=0.9513)
    ap.add_argument('--dedx', type=float, default=4.07)
    ap.add_argument('--axial', type=float, default=-0.5)
    ap.add_argument('--dirs', type=float, default=0.037)
    ap.add_argument('--phs', type=float, default=0.9225)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    thr = dict(gps=a.gps, pps=a.pps, dm=a.dm_halfwidth, dedx=a.dedx, axial=a.axial, dirs=a.dirs, phs=a.phs)
    want = [GPS, PPS, DM, NSH, 'true_category', 'reco_leading_primary_gOre_start_dedx',
            'reco_leading_primary_gOre_axial_spread', 'reco_leading_primary_gOre_directional_spread',
            'reco_leading_primary_gOre_photon_softmax', 'reco_leading_primary_gOre_shower_ke']

    full, pot_cv, _ = load(a.sys, 'full', want)
    off, _, lt_off = load(a.sys, 'offbeam', want)
    on, pot_on, lt_on = load(a.sys, 'onbeam', want)
    full['g'], full['w'] = group_of(full['true_category']), pot_on / pot_cv
    off['g'], off['w'] = COSMIC, lt_on / lt_off
    parts = [full, off]
    if a.dirt:
        dr, pot_dr, _ = load(a.dirt, 'dirt', want)
        c = dr['true_category'].to_numpy(dtype=float)
        dr = dr[np.isfinite(c) & (c != 9)].copy()
        dr['g'], dr['w'] = DIRT, pot_on / pot_dr
        parts.append(dr)
    else:
        print('WARNING: no --dirt: the out-of-FV neutrino background is missing (it matters most here)')
    mc = pd.concat(parts, ignore_index=True)
    print(f'{TREE}: onbeam {pot_on:.3g} POT; MC {mc["w"].sum():.1f} expected, data {len(on)}')

    cl = cuts(thr)
    rows = []
    for cname, n in CLASSES:
        m_mc = mc[NSH].to_numpy() == n
        m_d = on[NSH].to_numpy() == n
        sub, dsub = mc[m_mc], on[m_d]
        w = sub['w'].to_numpy()
        tot = w.sum()
        comp = {GROUPS[g][0]: sub.loc[sub['g'] == g, 'w'].sum() for g in range(len(GROUPS))}
        print(f'\n== {cname}: MC {tot:.1f}, data {len(dsub)}, data/MC {len(dsub) / tot if tot else float("nan"):.2f}')
        for k, v in comp.items():
            if v:
                print(f'    {k:<58} {v:8.1f} ({v / tot:.1%})')
        pass_mc = np.ones(len(sub), bool)
        pass_d = np.ones(len(dsub), bool)
        print(f'    {"cumulative cut":<44} {"MC pass":>14} {"data pass":>16} {"true 1g":>8} {"true pi0":>9}')
        for name, fn in cl:
            pass_mc &= fn(sub)
            pass_d &= fn(dsub)
            pm, em = frac(w, pass_mc)
            k, nd = int(pass_d.sum()), len(dsub)
            pd_, ed = (k / nd, np.sqrt(k / nd * (1 - k / nd) / nd)) if nd else (float('nan'), float('nan'))
            g = sub['g'].to_numpy()
            p1 = frac(w[g == 0], pass_mc[g == 0])[0]
            p2 = frac(w[g == 1], pass_mc[g == 1])[0]
            print(f'    {name:<44} {pm:7.3f}±{em:.3f} {pd_:8.3f}±{ed:.3f} {p1:8.3f} {p2:9.3f}')
            rows.append(dict(cls=cname, cut=name, mc=pm, mc_err=em, data=pd_, data_err=ed, data_n=k,
                             true_1g=p1, true_pi0=p2))
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(a.out, 'pass_fractions.csv'), index=False)
    print('\nseparation (cumulative, after the pi0 rejection = first 3 cuts, and after all):')
    for i, label in ((2, 'pi0 rejection'), (len(cl) - 1, 'pi0 rejection + e/gamma')):
        s1 = tab[(tab['cls'] == CLASSES[0][0])].iloc[i]
        s2 = tab[(tab['cls'] == CLASSES[1][0])].iloc[i]
        print(f'  {label:<26} keeps signal-like {s1["mc"]:.3f} (data {s1["data"]:.3f}) vs pi0-like {s2["mc"]:.3f} '
              f'(data {s2["data"]:.3f}): ratio MC {s1["mc"] / s2["mc"] if s2["mc"] else float("inf"):.1f}, '
              f'data {s1["data"] / s2["data"] if s2["data"] else float("inf"):.1f}')
    plots(mc, on, thr, a.out, pot_on)
    print(f'\nwrote {a.out}/pass_fractions.csv and plots')
    return 0


def plots(mc, on, thr, out, pot_on):
    try:
        import matplotlib
    except ImportError:
        print('no matplotlib: plots skipped')
        return
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    here = os.path.dirname(os.path.abspath(__file__))
    sheet = os.path.join(here, '..', 'spineplot', 'styles', 'muon2024.mplstyle')
    if os.path.exists(sheet):
        plt.style.use(sheet)
    variables = [(GPS, (0.99, 1.0), 50, 'Shower Primary Softmax', thr['gps']),
                 (DM, (950.0, 1600.0), 26, r'Reconstructed $M_\Delta$ [MeV/$c^2$]', None),
                 (PPS, (0.0, 1.0), 50, 'Leading Proton Primary Softmax', thr['pps']),
                 ('reco_leading_primary_gOre_shower_ke', (0.0, 1000.0), 25, 'Shower KE [MeV]', None)]
    for var, (lo, hi), nb, xl, cut in variables:
        fig, axs = plt.subplots(1, 2, figsize=(16, 6), sharey=False)
        edges = np.linspace(lo, hi, nb + 1)
        for ax, (cname, n) in zip(axs, CLASSES):
            sub, dsub = mc[mc[NSH] == n], on[on[NSH] == n]
            xs = [sub.loc[sub['g'] == g, var].to_numpy(dtype=float) for g in range(len(GROUPS))]
            ws = [sub.loc[sub['g'] == g, 'w'].to_numpy(dtype=float) for g in range(len(GROUPS))]
            keep = [i for i in range(len(GROUPS)) if ws[i].sum() > 0]
            ax.hist([xs[i] for i in keep], bins=edges, weights=[ws[i] for i in keep], stacked=True,
                    histtype='barstacked', color=[GROUPS[i][1] for i in keep],
                    label=[f'{GROUPS[i][0]} ({ws[i].sum():.1f})' for i in keep])
            k = np.histogram(dsub[var].to_numpy(dtype=float), edges)[0]
            ctr = 0.5 * (edges[1:] + edges[:-1])
            ax.errorbar(ctr, k, yerr=np.sqrt(k), fmt='o', color='k', ms=4, label=f'Data (BNB) ({k.sum()})')
            if cut is not None:
                ax.axvline(cut, color='r', ls='--', label='re-tuned cut')
            if var == DM:
                for x in (1232.0 - thr['dm'], 1232.0 + thr['dm']):
                    ax.axvline(x, color='r', ls='--')
            ax.set_xlim(lo, hi)
            ax.set_xlabel(xl)
            ax.set_ylabel('Candidates')
            ax.set_title(f'non-fid. sideband, {cname}', fontsize='medium')
            ax.text(0.02, 0.97, f'{pot_on:.2e} POT', transform=ax.transAxes, va='top', fontsize='small')
            ax.legend(loc='upper right', fontsize='x-small', frameon=False)
            top = max(ax.get_ylim()[1], 1.0)
            ax.set_ylim(0, 1.4 * top)
        fig.tight_layout()
        for ext in ('png', 'pdf'):
            fig.savefig(os.path.join(out, f'{var}.{ext}'))
        plt.close(fig)


if __name__ == '__main__':
    sys.exit(main())
