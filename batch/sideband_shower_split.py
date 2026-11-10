#!/usr/bin/env python3
"""
Shower-multiplicity split of a gOre sideband: do the pi0-rejection cuts tell
showers from the neutrino interaction apart from secondary (pi0-decay,
reinteraction) showers where nothing is ever selected?

The non-fiducial NgNp sideband (outside_signal_box + >=1 shower + >=1 proton,
gOre_1g1p_sidebands_osc_run2.toml) is disjoint from the 1g1p selection by
construction, so it can be looked at freely. Inside it:
  signal-like   reco_n_primary_showers == 1   (looks like a lone photon)
  pi0-like      reco_n_primary_showers == 2   (looks like a pi0 decay)
Each cut of the pi0 rejection (and then e/gamma separation) is applied to both
classes, at the thresholds given (default: the 2D re-tune), and the data and MC
pass fractions compared.

MC truth is grouped by the origin of the true showers (photons/electrons above
the shower threshold) of the matched interaction. A shower is "from the
neutrino interaction" when the particle is its own ancestor (a generator
final-state particle, CC or NC); it is "secondary" when it descends from another
final-state particle (pi0 decays, hadron reinteractions, delta rays). Branches
true_n_nu_photons, true_n_nu_electrons, true_n_secondary_true_showers:
  true single gamma primary      1 photon from the nu, no other true shower
  true single e primary          1 electron from the nu, no other true shower
  true multiple primary e/gamma  >=2 showers from the nu, or 1 plus secondaries
  true secondary showers         only secondary showers (pi0-like)
  other nu                       no true shower above threshold
  cosmics                        off-beam data, scaled by livetime
  dirt                           out-of-FV neutrinos (dirt MC, matched rows)
The cuts validate the pi0 rejection if they keep a much larger fraction of the
true-single-gamma (and single-e) groups than of the true-secondary group, and
the MC models them if its pass fractions match the data's in both classes.

Input: a selection-only medulla output (events/cv, offbeam, onbeam[, dirt]) or a
sys ROOT (events/full, offbeam, onbeam). The cv rows without a neutrino match
(overlay cosmics; category 9/NaN) are dropped, as run_systematics drops them
from events/full: the off-beam sample is the only cosmic source.

Usage
  sideband_shower_split.py build/output_gOre_1g1p_sidebands_osc_run2_origin.root --out build/sideband_split \\
      [--dirt <file with events/dirt, default: the input if it has one>] [--gps 0.9985 --dm-halfwidth 95]
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
ORIGIN = ['true_n_nu_photons', 'true_n_nu_electrons', 'true_n_secondary_true_showers']
GROUPS = [('True single $\\gamma$ primary', '#E69900'),
          ('True single $e$ primary', '#2CA02C'),
          ('True multiple primary $e/\\gamma$', '#1F77B4'),
          ('True secondary showers', '#F2E640'),
          ('Other $\\nu$', '#9980CC'),
          ('Cosmics', '#808080'),
          ('Dirt', '#8C564B')]
SINGLE_G, SINGLE_E, MULTI, SECONDARY, OTHER, COSMIC, DIRT = range(7)
CLASSES = [('signal-like (1 shower)', 1), ('pi0-like (2 showers)', 2)]
PRELIM = 'ICARUS Analysis In-Progress'
N_PI0_REJ = 3  # the first three cuts are the pi0 rejection


def plain(label):
    for tex, txt in (('$e/\\gamma$', 'e/gamma'), ('$\\gamma$', 'gamma'), ('$e$', 'e'), ('$\\nu$', 'nu')):
        label = label.replace(tex, txt)
    return label


def group_of(d):
    """Truth group per row from the shower-origin counts (NaN counts: no match -> other nu)."""
    ng, ne, ns = (np.nan_to_num(d[c].to_numpy(dtype=float), nan=0.0) for c in ORIGIN)
    nnu = ng + ne
    g = np.full(len(d), OTHER)
    g[(nnu == 0) & (ns > 0)] = SECONDARY
    g[(nnu >= 2) | ((nnu == 1) & (ns > 0))] = MULTI
    g[(ng == 1) & (ne == 0) & (ns == 0)] = SINGLE_G
    g[(ne == 1) & (ng == 0) & (ns == 0)] = SINGLE_E
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


def has(path, key):
    import uproot
    f = uproot.open(path)
    return 'events' in f and key in f['events']


def load(path, key, want):
    import uproot
    f = uproot.open(path)
    if 'events' not in f or key not in f['events']:
        raise SystemExit(f'{path}: no events/{key}')
    t = f[f'events/{key}/{TREE}']
    miss = [b for b in want if b not in t.keys()]
    if miss:
        hint = (' (rerun the sidebands selection with the shower-origin branches: see the module docstring)'
                if set(miss) & set(ORIGIN) else '')
        raise SystemExit(f'{path}: events/{key}/{TREE} lacks {miss}{hint}')
    d = f[f'events/{key}']
    expo = [float(np.sum(d[h].values())) if h in d else None for h in ('POT', 'Livetime')]
    return (t.arrays(want, library='pd'), *expo)


def exposure(val, path, key, hist):
    """`val` if the sample carried it, else events/<key>/<hist> from the --exposure-from file."""
    if val is not None:
        return val
    if not path:
        raise SystemExit(f'events/{key} has no {hist} histogram (sys ROOT made before the off-beam exposure fix?): '
                         f'pass --exposure-from <a ROOT with the same onbeam/offbeam samples>, e.g. the signal one')
    import uproot
    f = uproot.open(path)
    v = float(np.sum(f[f'events/{key}/{hist}'].values()))
    print(f'  {key} {hist} {v:.4g} taken from {path}')
    return v


def matched(d):
    """Rows with a neutrino match (drops overlay cosmics: category 9 or NaN)."""
    c = d['true_category'].to_numpy(dtype=float)
    return d[np.isfinite(c) & (c != 9)].copy()


def frac(w, m):
    """Weighted pass fraction and its error (binomial, effective entries)."""
    sw = w.sum()
    if sw <= 0:
        return float('nan'), float('nan')
    p = w[m].sum() / sw
    neff = sw ** 2 / (w ** 2).sum()
    return p, np.sqrt(max(p * (1 - p), 0) / neff)


def poisson_z(n, mu):
    """Signed significance of n observed vs mu expected (Poisson likelihood ratio)."""
    if mu <= 0:
        return float('nan')
    t = 2 * (mu - n + (n * np.log(n / mu) if n > 0 else 0.0))
    return float(np.sign(n - mu) * np.sqrt(max(t, 0.0)))


def build(a, want):
    import pandas as pd
    key = 'full' if has(a.sys, 'full') else 'cv'
    full, pot_cv, _ = load(a.sys, key, want)
    off, _, lt_off = load(a.sys, 'offbeam', want)
    on, pot_on, lt_on = load(a.sys, 'onbeam', want)
    pot_cv = exposure(pot_cv, a.exposure_from, key, 'POT')
    lt_off = exposure(lt_off, a.exposure_from, 'offbeam', 'Livetime')
    pot_on = exposure(pot_on, a.exposure_from, 'onbeam', 'POT')
    lt_on = exposure(lt_on, a.exposure_from, 'onbeam', 'Livetime')
    if not (lt_off > 0 and lt_on > 0):
        raise SystemExit(f'off-beam/on-beam livetime {lt_off}/{lt_on}: cannot scale the off-beam cosmics')
    n0 = len(full)
    full = matched(full)
    print(f'events/{key}: {len(full)} of {n0} rows have a neutrino match (the rest are overlay cosmics, '
          f'covered by off-beam); POT {pot_cv:.4g}; off-beam livetime {lt_off:.4g}, on-beam {lt_on:.4g}')
    full['g'], full['w'] = group_of(full), pot_on / pot_cv
    off['g'], off['w'] = COSMIC, lt_on / lt_off
    parts = [full, off]
    dirt = a.dirt or (a.sys if has(a.sys, 'dirt') else None)
    if dirt:
        dr, pot_dr, _ = load(dirt, 'dirt', want)
        dr = matched(dr)
        dr['g'], dr['w'] = DIRT, pot_on / pot_dr
        parts.append(dr)
        print(f'dirt from {dirt}: {len(dr)} matched rows, POT {pot_dr:.4g}')
    else:
        print('WARNING: no dirt: the out-of-FV neutrino background is missing (it matters most here)')
    return pd.concat(parts, ignore_index=True), on, pot_on


def main(argv=None):
    import pandas as pd
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('sys', help='sidebands medulla output: events/{cv|full,offbeam,onbeam[,dirt]}')
    ap.add_argument('--dirt', help='file with events/dirt (default: the input file, if it has one)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--exposure-from', help='ROOT file to take missing POT/Livetime histograms from (same onbeam/'
                                             'offbeam samples; only needed for old sys ROOTs)')
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
            'reco_leading_primary_gOre_photon_softmax', 'reco_leading_primary_gOre_shower_ke'] + ORIGIN
    mc, on, pot_on = build(a, want)
    print(f'{TREE}: onbeam {pot_on:.3g} POT; MC {mc["w"].sum():.1f} expected, data {len(on)}')

    cl = cuts(thr)
    rows = []
    for cname, n in CLASSES:
        sub, dsub = mc[mc[NSH].to_numpy() == n], on[on[NSH].to_numpy() == n]
        w, g = sub['w'].to_numpy(), sub['g'].to_numpy()
        tot = w.sum()
        print(f'\n== {cname}: MC {tot:.1f}, data {len(dsub)}, data/MC {len(dsub) / tot if tot else float("nan"):.2f}')
        for i, (label, _) in enumerate(GROUPS):
            v, raw = w[g == i].sum(), int((g == i).sum())
            if raw:
                print(f'    {plain(label):<32} {v:8.1f} ({v / tot:5.1%}; {raw} MC rows)')
        pass_mc = np.ones(len(sub), bool)
        pass_d = np.ones(len(dsub), bool)
        print(f'    {"cumulative cut":<40} {"MC exp":>7} {"data":>5} {"z":>5} {"MC pass":>14} {"data pass":>15} '
              f'{"1g":>6} {"1e":>6} {"sec":>6}')
        for name, fn in cl:
            pass_mc &= fn(sub)
            pass_d &= fn(dsub)
            pm, em = frac(w, pass_mc)
            k, nd = int(pass_d.sum()), len(dsub)
            pd_, ed = (k / nd, np.sqrt(k / nd * (1 - k / nd) / nd)) if nd else (float('nan'), float('nan'))
            mu = w[pass_mc].sum()
            p_g = {i: frac(w[g == i], pass_mc[g == i])[0] for i in (SINGLE_G, SINGLE_E, SECONDARY)}
            print(f'    {name:<40} {mu:7.1f} {k:5d} {poisson_z(k, mu):5.1f} {pm:7.4f}±{em:.4f} {pd_:7.4f}±{ed:.4f} '
                  f'{p_g[SINGLE_G]:6.3f} {p_g[SINGLE_E]:6.3f} {p_g[SECONDARY]:6.3f}')
            rows.append(dict(cls=cname, cut=name, mc_expected=mu, data_n=k, z=poisson_z(k, mu), mc=pm, mc_err=em,
                             data=pd_, data_err=ed, single_gamma=p_g[SINGLE_G], single_e=p_g[SINGLE_E],
                             secondary=p_g[SECONDARY]))
    pd.DataFrame(rows).to_csv(os.path.join(a.out, 'pass_fractions.csv'), index=False)
    separation(mc, cl)
    plots(mc, on, thr, a.out, pot_on)
    print(f'\nwrote {a.out}/pass_fractions.csv and plots')
    return 0


def separation(mc, cl):
    """True-origin survival through the pi0 rejection (and all cuts), both reco classes together (MC)."""
    sub = mc[mc[NSH].isin([n for _, n in CLASSES])]
    w, g = sub['w'].to_numpy(), sub['g'].to_numpy()
    print('\nseparation by truth (MC, 1- and 2-shower classes together): fraction of each truth group kept')
    for upto, label in ((N_PI0_REJ, 'pi0 rejection (first 3 cuts)'), (len(cl), 'pi0 rejection + e/gamma')):
        m = np.ones(len(sub), bool)
        for _, fn in cl[:upto]:
            m &= fn(sub)
        f = {i: frac(w[g == i], m[g == i]) for i in (SINGLE_G, SINGLE_E, MULTI, SECONDARY)}
        raw = {i: int(((g == i) & m).sum()) for i in f}
        sec = f[SECONDARY][0]
        print(f'  {label}:')
        for i, (p, e) in f.items():
            if not np.isfinite(p):
                print(f'    {plain(GROUPS[i][0]):<32} no MC rows')
                continue
            r = f'  x{p / sec:.1f} the secondary-shower rate' if sec and i != SECONDARY else ''
            print(f'    {plain(GROUPS[i][0]):<32} {p:.4f} ± {e:.4f} ({raw[i]} MC rows pass){r}')


def plots(mc, on, thr, out, pot_on):
    try:
        import matplotlib
    except ImportError:
        print('no matplotlib: plots skipped')
        return
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.transforms as mtransforms
    here = os.path.dirname(os.path.abspath(__file__))
    sheet = os.path.join(here, '..', 'spineplot', 'styles', 'muon2024.mplstyle')
    if os.path.exists(sheet):
        plt.style.use(sheet)
    # (variable, range, bins, label, [(cut value, side kept: +1 above / -1 below)])
    variables = [(GPS, (0.99, 1.0), 50, 'Shower Primary Softmax', [(thr['gps'], +1)]),
                 (DM, (950.0, 1600.0), 26, r'Reconstructed $M_\Delta$ [MeV/$c^2$]',
                  [(1232.0 - thr['dm'], +1), (1232.0 + thr['dm'], -1)]),
                 (PPS, (0.0, 1.0), 50, 'Leading Proton Primary Softmax', [(thr['pps'], +1)]),
                 ('reco_leading_primary_gOre_shower_ke', (0.0, 1000.0), 25, 'Shower KE [MeV]', [])]
    for var, (lo, hi), nb, xl, cutlines in variables:
        fig, axs = plt.subplots(1, 2, figsize=(16, 6.5))
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
            ax.set_xlim(lo, hi)
            ax.set_ylim(0, 1.45 * max(ax.get_ylim()[1], 1.0))
            trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
            span = 0.06 * (hi - lo)
            for j, (x, side) in enumerate(cutlines):
                ax.axvline(x, color='r', ls='--', label='cut (arrow: kept)' if j == 0 else None)
                ax.annotate('', xy=(x + side * span, 0.6), xytext=(x, 0.6), xycoords=trans, textcoords=trans,
                            arrowprops=dict(arrowstyle='->', color='r', lw=2))
            ax.set_xlabel(xl)
            ax.set_ylabel('Candidates')
            ax.text(0.0, 1.01, PRELIM, transform=ax.transAxes, ha='left', va='bottom', color='red',
                    fontweight='bold', fontsize='medium')
            ax.text(1.0, 1.01, f'{pot_on:.2e} POT', transform=ax.transAxes, ha='right', va='bottom',
                    fontsize='medium')
            ax.legend(loc='upper right', fontsize='x-small', frameon=False,
                      title=f'non-fid. sideband, {cname}', title_fontsize='small')
        fig.tight_layout()
        for ext in ('png', 'pdf'):
            fig.savefig(os.path.join(out, f'{var}.{ext}'))
        plt.close(fig)


if __name__ == '__main__':
    sys.exit(main())
