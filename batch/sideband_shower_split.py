#!/usr/bin/env python3
"""
Validate the 1g1p cuts in the non-fiducial sideband.

A sideband is a selection disjoint from the signal selection by construction:
here the non-fiducial NgNp preselection (outside_signal_box + >=1 shower +
>=1 proton, gOre_1g1p_sidebands_osc_run2.toml) restricted to the signal's
topology, exactly one reco primary shower. Nothing in it can be selected by the
analysis, so it can be looked at freely, and the 1g1p cuts are applied to it as
they are to the signal region. They are validated if they keep the "signal-like"
events (a shower from the neutrino interaction) while rejecting the "pi0-like"
ones (showers from pi0 decays/reinteractions), and if the data follow the MC
through the cuts.

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

Output (--out):
  cut flow      expected events per truth group, data, Poisson z, and the
                fraction of true single gamma vs secondary showers kept
  N-1 table     with all other cuts applied, the fraction of each group the
                left-out cut rejects (its own contribution), MC and data
  before/       each cut variable at the preselection
  after/        each cut variable after all cuts
  nm1/          each cut variable with all the other cuts applied
  cutflow.csv

Input: a selection-only medulla output (events/cv, offbeam, onbeam[, dirt]) or a
sys ROOT (events/full, offbeam, onbeam). The cv rows without a neutrino match
(overlay cosmics; category 9/NaN) are dropped, as run_systematics drops them
from events/full: the off-beam sample is the only cosmic source
(--overlay-cosmics keeps them as their own category).

Usage
  sideband_shower_split.py build/output_gOre_1g1p_sidebands_osc_run2_origin.root --out build/sideband_cuts \\
      [--gps 0.9985 --dm-halfwidth 95] [--overlay-cosmics]
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
          ('Dirt', '#8C564B'),
          ('Cosmics in $\\nu$ events (MC overlay)', '#4D4D4D')]
SINGLE_G, SINGLE_E, MULTI, SECONDARY, OTHER, COSMIC, DIRT, OVERLAY = range(8)
SHOWN = list(range(len(GROUPS)))
SHORT = ['1g', '1e', 'multi', 'sec', 'oth nu', 'cosmic', 'dirt', 'ovl cos']
KE = 'reco_leading_primary_gOre_shower_ke'
PRELIM = 'ICARUS Analysis In-Progress'


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
    """The 1g1p cuts in order (pi0 rejection, then e/gamma separation): name, variable, pass function,
    cut lines [(value, kept side +1 above / -1 below)] and the plot binning (spineplot's, finer where a cut needs it)."""
    def gt(col, t):
        return lambda d: np.nan_to_num(d[col].to_numpy(dtype=float), nan=-np.inf) > t

    def lt(col, t):
        return lambda d: np.nan_to_num(d[col].to_numpy(dtype=float), nan=np.inf) < t
    dm = lambda d: np.nan_to_num(np.abs(d[DM].to_numpy(dtype=float) - 1232.0), nan=np.inf) < thr['dm']
    c = lambda name, var, fn, lines, rng, nb, xl: dict(name=name, var=var, fn=fn, lines=lines, range=rng,
                                                      nbins=nb, xlabel=xl)
    return [c(f'gamma primary softmax > {thr["gps"]}', GPS, gt(GPS, thr['gps']), [(thr['gps'], +1)],
              (0.99, 1.0), 20, 'Shower Primary Softmax'),
            c(f'proton primary softmax > {thr["pps"]}', PPS, gt(PPS, thr['pps']), [(thr['pps'], +1)],
              (0.0, 1.0), 50, 'Leading Proton Primary Softmax'),
            c(f'|M_Delta - 1232| < {thr["dm"]} MeV', DM, dm, [(1232.0 - thr['dm'], +1), (1232.0 + thr['dm'], -1)],
              (950.0, 1600.0), 26, r'Reconstructed $M_\Delta$ Peak [MeV/$c^2$]'),
            c(f'start dE/dx > {thr["dedx"]}', 'reco_leading_primary_gOre_start_dedx',
              gt('reco_leading_primary_gOre_start_dedx', thr['dedx']), [(thr['dedx'], +1)],
              (0.0, 25.0), 50, 'Shower Start dE/dx [MeV/cm]'),
            c(f'axial spread > {thr["axial"]}', 'reco_leading_primary_gOre_axial_spread',
              gt('reco_leading_primary_gOre_axial_spread', thr['axial']), [(thr['axial'], +1)],
              (-1.0, 1.0), 40, 'Shower Axial Spread'),
            c(f'directional spread < {thr["dirs"]}', 'reco_leading_primary_gOre_directional_spread',
              lt('reco_leading_primary_gOre_directional_spread', thr['dirs']), [(thr['dirs'], -1)],
              (0.0, 0.2), 40, 'Shower Directional Spread'),
            c(f'photon softmax > {thr["phs"]}', 'reco_leading_primary_gOre_photon_softmax',
              gt('reco_leading_primary_gOre_photon_softmax', thr['phs']), [(thr['phs'], +1)],
              (0.0, 1.0), 50, 'Shower Photon Softmax')]


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
    keep = matched(full)
    overlay = full.drop(keep.index)
    full = keep
    note = ('kept as "cosmics in nu events" (--overlay-cosmics)' if a.overlay_cosmics else
            'dropped: off-beam is the only cosmic source; --overlay-cosmics keeps them')
    print(f'events/{key}: {len(full)} of {n0} rows have a neutrino match; the other {len(overlay)} '
          f'(overlay cosmics, {len(overlay) * pot_on / pot_cv:.1f} at on-beam POT) {note}; '
          f'POT {pot_cv:.4g}; off-beam livetime {lt_off:.4g}, on-beam {lt_on:.4g}')
    full['g'], full['w'] = group_of(full), pot_on / pot_cv
    off['g'], off['w'] = COSMIC, lt_on / lt_off
    parts = [full, off]
    if a.overlay_cosmics:
        if key == 'full':
            raise SystemExit('--overlay-cosmics needs a selection output (events/cv): a sys ROOT keeps them '
                             'in *_nonmatched')
        overlay = overlay.copy()
        overlay['g'], overlay['w'] = OVERLAY, pot_on / pot_cv
        parts.append(overlay)
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
    ap.add_argument('--n-showers', type=int, default=1,
                    help='reco primary showers defining the sideband region (default 1, as in the 1g1p selection)')
    ap.add_argument('--overlay-cosmics', action='store_true',
                    help='also keep the cv rows without a neutrino match (cosmics recorded in neutrino-triggered '
                         'readouts) as their own category, POT-scaled; off-beam covers only beam-off triggers')
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
    cl = cuts(thr)
    want = sorted({c['var'] for c in cl} | {NSH, 'true_category', KE} | set(ORIGIN))
    mc, on, pot_on = build(a, want)
    region = f'non-fid. sideband, {a.n_showers} reco shower' + ('s' if a.n_showers != 1 else '')
    mc = mc[mc[NSH].to_numpy() == a.n_showers].reset_index(drop=True)
    on = on[on[NSH].to_numpy() == a.n_showers].reset_index(drop=True)
    print(f'{TREE}, {region}: onbeam {pot_on:.3g} POT; MC {mc["w"].sum():.1f} expected, data {len(on)}, '
          f'data/MC {len(on) / mc["w"].sum():.2f}')

    p_mc = np.column_stack([c['fn'](mc) for c in cl])
    p_on = np.column_stack([c['fn'](on) for c in cl])
    rows = []

    def report(stage, m_mc, m_on):
        w, g = mc['w'].to_numpy(), mc['g'].to_numpy()
        mu, k = w[m_mc].sum(), int(m_on.sum())
        row = dict(stage=stage, mc_expected=mu, data_n=k, z=poisson_z(k, mu))
        for i in SHOWN:
            row[plain(GROUPS[i][0])] = w[m_mc & (g == i)].sum()
        rows.append(row)
        return row

    print('\ncut flow (cumulative): expected events per truth group, and the fraction of each truth group kept')
    hdr = ''.join(f'{SHORT[i]:>9}' for i in SHOWN)
    print(f'  {"stage":<38} {"MC":>8} {"data":>6} {"z":>5} |{hdr} | {"1g kept":>8} {"sec kept":>8} {"1g/sec":>7}')
    w, g = mc['w'].to_numpy(), mc['g'].to_numpy()
    tot = {i: w[g == i].sum() for i in SHOWN}
    m_mc, m_on = np.ones(len(mc), bool), np.ones(len(on), bool)
    stages = [('preselection', None)] + [(c['name'], j) for j, c in enumerate(cl)]
    for name, j in stages:
        if j is not None:
            m_mc &= p_mc[:, j]
            m_on &= p_on[:, j]
        r = report(name if j is None else '+ ' + name, m_mc, m_on)
        kg = r[plain(GROUPS[SINGLE_G][0])] / tot[SINGLE_G] if tot[SINGLE_G] else float('nan')
        ks = r[plain(GROUPS[SECONDARY][0])] / tot[SECONDARY] if tot[SECONDARY] else float('nan')
        cols = ''.join(f'{r[plain(GROUPS[i][0])]:9.1f}' for i in SHOWN)
        print(f'  {r["stage"]:<38} {r["mc_expected"]:8.1f} {r["data_n"]:6d} {r["z"]:5.1f} |{cols} | '
              f'{kg:8.3f} {ks:8.4f} {kg / ks if ks else float("inf"):7.1f}')

    print('\nN-1: all other cuts applied; what this cut alone then removes (fraction of each group it rejects)')
    print(f'  {"cut left out":<38} {"MC":>8} {"data":>6} {"z":>5} | {"rejects: 1g":>11} {"sec":>6} '
          f'{"cosmic":>7} {"other nu":>8} {"MC all":>7} {"data":>6}')
    for j, c in enumerate(cl):
        others_mc = np.all(np.delete(p_mc, j, axis=1), axis=1)
        others_on = np.all(np.delete(p_on, j, axis=1), axis=1)
        r = report('N-1 ' + c['name'], others_mc, others_on)

        def rej(sel):
            sw = w[others_mc & sel].sum()
            return 1 - w[others_mc & sel & p_mc[:, j]].sum() / sw if sw else float('nan')
        nd = int(others_on.sum())
        rd = 1 - int((others_on & p_on[:, j]).sum()) / nd if nd else float('nan')
        print(f'  {c["name"]:<38} {r["mc_expected"]:8.1f} {r["data_n"]:6d} {r["z"]:5.1f} | {rej(g == SINGLE_G):11.3f} '
              f'{rej(g == SECONDARY):6.3f} {rej(g == COSMIC):7.3f} {rej(g == OTHER):8.3f} '
              f'{rej(np.ones(len(mc), bool)):7.3f} {rd:6.3f}')
    pd.DataFrame(rows).to_csv(os.path.join(a.out, 'cutflow.csv'), index=False)

    plots(mc, on, p_mc, p_on, cl, a.out, pot_on, region)
    print(f'\nwrote {a.out}/cutflow.csv and {a.out}/{{before,after,nm1}}/*.png,pdf')
    return 0


def plots(mc, on, p_mc, p_on, cl, out, pot_on, region):
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
    all_mc, all_on = np.all(p_mc, axis=1), np.all(p_on, axis=1)
    ke = dict(var=KE, range=(0.0, 2000.0), nbins=20, xlabel='Shower KE [MeV]', lines=[], name='shower KE')
    jobs = []
    for c in cl + [ke]:
        lines = c['lines'] if c is not ke else []
        jobs.append(('before', c, np.ones(len(mc), bool), np.ones(len(on), bool), lines, 'preselection'))
        jobs.append(('after', c, all_mc, all_on, lines, 'all cuts'))
    for j, c in enumerate(cl):
        jobs.append(('nm1', c, np.all(np.delete(p_mc, j, axis=1), axis=1),
                     np.all(np.delete(p_on, j, axis=1), axis=1), c['lines'], f'all cuts but {c["name"]}'))
    for sub, c, m_mc, m_on, lines, stage in jobs:
        d = os.path.join(out, sub)
        os.makedirs(d, exist_ok=True)
        fig, ax = plt.subplots(figsize=(10, 7.5))
        lo, hi = c['range']
        edges = np.linspace(lo, hi, c['nbins'] + 1)
        sel, dsel = mc[m_mc], on[m_on]
        xs, ws = [], []
        for i in range(len(GROUPS)):
            s = sel[sel['g'] == i]
            xs.append(np.clip(s[c['var']].to_numpy(dtype=float), lo, hi - 1e-9 * (hi - lo)))
            ws.append(s['w'].to_numpy(dtype=float))
        keep = [i for i in range(len(GROUPS)) if ws[i].sum() > 0]
        if keep:
            ax.hist([xs[i] for i in keep], bins=edges, weights=[ws[i] for i in keep], stacked=True,
                    histtype='barstacked', color=[GROUPS[i][1] for i in keep],
                    label=[f'{GROUPS[i][0]} ({ws[i].sum():.1f})' for i in keep])
        x_on = np.clip(dsel[c['var']].to_numpy(dtype=float), lo, hi - 1e-9 * (hi - lo))
        k = np.histogram(x_on, edges)[0]
        ctr = 0.5 * (edges[1:] + edges[:-1])
        ax.errorbar(ctr, k, yerr=np.sqrt(k), fmt='o', color='k', ms=4, label=f'Data (BNB) ({k.sum()})')
        ax.set_xlim(lo, hi)
        ax.set_ylim(0, 1.5 * max(ax.get_ylim()[1], 1.0))
        trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        span = 0.06 * (hi - lo)
        for n, (x, side) in enumerate(lines):
            ax.axvline(x, color='r', ls='--', label='cut (arrow: kept)' if n == 0 else None)
            ax.annotate('', xy=(x + side * span, 0.6), xytext=(x, 0.6), xycoords=trans, textcoords=trans,
                        arrowprops=dict(arrowstyle='->', color='r', lw=2))
        ax.set_xlabel(c['xlabel'])
        ax.set_ylabel('Candidates')
        ax.text(0.0, 1.01, PRELIM, transform=ax.transAxes, ha='left', va='bottom', color='red',
                fontweight='bold', fontsize='medium')
        ax.text(1.0, 1.01, f'{pot_on:.2e} POT', transform=ax.transAxes, ha='right', va='bottom', fontsize='medium')
        ax.legend(loc='upper right', fontsize='x-small', frameon=False, title=f'{region}\n{stage}',
                  title_fontsize='x-small')
        fig.tight_layout()
        for ext in ('png', 'pdf'):
            fig.savefig(os.path.join(d, f'{c["var"]}.{ext}'))
        plt.close(fig)


if __name__ == '__main__':
    sys.exit(main())
