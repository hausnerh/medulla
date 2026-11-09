#!/usr/bin/env python3
"""
gOre data/MC plots with the dedicated NC Δ radiative sample mixed into the CV.

Reads the same spineplot configs as the standard plots (categories, labels,
colours, variables, binning, figure list, style) and draws the same stacked
spectra, but with the CV + ncdeltarad combination of batch/gore_cut_opt.py
(per category: CV NC RES rows x alpha, ncdeltarad rows x sidecar weight x POT
scale x (1 - alpha); alpha from the stage1 tree, fixed for every stage) and
optionally the dirt MC as its own category.

Errors: the hatched band is the MC statistical error, sqrt(sum w^2) per bin;
NO systematic band yet (the ncdeltarad rows carry no universe weights until a
systematics pass is run on them). Data points: sqrt(N).

By default MC is scaled to the onbeam (10% data) POT and the data overlaid;
--pot X scales MC to X POT and drops the data (e.g. the full analysis POT).

Usage
  plot_combined.py --sys build/output_gOre_1g1p_signal_osc_run2_sys.root \\
      --ncdr build/gOre_1g1p_signal_ncdeltarad_run2.root --ncdr-weights <base>/gen_sidecar/ncdeltarad_weights.csv \\
      [--dirt build/gOre_1g1p_signal_osc_run2_dirt.root] --out build/plots_ncdr \\
      spineplot/configurations/analyses/icarus/gOre_1g1p_stage{1,2,3}_datamc.toml
"""
import argparse
import os
import sys

import numpy as np

try:
    import tomllib

    def toml_loads(s):
        return tomllib.loads(s)
except ImportError:                              # python < 3.11 (the gpvm's 3.9)
    import toml

    def toml_loads(s):
        return toml.loads(s)

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
SPINEPLOT = os.path.join(HERE, '..', '..', 'spineplot')
DIRT_LABEL, DIRT_COLOR = r'Dirt (out-of-FV $\nu$)', '#8C564B'


def read_cfg(path):
    cfg = toml_loads(open(path).read())
    styles = {}
    for inc in cfg.get('this_includes', []):
        f = os.path.join(SPINEPLOT, inc['file'])
        styles.update({s['name']: s for s in toml_loads(open(f).read()).get('style', [])})
    cfg['_styles'] = styles
    return cfg


def tree_df(f, key, tree):
    t = f[f'events/{key}/{tree}']
    return t.arrays([b for b in t.keys() if not b.endswith('_exposure')], library='pd')


def precompute(df, exprs):
    for name, expr in (exprs or {}).items():
        if name != 'true_category':
            df[name] = df.eval(expr)
    return df


def build(f, nd_file, dirt_file, side, tree, exprs, pots, alpha):
    """All rows of one tree with weights at the onbeam POT and a category."""
    import pandas as pd
    import gore_cut_opt as G
    pot_on, pot_cv, pot_nd, lt_on, lt_off = (pots[k] for k in ('on', 'cv', 'nd', 'lt_on', 'lt_off'))
    cv = tree_df(f, 'full', tree)
    cv['true_category'] = cv['true_category'].fillna(9).astype(int)
    cv['w'] = pot_on / pot_cv
    nd = tree_df(nd_file, 'ncdeltarad', tree)
    cv, nd, _ = G.combine(cv, nd, side, pot_cv, pot_nd, pot_on, alpha=alpha)
    parts = [cv, nd]
    if dirt_file is not None:
        dr = tree_df(dirt_file, 'dirt', tree)
        c = dr['true_category'].to_numpy(dtype=float)
        dr = dr[np.isfinite(c) & (c != 9)].copy()
        dr['true_category'] = G.DIRT
        dr['w'] = pot_on / pots['dirt']
        parts.append(dr)
    off = tree_df(f, 'offbeam', tree)
    off['true_category'] = 9
    off['w'] = lt_on / lt_off
    parts.append(off)
    mc = pd.concat([precompute(p, exprs) for p in parts], ignore_index=True)
    data = precompute(tree_df(f, 'onbeam', tree), exprs)
    return mc, data


def draw(mc, data, var, vcfg, labels, colors, style, pot, out, show_data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import gore_cut_opt as G
    sheet = os.path.join(SPINEPLOT, style.get('style_sheet', 'styles/muon2024.mplstyle'))
    plt.style.use(sheet if os.path.exists(sheet) else 'default')
    lo, hi = vcfg['range']
    edges = np.linspace(lo, hi, int(vcfg['nbins']) + 1)
    cats = [c for c in list(range(10)) + [G.DIRT] if (mc['true_category'] == c).any()]
    xs = [mc.loc[mc['true_category'] == c, var].to_numpy(dtype=float) for c in cats]
    ws = [mc.loc[mc['true_category'] == c, 'w'].to_numpy(dtype=float) for c in cats]
    tot = sum(np.histogram(x, edges, weights=w)[0] for x, w in zip(xs, ws))
    err = np.sqrt(sum(np.histogram(x, edges, weights=w ** 2)[0] for x, w in zip(xs, ws)))
    den = tot.sum() or 1.0
    lab = [f'{DIRT_LABEL if c == G.DIRT else labels[c]} '
           f'({np.histogram(x, edges, weights=w)[0].sum():.1f}, {np.histogram(x, edges, weights=w)[0].sum() / den:.1%})'
           for c, x, w in zip(cats, xs, ws)]
    col = [DIRT_COLOR if c == G.DIRT else colors[c] for c in cats]
    fig, ax = plt.subplots(figsize=style.get('default_figsize', [9.0, 6.0]))
    ax.hist(xs, bins=edges, weights=ws, stacked=True, histtype='barstacked', color=col, label=lab)
    ax.fill_between(edges, np.r_[tot - err, (tot - err)[-1]], np.r_[tot + err, (tot + err)[-1]], step='post',
                    facecolor='none', edgecolor='k', hatch='///', lw=0, label='MC stat.')
    ymax = tot.max() if len(tot) else 1.0
    if show_data:
        n = np.histogram(data[var].to_numpy(dtype=float), edges)[0]
        ctr = 0.5 * (edges[1:] + edges[:-1])
        ax.errorbar(ctr, n, yerr=np.sqrt(n), fmt='o', color='k', ms=4, label=f'{labels[10]} ({n.sum()})')
        ymax = max(ymax, (n + np.sqrt(n)).max())
    ax.set_xlim(lo, hi)
    ax.set_ylim(0, 1.5 * ymax if ymax > 0 else 1)
    ax.set_xlabel(vcfg.get('xlabel', var))
    ax.set_ylabel('Candidates')
    ax.set_title(style.get('title', ''))
    ax.text(0.02, 0.97, f'{pot:.2e} POT', transform=ax.transAxes, va='top', fontsize='small')
    if style.get('mark_preliminary'):
        ax.text(0.02, 0.90, style['mark_preliminary'], transform=ax.transAxes, va='top', fontsize='small')
    ax.legend(loc='upper right', fontsize='x-small', ncol=2, frameon=False)
    fig.tight_layout()
    for ext in ('png', 'pdf'):
        fig.savefig(f'{out}.{ext}')
    plt.close(fig)


def main(argv=None):
    import uproot
    import gore_cut_opt as G
    import sample_check as SC
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('configs', nargs='+', help='spineplot data/MC configs (their mc tree is plotted)')
    ap.add_argument('--sys', required=True, help='sys ROOT with events/{full,offbeam,onbeam}')
    ap.add_argument('--ncdr', required=True)
    ap.add_argument('--ncdr-weights', required=True)
    ap.add_argument('--dirt', help='medulla output with the dirt MC (events/dirt)')
    ap.add_argument('--pot', type=float, help='scale MC to this POT and drop the data (default: onbeam POT + data)')
    ap.add_argument('--out', required=True)
    a = ap.parse_args(argv)
    if not a.dirt:
        print('WARNING: no --dirt: the out-of-FV background is missing from these plots')

    f, nd_file = uproot.open(a.sys), uproot.open(a.ncdr)
    dirt_file = uproot.open(a.dirt) if a.dirt else None
    pots = dict(on=G.hist_sum(f, 'events/onbeam/POT'), cv=G.hist_sum(f, 'events/full/POT'),
                nd=G.hist_sum(nd_file, 'events/ncdeltarad/POT'),
                lt_on=G.hist_sum(f, 'events/onbeam/Livetime'), lt_off=G.hist_sum(f, 'events/offbeam/Livetime'))
    if dirt_file is not None:
        pots['dirt'] = G.hist_sum(dirt_file, 'events/dirt/POT')
    side = SC.read_weights(a.ncdr_weights)

    # alpha from stage1 (most statistics), then fixed for every tree
    cv1 = tree_df(f, 'full', 'selected_1g1p_stage1')
    cv1['true_category'] = cv1['true_category'].fillna(9).astype(int)
    cv1['w'] = 1.0
    _, _, tab = G.combine(cv1, tree_df(nd_file, 'ncdeltarad', 'selected_1g1p_stage1'), side,
                          pots['cv'], pots['nd'], pots['on'])
    alpha = tab['alpha'].to_numpy()
    print('alpha per category (stage1):', dict(zip(G.CATS, np.round(alpha, 3))))

    scale = (a.pot / pots['on']) if a.pot else 1.0
    pot_shown = a.pot or pots['on']
    for path in a.configs:
        cfg = read_cfg(path)
        an = cfg['analysis']
        tree = cfg['samples']['mc']['trees'][0]
        mc, data = build(f, nd_file, dirt_file, side, tree, cfg['samples']['mc'].get('precompute'), pots, alpha)
        mc['w'] *= scale
        stem = os.path.join(a.out, os.path.splitext(os.path.basename(path))[0])
        os.makedirs(stem, exist_ok=True)
        y = mc.groupby('true_category')['w'].sum()
        print(f'\n{os.path.basename(path)} ({tree}) at {pot_shown:.3g} POT: MC {mc["w"].sum():.2f} '
              f'(signal {y.get(0, 0.0):.2f}), data {len(data) if not a.pot else "-"}')
        for c, v in y.items():
            print(f'    {G.CATS[c]:<24} {v:8.2f}')
        for fig in cfg.get('figure', []):
            var = fig['artists'][0]['variable']
            if var not in cfg['variables'] or var not in mc.columns:
                print(f'  skip {var}: not in the variables/tree')
                continue
            draw(mc, data, var, cfg['variables'][var], an['category_labels'], an['category_colors'],
                 cfg['_styles'].get(fig.get('style'), {}), pot_shown, os.path.join(stem, fig['name']), not a.pot)
        print(f'  wrote {stem}/*.png, *.pdf')
    return 0


if __name__ == '__main__':
    sys.exit(main())
