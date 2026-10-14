#!/usr/bin/env python3
"""
S/sqrt(S+B) re-optimization of the gOre 1γ1p stage2 (pi0_rejection) and
stage3 (egamma_separation) thresholds, from a merged signal sys output.

The stage2/3 cuts are emulated on the stage1 (preselection) tree, mirroring
selection/include/gOre/cuts_gOre.h (strict inequalities, NaN fails, no
subleading shower); the emulation is checked against the stage2/3 tree
counts before anything is optimized.

S = true_category 0 (NC Δ→Nγ 1γ1p). B = every other MC category plus the
off-beam cosmics. Weights match spineplot: MC scaled by POT and off-beam by
livetime to the onbeam exposure. The argmax does not depend on the exposure.

Methods
  seq   the analysisTools procedure (run_analysis_gOre.cc): one cut at a time
        in the original order and ranges, 1000 steps, each on top of the
        previous ones, starting from preselection + no subleading shower.
  iter  coordinate ascent: each cut re-optimized with all the others at their
        current values (N-1), swept until nothing moves.
  2d    (γ primary softmax, Δ-mass halfwidth) grid with the other five cuts at
        their current values; FOM 3x3-smoothed before the argmax.
On a 1D FOM plateau the plateau midpoint is taken rather than its edge.

Every method is validated by random split-halves: optimize on one half, then
evaluate on the other, relative to the current cuts on the same half. With
O(100) stage1 signal MC events this is what separates real gain from
overtraining.

Usage
  python batch/gore_cut_opt.py output_gOre_1g1p_sys.root --out cutopt
  python batch/gore_cut_opt.py ... --fom purity   # sqrt(purity)*S/sqrt(S+B)
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
import uproot

CATS = ['NC Δ→Nγ (1γ1p)', 'NC Δ→Nγ (other topo)', 'Other NC 1γXp post-FSI',
        'NC π0 (Δ res)', 'NC π0 (other)', 'NC π±', 'Other NC', 'CC νe',
        'Other CC', 'Cosmic']

# (name, branch, kind, current threshold, analysisTools scan range);
# kind 'lo' passes x > t, 'hi' passes x < t. Order = analysisTools order.
CUTS = [
    ('gamma_primary_softmax',  'reco_leading_primary_gOre_primary_softmax',    'lo', 0.999100, (0.99, 1.0)),
    ('proton_primary_softmax', 'reco_leading_primary_proton_primary_softmax',  'lo', 0.951300, (0.9, 1.0)),
    ('delta_mass_halfwidth',   'dm_dev',                                       'hi', 59.7,     (0.0, 300.0)),
    ('start_dedx',             'reco_leading_primary_gOre_start_dedx',         'lo', 4.07,     (0.0, 5.0)),
    ('axial_spread',           'reco_leading_primary_gOre_axial_spread',       'lo', -0.5,     (-0.5, 1.0)),
    ('directional_spread',     'reco_leading_primary_gOre_directional_spread', 'hi', 0.037,    (0.0, 1.0)),
    ('photon_softmax',         'reco_leading_primary_gOre_photon_softmax',     'lo', 0.9225,   (0.9, 1.0)),
]
NAMES = [c[0] for c in CUTS]
STAGE2 = NAMES[:3]
NOMINAL = {c[0]: c[3] for c in CUTS}
OPEN = {c[0]: (-1e9 if c[2] == 'lo' else 1e9) for c in CUTS}
NSTEP = 1000
G2D = np.round(np.arange(0.9960, 0.99981, 0.0001), 5)
H2D = np.arange(30.0, 201.0, 5.0)


def hist_sum(f, key):
    return float(f[key].values().sum())


def load(path):
    f = uproot.open(path)
    pot_on = hist_sum(f, 'events/onbeam/POT')
    lt_on = hist_sum(f, 'events/onbeam/Livetime')
    want = {c[1] for c in CUTS if c[1] != 'dm_dev'} | {
        'reco_delta_mass_P', 'reco_subleading_primary_gOre_shower_ke', 'true_category'}
    out = {}
    for s, scale in (('full', pot_on / hist_sum(f, 'events/full/POT')),
                     ('offbeam', lt_on / hist_sum(f, 'events/offbeam/Livetime')),
                     ('onbeam', 1.0)):
        t = f[f'events/{s}/selected_1g1p_stage1']
        df = t.arrays([b for b in want if b in t.keys()], library='pd')
        if s != 'full':
            df['true_category'] = 9 if s == 'offbeam' else -1
        df['true_category'] = df['true_category'].fillna(9).astype(int)
        df['dm_dev'] = (df['reco_delta_mass_P'] - 1232.0).abs()
        df['no_sub'] = ~(df['reco_subleading_primary_gOre_shower_ke'] > 0)
        df['w'] = scale
        out[s] = df
    return f, out, pot_on


def passes(df, thr, names=NAMES):
    m = df['no_sub'].to_numpy().copy()
    for name, br, kind, *_ in CUTS:
        if name not in names:
            continue
        x = df[br].to_numpy()
        with np.errstate(invalid='ignore'):
            m &= (x > thr[name]) if kind == 'lo' else (x < thr[name])
    return m


class Optimizer:
    def __init__(self, ev, sig, fom_kind):
        self.ev, self.sig, self.kind = ev, sig, fom_kind

    def fom(self, s, b):
        s, b = np.asarray(s, float), np.asarray(b, float)
        f = np.where(s + b > 0, s / np.sqrt(np.maximum(s + b, 1e-300)), 0.0)
        if self.kind == 'purity':
            f = f * np.sqrt(np.where(s + b > 0, s / np.maximum(s + b, 1e-300), 0.0))
        return f

    def summary(self, w, thr):
        m = passes(self.ev, thr)
        S, B = float((w * self.sig)[m].sum()), float((w * ~self.sig)[m].sum())
        return dict(S=S, B=B, fom=float(self.fom(S, B)), eff=S / float((w * self.sig).sum()),
                    purity=S / (S + B) if S + B else 0.0,
                    nS=int((self.sig & m & (w > 0)).sum()), nB=int((~self.sig & m & (w > 0)).sum()))

    def scan(self, w, thr, name, active):
        base = passes(self.ev, thr, [n for n in active if n != name])
        _, br, kind, _, (lo, hi) = next(c for c in CUTS if c[0] == name)
        x = self.ev[br].to_numpy()
        grid = lo + (hi - lo) * np.arange(NSTEP) / NSTEP
        with np.errstate(invalid='ignore'):
            P = (x[None, :] > grid[:, None]) if kind == 'lo' else (x[None, :] < grid[:, None])
        P &= base[None, :]
        S = P @ (w * self.sig)
        B = P @ (w * ~self.sig)
        return grid, S, B, self.fom(S, B)

    @staticmethod
    def best(grid, F):
        i = j = int(np.argmax(F))
        while j + 1 < len(F) and np.isclose(F[j + 1], F[i], rtol=1e-12, atol=0):
            j += 1
        return 0.5 * (grid[i] + grid[j])

    def seq(self, w):
        thr, active = dict(OPEN), []
        for name in NAMES:
            active.append(name)
            thr[name] = self.best(*self.scan(w, thr, name, active)[::3])
        return thr

    def iter(self, w, start=NOMINAL, maxit=20):
        thr = dict(start)
        for _ in range(maxit):
            moved = False
            for name in NAMES:
                t = self.best(*self.scan(w, thr, name, NAMES)[::3])
                if not np.isclose(t, thr[name]):
                    thr[name], moved = t, True
            if not moved:
                break
        return thr

    def grid2d(self, w):
        others = [n for n in NAMES if n not in ('gamma_primary_softmax', 'delta_mass_halfwidth')]
        base = passes(self.ev, NOMINAL, others)
        gps = self.ev['reco_leading_primary_gOre_primary_softmax'].to_numpy()
        dmd = self.ev['dm_dev'].to_numpy()
        with np.errstate(invalid='ignore'):
            P = base[None, None, :] & (gps[None, None, :] > G2D[:, None, None]) & (dmd[None, None, :] < H2D[None, :, None])
        S, B = P @ (w * self.sig), P @ (w * ~self.sig)
        return S, B, self.fom(S, B)

    def two_d(self, w):
        F = self.grid2d(w)[2]
        Pd = np.pad(F, 1, mode='edge')
        Fs = sum(Pd[1 + a:1 + a + F.shape[0], 1 + b:1 + b + F.shape[1]] for a in (-1, 0, 1) for b in (-1, 0, 1)) / 9
        i, j = np.unravel_index(np.argmax(Fs), Fs.shape)
        thr = dict(NOMINAL)
        thr['gamma_primary_softmax'], thr['delta_mass_halfwidth'] = float(G2D[i]), float(H2D[j])
        return thr


def composition(ev, w, thr):
    m = passes(ev, thr)
    out = []
    for c in range(len(CATS)):
        mc = m & (ev['true_category'] == c).to_numpy()
        if mc.any():
            out.append(f'{CATS[c]} {w[mc].sum():.3f} ({int(mc.sum())})')
    return '; '.join(out)


def plots(opt, w, sets, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for setname in ('nominal', 'two_d'):
        thr = sets[setname]
        fig, axs = plt.subplots(2, 4, figsize=(20, 9))
        for ax, name in zip(axs.flat, NAMES):
            g, S, B, F = opt.scan(w, thr, name, NAMES)
            ax.plot(g, S / S.max() if S.max() else S, label='Signal (N-1 rel.)', color='tab:orange')
            ax.plot(g, B / B.max() if B.max() else B, label='Background (N-1 rel.)', color='tab:blue')
            ax.plot(g, np.where(S + B > 0, S / np.maximum(S + B, 1e-300), np.nan), label='Purity', color='tab:green')
            ax.plot(g, F / F.max() if F.max() else F, label='FOM (norm.)', color='black', lw=2)
            ax.axvline(NOMINAL[name], color='gray', ls='--', label='current')
            ax.axvline(sets['two_d'][name], color='red', ls=':', label='2D re-tune')
            ax.set_title(name)
            ax.set_ylim(0, 1.05)
        axs.flat[-1].axis('off')
        axs.flat[0].legend(fontsize=8)
        fig.suptitle(f'N-1 scans, other cuts at the {setname} thresholds')
        fig.tight_layout()
        fig.savefig(f'{out}/nm1_scans_{setname}.png', dpi=90)
        plt.close(fig)
    S, B, F = opt.grid2d(w)
    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.pcolormesh(H2D, G2D, F, shading='nearest')
    fig.colorbar(im, label='FOM at the onbeam exposure')
    ax.plot(NOMINAL['delta_mass_halfwidth'], NOMINAL['gamma_primary_softmax'], 'wx', ms=10, label='current')
    ax.plot(sets['two_d']['delta_mass_halfwidth'], sets['two_d']['gamma_primary_softmax'], 'r*', ms=14, label='2D re-tune')
    ax.set_xlabel('|M_Δ − 1232| halfwidth [MeV]')
    ax.set_ylabel('γ primary softmax lower bound')
    ax.legend()
    fig.tight_layout()
    fig.savefig(f'{out}/fom2d_softmax_dmass.png', dpi=90)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('sysroot', help='merged signal output (events/{full,offbeam,onbeam}/selected_1g1p_stage*)')
    ap.add_argument('--out', default='cutopt')
    ap.add_argument('--fom', choices=['ssb', 'purity'], default='ssb',
                    help="ssb: S/sqrt(S+B); purity: sqrt(purity)*S/sqrt(S+B) (last analysisTools get_fom)")
    ap.add_argument('--nsplit', type=int, default=20, help='random split-halves for validation')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    pd.set_option('display.width', 250)
    pd.set_option('display.max_columns', 30)

    f, d, pot_on = load(a.sysroot)
    print(f'onbeam exposure {pot_on:.4g} POT')
    for s, df in d.items():
        n2, n3 = (f[f'events/{s}/selected_1g1p_stage{k}'].num_entries for k in (2, 3))
        e2, e3 = int(passes(df, NOMINAL, STAGE2).sum()), int(passes(df, NOMINAL).sum())
        ok = (e2, e3) == (n2, n3)
        print(f'  emulation {s:8s}: stage1 {len(df)}, stage2 {e2}/{n2}, stage3 {e3}/{n3}  {"OK" if ok else "MISMATCH"}')
        if not ok:
            raise SystemExit('cut emulation does not reproduce the stage trees -- update CUTS to the selection toml')

    ev = pd.concat([d['full'], d['offbeam']], ignore_index=True)
    sig = (ev['true_category'] == 0).to_numpy()
    w = ev['w'].to_numpy()
    opt = Optimizer(ev, sig, a.fom)

    sets = {'nominal': dict(NOMINAL), 'seq': opt.seq(w), 'iter': opt.iter(w), 'two_d': opt.two_d(w)}
    tab = pd.DataFrame([{'set': k, **t, **opt.summary(w, t)} for k, t in sets.items()])
    print(f'\nthresholds and yields at {pot_on:.3g} POT (FOM = {a.fom}):')
    print(tab.to_string(index=False, float_format=lambda v: f'{v:.5g}'))
    print('\ncomposition, weighted (raw MC):')
    for k, t in sets.items():
        print(f'  {k:8s}: {composition(ev, w, t)}')

    rng = np.random.default_rng(7)
    gains = {k: [] for k in ('seq', 'iter', 'two_d')}
    for _ in range(a.nsplit):
        tr = rng.random(len(ev)) < 0.5
        wt, we = np.where(tr, 2 * w, 0.0), np.where(~tr, 2 * w, 0.0)
        nom = opt.summary(we, NOMINAL)['fom']
        for k, fn in (('seq', opt.seq), ('iter', opt.iter), ('two_d', opt.two_d)):
            gains[k].append(opt.summary(we, fn(wt))['fom'] / nom if nom else np.nan)
    print(f'\nsplit-half validation ({a.nsplit} random halves): held-out FOM / current cuts on the same half')
    for k, g in gains.items():
        q = np.nanpercentile(g, [16, 50, 84])
        print(f'  {k:6s} median {q[1]:.3f}  [16%, 84%] = [{q[0]:.3f}, {q[2]:.3f}]  worse than current in {np.mean(np.array(g) < 1):.0%}')

    tab.to_csv(f'{a.out}/thresholds.csv', index=False)
    json.dump({'fom': a.fom, 'sets': sets, 'split_gain': gains}, open(f'{a.out}/thresholds.json', 'w'), indent=1)
    plots(opt, w, sets, a.out)
    print(f'\nwrote {a.out}/thresholds.csv, thresholds.json, nm1_scans_*.png, fom2d_softmax_dmass.png')


if __name__ == '__main__':
    main()
