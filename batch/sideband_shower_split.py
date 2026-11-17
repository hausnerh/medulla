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
true_n_nu_photons, true_n_nu_electrons, true_n_secondary_true_showers,
true_n_secondary_true_showers_in_av, and the pi0 photons among the secondaries
(parent pdg 111): true_n_pi0_true_showers, true_n_pi0_true_showers_escaping
(starts outside the active volume or not contained), true_n_pi0_lone_photons
(a pi0 with only one decay photon in the truth record):
  true single gamma primary      1 photon from the nu, no other true shower
  true single e primary          1 electron from the nu, no other true shower
  true multiple primary e/gamma  >=2 showers from the nu, or 1 plus secondaries
  true pi0, all photons in argon only secondary showers, >=1 from a pi0, every
                                 pi0 photon starting and contained in the argon
  true pi0, a photon escapes     as above, but a pi0 photon starts outside or
                                 leaves the argon, or is missing from the record
  true other secondary showers   only secondary showers, none from a pi0
                                 (delta rays, Michel electrons, brems, ...)
  other nu                       no true shower above threshold
Files without the pi0 branches fall back to the split by start point only
(secondary, all start in argon / >=1 starts outside) with no "other secondary".
MC rows of readouts that fail the trigger emulation (event_trigger_emulated = 0)
are dropped: such a readout exists in data only if a cosmic fired the trigger,
which the off-beam sample already counts (--no-trigger-emulation keeps them).
  dirt                           out-of-FV neutrinos (dirt MC, matched rows)
  cosmics                        one category from two sources (the console
                                 prints the split): off-beam data (beam-off
                                 triggers, scaled by livetime) and the cv rows
                                 without a neutrino match (cosmics recorded in
                                 neutrino-triggered readouts, which off-beam
                                 cannot represent; POT-scaled). This sideband has
                                 no flash cut, so out-of-time cosmics enter it.
                                 --no-overlay-cosmics drops the second source
                                 (the old matched-only scheme).

Output (--out):
  cut flow      expected events per truth group, data, Poisson z, and the
                fraction of true single gamma vs secondary showers kept
                (plot binning: the cut values fall on bin edges)
  N-1 table     with all other cuts applied, the fraction of each group the
                left-out cut rejects (its own contribution), MC and data
  before/       each cut variable at the preselection
  after/        each cut variable after all cuts
  nm1/          each cut variable with all the other cuts applied
  cutflow.csv

Input: a selection-only medulla output (events/cv, offbeam, onbeam[, dirt]); a
sys ROOT (events/full, offbeam, onbeam) works with --no-overlay-cosmics, since
run_systematics moves the overlay cosmics to *_nonmatched.

Usage
  sideband_shower_split.py build/output_gOre_1g1p_sidebands_osc_run2_origin.root --out build/sideband_cuts \\
      [--gps 0.9985 --dm-halfwidth 95] [--no-overlay-cosmics] \\
      [--weights build/output_gOre_1g1p_sidebands_osc_run2_origin_sys.root]
--weights takes the run_systematics output of the same selection (systematics/toml/
gOre_1g1p_sidebands_osc_run2_origin.toml, `run_cc_sideband.sh --sys-only`) and adds the flux +
cross section band (spineplot's method: per-knob variance around the CV, summed) and MC statistics
to every cut-flow stage, with the data/MC significance with and without them (sysband.csv).
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
ORIGIN = ['true_n_nu_photons', 'true_n_nu_electrons', 'true_n_secondary_true_showers',
          'true_n_secondary_true_showers_in_av']
PI0 = ['true_n_pi0_true_showers', 'true_n_pi0_true_showers_escaping', 'true_n_pi0_lone_photons']
TRIG = ('event_trigger_emulated', 'trigger_emulated')
GROUPS = [('True single $\\gamma$ primary', '#E69900'),
          ('True single $e$ primary', '#2CA02C'),
          ('True multiple primary $e/\\gamma$', '#1F77B4'),
          ('True $\\pi^0$ showers, all photons in argon', '#F2E640'),
          ('True $\\pi^0$ showers, a photon escapes the argon', '#17BECF'),
          ('True other secondary showers ($\\delta$-ray, Michel, ...)', '#E377C2'),
          ('Other $\\nu$', '#9980CC'),
          ('Cosmics', '#808080'),
          ('Dirt', '#8C564B')]
SINGLE_G, SINGLE_E, MULTI, SEC_IN, SEC_OUT, SEC_OTHER, OTHER, COSMIC, DIRT = range(9)
LEGACY_LABELS = {SEC_IN: 'True secondary showers, all start in argon',
                 SEC_OUT: 'True secondary showers, $\\geq$1 starts outside argon'}
PI0_LABELS = {i: GROUPS[i][0] for i in LEGACY_LABELS}
SHOWN = list(range(len(GROUPS)))
SHORT = ['1g', '1e', 'multi', 'pi0 in', 'pi0 out', 'oth sec', 'oth nu', 'cosmic', 'dirt']
KE = 'reco_leading_primary_gOre_shower_ke'
PRELIM = 'ICARUS Analysis In-Progress'


def plain(label):
    for tex, txt in (('$\\geq$', '>='), ('$\\pi^0$', 'pi0'), ('$\\delta$', 'delta'), ('$e/\\gamma$', 'e/gamma'), ('$\\gamma$', 'gamma'), ('$e$', 'e'), ('$\\nu$', 'nu')):
        label = label.replace(tex, txt)
    return label


def group_of(d):
    """Truth group per row from the shower-origin counts (NaN counts: no match -> other nu).
    With the pi0 branches the secondary-only rows split into pi0 in argon / pi0 escaping / other secondary;
    without them (older files), by whether every secondary shower starts in the argon."""
    ng, ne, ns, ns_in = (np.nan_to_num(d[c].to_numpy(dtype=float), nan=0.0) for c in ORIGIN)
    nnu = ng + ne
    sec = (nnu == 0) & (ns > 0)
    g = np.full(len(d), OTHER)
    if all(c in d for c in PI0):
        npi, esc, lone = (np.nan_to_num(d[c].to_numpy(dtype=float), nan=0.0) for c in PI0)
        g[sec & (npi == 0)] = SEC_OTHER
        g[sec & (npi > 0)] = SEC_IN
        g[sec & (npi > 0) & ((esc > 0) | (lone > 0))] = SEC_OUT
    else:
        g[sec & (ns_in >= ns)] = SEC_IN
        g[sec & (ns_in < ns)] = SEC_OUT
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


def aligned_edges(lo, hi, nbins, cuts):
    """Bin edges close to linspace(lo, hi, nbins + 1) with every cut value on an edge. One cut: the width is
    adjusted so lo and the cut are both edges; several: the width divides the gaps between the cuts and the
    range is extended to the next edge below lo / above hi."""
    w0 = (hi - lo) / nbins
    cuts = sorted(c for c in cuts if lo < c < hi)
    if not cuts:
        return np.linspace(lo, hi, nbins + 1)
    if len(cuts) == 1:
        w = (cuts[0] - lo) / max(1, round((cuts[0] - lo) / w0))
        start = lo
    else:
        gaps = np.diff(cuts)
        k = max(1, round(gaps[0] / w0))
        w = gaps[0] / k
        start = cuts[0] - np.ceil((cuts[0] - lo) / w - 1e-9) * w
    n = int(np.ceil((hi - start) / w - 1e-9))
    return start + w * np.arange(n + 1)


def has(path, key):
    import uproot
    f = uproot.open(path)
    return 'events' in f and key in f['events']


def load(path, key, want, optional=()):
    """events/<key>/TREE as a DataFrame (`want` required, `optional` read when present), POT, Livetime."""
    import uproot
    f = uproot.open(path)
    if 'events' not in f or key not in f['events']:
        raise SystemExit(f'{path}: no events/{key}')
    t = f[f'events/{key}/{TREE}']
    want = list(want) + [b for b in optional if b in t.keys() and b not in want]
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


def triggered(d, label, keep_all):
    """Drop MC rows whose readout failed the trigger emulation (prints how many); unchanged without the branch."""
    col = next((c for c in TRIG if c in d), None)
    if col is None:
        print(f'  {label}: no trigger-emulation branch ({TRIG[0]}): untriggered MC readouts are kept')
        return d
    fail = ~(d[col].to_numpy(dtype=float) == 1)
    print(f'  {label}: {int(fail.sum())} of {len(d)} rows ({fail.mean() if len(d) else 0:.1%}) fail the trigger '
          f'emulation' + (' (kept: --no-trigger-emulation)' if keep_all else ' (dropped)'))
    return d if keep_all else d[~fail].copy()


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


def poisson_tail(n, mu, upper):
    """P(N >= n) if `upper`, else P(N <= n), for N ~ Poisson(mu)."""
    from math import exp, lgamma, log
    if mu <= 0:
        return float(n == 0) if upper else 1.0
    pmf = lambda k: exp(k * log(mu) - mu - lgamma(k + 1))
    if upper:
        return max(0.0, 1.0 - sum(pmf(k) for k in range(n)))
    return min(1.0, sum(pmf(k) for k in range(n + 1)))


def marginal_z(n, mu, sigma):
    """Signed significance of n observed vs mu expected with a Gaussian uncertainty sigma on mu
    (the Poisson tail averaged over mu ~ Normal(mu, sigma) truncated at 0; numerical quadrature)."""
    from statistics import NormalDist
    if mu <= 0:
        return float('nan')
    x = np.linspace(-6, 6, 481)
    pdf = np.exp(-0.5 * x ** 2)
    mus = mu + sigma * x
    ok = mus > 0
    up = n > mu                                               # the tail on the side of the observation
    p = np.sum(pdf[ok] * np.array([poisson_tail(n, m, up) for m in mus[ok]])) / np.sum(pdf[ok])
    p = min(max(p, 1e-300), 0.5)
    return float(np.sign(n - mu) * NormalDist().inv_cdf(1 - p))


SIGMA_LEVELS = {6: np.array([-1, 1, -2, 2, -3, 3], float), 7: np.linspace(-3, 3, 7)}
WEIGHT_CLIP = (0.0, 10.0)        # upstream's AR23+ template (general.weight_clip), applied again on read


def sigma_draws(npts, nuniv, rng):
    """The sigmas a multisigma knob (6 or 7 weight points) is evaluated at; None for multisim universes."""
    return rng.normal(0, 1, nuniv) if npts in SIGMA_LEVELS else None


def interpolate_universes(sums, s):
    """Universe counts of one knob from its selected-event sums per weight point (... x points).
    Multisigma knobs (GENIE -1,+1,-2,+2,-3,+3 or detector -3..+3) are linearly interpolated at the sigmas
    `s` (drawn once per knob, shared by all events and bins), clamped at the end points, as
    spineplot/systematic.py does per event; the interpolation is linear in the weights, so doing it on
    the sums is exact. Without `s` the points are the universes."""
    if s is None:
        return sums
    lv = SIGMA_LEVELS[sums.shape[-1]]
    order = np.argsort(lv)
    lv, sums = lv[order], sums[..., order]
    i = np.clip(np.searchsorted(lv, s) - 1, 0, len(lv) - 2)
    t = np.clip((s - lv[i]) / (lv[i + 1] - lv[i]), 0, 1)
    return sums[..., i] * (1 - t) + sums[..., i + 1] * t


def knob_universe_sums(W, S, nuniv, rng):
    """Selected-event sums per universe for one knob: S (stages x events, 0/1) @ weights (events x points)."""
    return interpolate_universes(S @ W, sigma_draws(W.shape[1], nuniv, rng))


def job_mask(kind, j, P):
    """Rows of a (rows x cuts) pass matrix in a cut-flow stage ('stage', j cuts applied), before all cuts,
    after all cuts, or with all cuts but cut j ('nm1')."""
    if kind == 'stage':
        return np.all(P[:, :j], axis=1)
    if kind == 'before':
        return np.ones(len(P), bool)
    if kind == 'after':
        return np.all(P, axis=1)
    return np.all(np.delete(P, j, axis=1), axis=1)


def bin_index(x, edges):
    """Histogram bin of each value as the plots fill it (values outside the range go to the edge bins,
    NaN to -1)."""
    x = np.asarray(x, float)
    lo, hi = edges[0], edges[-1]
    i = np.clip(np.searchsorted(edges, np.clip(x, lo, hi - 1e-9 * (hi - lo)), side='right') - 1, 0, len(edges) - 2)
    i[~np.isfinite(x)] = -1
    return i


def stage_jobs(cl):
    """One single-bin band job per cut-flow stage (preselection + each cut, cumulative)."""
    return [dict(key=('stage', k), kind='stage', j=k, var=None, edges=None) for k in range(len(cl) + 1)]


def sys_band(path, cl, region_fn, jobs, scale, clip=WEIGHT_CLIP, nuniv=1000, seed=12345, step=4000):
    """Flux + cross section uncertainty on the expected neutrino counts of a set of histograms.

    `path` is a run_systematics output: events/full/TREE plus TREE_{multisim,multisigma}Tree, filled row by
    row with it. Each job (key, kind, j, var, edges) is a histogram of `var` (one bin if None) over the
    rows job_mask(kind, j, cuts) selects within `region_fn`. For each knob and bin, the variance of `scale`
    x (universe count - CV count) over its universes; knobs add in quadrature (spineplot's total band).
    Weights are clamped to `clip`, NaN set to 1. The weights are read `step` rows at a time.
    Returns dict(cv={key: counts}, var={key: variance}, per_knob={knob: variance per stage job},
    points={knob: point sums / CV after all cuts}, n_knobs)."""
    import uproot
    import awkward as ak
    f = uproot.open(path)
    d = f['events/full']
    t = d[TREE]
    names = [k.split(';')[0] for k in d.keys()]
    cols = sorted({c['var'] for c in cl} | {NSH} | {jb['var'] for jb in jobs if jb['var']})
    df = t.arrays([c for c in cols if c in t.keys()] + [c for c in TRIG if c in t.keys()], library='pd')
    trig = next((c for c in TRIG if c in df), None)
    keep = np.asarray(region_fn(df), bool)
    if trig:
        keep &= df[trig].to_numpy(dtype=float) == 1
    P = np.column_stack([c['fn'](df) for c in cl])
    idx, nb = {}, {}
    missing = sorted({jb['var'] for jb in jobs if jb['var'] and jb['var'] not in df})
    if missing:
        print(f'  NOTE: {path} lacks {missing}: those histograms get the MC-statistics band only')
    for jb in jobs:
        if jb['var'] in missing:
            continue
        m = job_mask(jb['kind'], jb['j'], P) & keep
        k = jb['key']
        if jb['var'] is None:
            idx[k], nb[k] = np.where(m, 0, -1), 1
        else:
            idx[k], nb[k] = np.where(m, bin_index(df[jb['var']].to_numpy(), jb['edges']), -1), len(jb['edges']) - 1
    cv = {k: np.bincount(i[i >= 0], minlength=nb[k]).astype(float) for k, i in idx.items()}
    stages = [jb['key'] for jb in jobs if jb['kind'] == 'stage']
    final = max(stages, key=lambda k: k[1]) if stages else None
    rng = np.random.default_rng(seed)
    var = {k: np.zeros(nb[k]) for k in idx}
    per_knob, points, clipped = {}, {}, {}
    for suffix in ('multisimTree', 'multisigmaTree'):
        tname = f'{TREE}_{suffix}'
        if tname not in names:
            print(f'  {path}: no events/full/{tname}')
            continue
        wt = d[tname]
        if wt.num_entries != t.num_entries:
            raise SystemExit(f'{tname} has {wt.num_entries} rows, {TREE} {t.num_entries}: not row-aligned')
        for knob in wt.keys():
            if knob in ('Run', 'Subrun', 'Evt') or knob.endswith('_sigma'):
                continue
            sums, bad = {}, None
            for a in range(0, t.num_entries, step):
                b = min(t.num_entries, a + step)
                kc = keep[a:b]
                if not kc.any():
                    continue
                arr = wt[knob].array(entry_start=a, entry_stop=b, library='ak')
                if arr.ndim < 2:
                    bad = 'scalar'
                    break
                arr = arr[kc]
                lens = np.unique(ak.to_numpy(ak.num(arr, axis=1)))
                npts = next(iter(sums.values())).shape[1] if sums else None
                if len(lens) != 1 or lens[0] == 0 or (npts is not None and lens[0] != npts):
                    bad = f'ragged/empty weight vectors (lengths {lens[:5]})'
                    break
                W = ak.to_numpy(ak.to_regular(arr, axis=1)).astype(float)
                fin = np.isfinite(W)
                W = np.where(fin, W, 1.0)
                out_ = ((W < clip[0]) | (W > clip[1])).any(axis=1) | ~fin.all(axis=1)
                if out_.any():                                # what the clamp touches (reported, not hidden)
                    cr = clipped.setdefault(knob, [0, 0, -np.inf])
                    cr[0] += int(out_.sum())
                    if final is not None:
                        cr[1] += int((out_ & (idx[final][a:b][kc] >= 0)).sum())
                    cr[2] = max(cr[2], float(W.max()))
                W = np.clip(W, *clip)
                for k, i_all in idx.items():
                    ii = i_all[a:b][kc]
                    sel = ii >= 0
                    if not sel.any():
                        continue
                    order = np.argsort(ii[sel], kind='stable')
                    bins, starts = np.unique(ii[sel][order], return_index=True)
                    if k not in sums:
                        sums[k] = np.zeros((nb[k], W.shape[1]))
                    sums[k][bins] += np.add.reduceat(W[sel][order], starts, axis=0)
            if bad == 'scalar' or not sums:
                continue
            if bad:
                print(f'  WARNING {knob}: {bad}: skipped')
                continue
            npts = next(iter(sums.values())).shape[1]
            s = sigma_draws(npts, nuniv, rng)
            kv = {}
            for k in idx:
                N = interpolate_universes(sums.get(k, np.zeros((nb[k], npts))), s)
                kv[k] = (scale ** 2) * np.mean((N - cv[k][:, None]) ** 2, axis=1)
                var[k] += kv[k]
            per_knob[knob] = np.array([kv[k][0] for k in stages])
            if final is not None and cv[final][0] > 0:
                points[knob] = sums.get(final, np.zeros((1, npts)))[0] / cv[final][0]
    return dict(cv={k: scale * v for k, v in cv.items()}, var=var, per_knob=per_knob, points=points,
                n_knobs=len(per_knob), clipped=clipped)


def build(a, want):
    import pandas as pd
    key = 'full' if has(a.sys, 'full') else 'cv'
    opt = PI0 + list(TRIG)
    full, pot_cv, _ = load(a.sys, key, want, opt)
    off, _, lt_off = load(a.sys, 'offbeam', want, opt)
    on, pot_on, lt_on = load(a.sys, 'onbeam', want, opt)
    legacy = not all(c in full for c in PI0)
    if legacy:
        print(f'NOTE: no pi0 branches ({PI0[0]}, ...): secondary showers split by start point only')
    for i, lab in (LEGACY_LABELS if legacy else PI0_LABELS).items():
        GROUPS[i] = (lab, GROUPS[i][1])
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
    note = ('dropped (--no-overlay-cosmics): off-beam is the only cosmic source' if a.no_overlay_cosmics else
            'kept as "cosmics in nu readouts"')
    print(f'events/{key}: {len(full)} of {n0} rows have a neutrino match; the other {len(overlay)} '
          f'(overlay cosmics, {len(overlay) * pot_on / pot_cv:.1f} at on-beam POT) {note}; '
          f'POT {pot_cv:.4g}; off-beam livetime {lt_off:.4g}, on-beam {lt_on:.4g}')
    full = triggered(full, f'events/{key} (nu matched)', a.no_trigger_emulation)
    overlay = triggered(overlay, f'events/{key} (overlay cosmics)', a.no_trigger_emulation)
    full['g'], full['w'], full['src'] = group_of(full), pot_on / pot_cv, 'nu'
    off['g'], off['w'], off['src'] = COSMIC, lt_on / lt_off, 'offbeam'
    parts = [full, off]
    if not a.no_overlay_cosmics:
        if key == 'full':
            raise SystemExit('a sys ROOT keeps the overlay cosmics in *_nonmatched: pass a selection output '
                             '(events/cv) or --no-overlay-cosmics')
        overlay = overlay.copy()
        overlay['g'], overlay['w'], overlay['src'] = COSMIC, pot_on / pot_cv, 'overlay'
        parts.append(overlay)
    dirt = a.dirt or (a.sys if has(a.sys, 'dirt') else None)
    if dirt:
        dr, pot_dr, _ = load(dirt, 'dirt', want, opt)
        dr = triggered(matched(dr), 'events/dirt', a.no_trigger_emulation)
        dr['g'], dr['w'], dr['src'] = DIRT, pot_on / pot_dr, 'dirt'
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
    ap.add_argument('--no-overlay-cosmics', action='store_true',
                    help='drop the cv rows without a neutrino match (cosmics recorded in neutrino-triggered '
                         'readouts); off-beam then is the only cosmic source (the old matched-only scheme)')
    ap.add_argument('--no-trigger-emulation', action='store_true',
                    help='keep MC rows of readouts that fail the trigger emulation (data have them only when a '
                         'cosmic fired the trigger, which off-beam already counts)')
    ap.add_argument('--exposure-from', help='ROOT file to take missing POT/Livetime histograms from (same onbeam/'
                                             'offbeam samples; only needed for old sys ROOTs)')
    ap.add_argument('--weights', help='run_systematics output for this selection (events/full/TREE + its '
                                      'multisim/multisigma trees): adds the flux + xsec band and MC statistics '
                                      'to the cut flow, with the data/MC significance (sysband.csv)')
    ap.add_argument('--weight-clip', type=float, nargs=2, default=list(WEIGHT_CLIP), metavar=('LO', 'HI'),
                    help='clamp every systematic weight to [LO, HI] on read (default: upstream general.weight_clip)')
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
    cos = mc['g'].to_numpy() == COSMIC
    src = mc['src'].to_numpy()
    print('  cosmics: ' + ', '.join(f'{lab} {mc["w"].to_numpy()[cos & (src == s)].sum():.1f}'
                                    for s, lab in (('offbeam', 'off-beam data'), ('overlay', 'overlay (nu readouts)'))))

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
    print(f'  {"stage":<38} {"MC":>8} {"data":>6} {"z":>5} |{hdr} | {"1g kept":>8} {"pi0 in":>8} {"pi0 out":>8} '
          f'{"oth sec":>8} {"1g/pi0":>7}')
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
        kin = r[plain(GROUPS[SEC_IN][0])] / tot[SEC_IN] if tot[SEC_IN] else float('nan')
        kout = r[plain(GROUPS[SEC_OUT][0])] / tot[SEC_OUT] if tot[SEC_OUT] else float('nan')
        koth = r[plain(GROUPS[SEC_OTHER][0])] / tot[SEC_OTHER] if tot[SEC_OTHER] else float('nan')
        tsec = tot[SEC_IN] + tot[SEC_OUT]
        ks = (r[plain(GROUPS[SEC_IN][0])] + r[plain(GROUPS[SEC_OUT][0])]) / tsec if tsec else float('nan')
        cols = ''.join(f'{r[plain(GROUPS[i][0])]:9.1f}' for i in SHOWN)
        print(f'  {r["stage"]:<38} {r["mc_expected"]:8.1f} {r["data_n"]:6d} {r["z"]:5.1f} |{cols} | '
              f'{kg:8.3f} {kin:8.4f} {kout:8.4f} {koth:8.4f} {kg / ks if ks else float("inf"):7.1f}')

    print('\nN-1: all other cuts applied; what this cut alone then removes (fraction of each group it rejects)')
    print(f'  {"cut left out":<38} {"MC":>8} {"data":>6} {"z":>5} | {"rejects: 1g":>11} {"pi0 in":>6} {"pi0 out":>7} '
          f'{"oth sec":>7} {"cosmic":>7} {"other nu":>8} {"MC all":>7} {"data":>6}')
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
              f'{rej(g == SEC_IN):6.3f} {rej(g == SEC_OUT):7.3f} {rej(g == SEC_OTHER):7.3f} {rej(g == COSMIC):7.3f} {rej(g == OTHER):8.3f} '
              f'{rej(np.ones(len(mc), bool)):7.3f} {rd:6.3f}')
    pd.DataFrame(rows).to_csv(os.path.join(a.out, 'cutflow.csv'), index=False)
    specs = plot_specs(cl)
    band = None
    if a.weights:
        nu = mc['src'].to_numpy() == 'nu'
        scale = float(mc['w'].to_numpy()[nu][0]) if nu.any() else 1.0
        jobs = stage_jobs(cl) + [dict(key=(s['sub'], s['c']['var']), kind=s['kind'], j=s['j'], var=s['c']['var'],
                                      edges=s['edges']) for s in specs]
        band = sys_band(a.weights, cl, lambda d: d[NSH].to_numpy() == a.n_showers, jobs, scale,
                        clip=tuple(a.weight_clip))
        uncertainty_table(a, mc, on, p_mc, p_on, cl, stages, band)

    plots(mc, on, p_mc, p_on, specs, a.out, pot_on, region, band)
    print(f'\nwrote {a.out}/cutflow.csv and {a.out}/{{before,after,nm1}}/*.png,pdf')
    return 0


def uncertainty_table(a, mc, on, p_mc, p_on, cl, stages, band):
    """Cut flow with its uncertainties: MC statistics, flux + cross section (from the --weights sys ROOT),
    and the data/MC significance with and without them. Writes sysband.csv."""
    import pandas as pd
    src, w = mc['src'].to_numpy(), mc['w'].to_numpy()
    nu = src == 'nu'
    cv = np.array([band['cv'][('stage', k)][0] for k in range(len(stages))])
    var = np.array([band['var'][('stage', k)][0] for k in range(len(stages))])
    m_mc = np.vstack([job_mask('stage', k, p_mc) for k in range(len(stages))])
    m_on = np.vstack([job_mask('stage', k, p_on) for k in range(len(stages))])
    nu_sel = np.array([w[m & nu].sum() for m in m_mc])
    print(f'\nuncertainties: flux + cross section from {a.weights} ({band["n_knobs"]} knobs, weights clamped to '
          f'{list(a.weight_clip)}; no detector variations in this production), MC statistics sqrt(sum w^2) over '
          f'all MC terms')
    lost = np.abs(cv - nu_sel) > 0.01 * np.maximum(nu_sel, 1e-9) + 0.05
    if lost.any():
        k = int(np.argmax(lost))
        print(f'  NOTE: at "{stages[k][0]}" the weighted tree holds {cv[k]:.2f} nu events, the selection file '
              f'{nu_sel[k]:.2f}: the band is taken as a fraction from the weighted rows (the difference is rows '
              f'run_systematics could not match to a neutrino in the CAFs, or a different selection file)')
    print(f'  {"stage":<38} {"MC":>8} {"nu":>7} {"stat":>6} {"syst":>6} {"syst%":>6} {"data":>6} '
          f'{"z stat":>7} {"z all":>7}')
    rows = []
    for k, (name, j) in enumerate(stages):
        mu = w[m_mc[k]].sum()
        st = np.sqrt((w[m_mc[k]] ** 2).sum())
        sy = np.sqrt(var[k]) * (nu_sel[k] / cv[k] if cv[k] > 0 else 1.0)   # fractional band on the nu part
        n = int(m_on[k].sum())
        z0 = marginal_z(n, mu, st)
        z1 = marginal_z(n, mu, np.hypot(st, sy))
        label = name if j is None else '+ ' + name
        rows.append(dict(stage=label, mc_expected=mu, nu_expected=nu_sel[k], mc_stat=st, syst=sy,
                         syst_frac=sy / nu_sel[k] if nu_sel[k] else float('nan'), data_n=n,
                         z_stat=z0, z_stat_syst=z1))
        print(f'  {label:<38} {mu:8.2f} {nu_sel[k]:7.2f} {st:6.2f} {sy:6.2f} '
              f'{100 * rows[-1]["syst_frac"]:5.0f}% {n:6d} {z0:7.2f} {z1:7.2f}')
    last = len(stages) - 1
    top = sorted(band['per_knob'].items(), key=lambda kv: kv[1][last], reverse=True)[:8]
    if top and var[last] > 0:
        print('  largest knobs after all cuts: share of the flux + xsec variance, and the selected count in the '
              'knob\'s variations relative to the CV')
        for knob, v in top:
            p = band['points'].get(knob)
            if p is None:
                how = ''
            elif len(p) == 6:
                how = f'-1 sigma x{p[0]:.2f}, +1 sigma x{p[1]:.2f}'
            elif len(p) == 7:
                how = f'-1 sigma x{p[2]:.2f}, +1 sigma x{p[4]:.2f}'
            else:
                how = f'{len(p)} universes: mean x{p.mean():.2f}, sd {p.std():.2f}'
            print(f'    {knob:<60} {100 * v[last] / var[last]:5.1f}%   {how}')
    cr = band.get('clipped', {})
    nrow = int(round(cv[0] / max(w[nu][0], 1e-300))) if nu.any() else 0
    if cr:
        print(f'  weight clamp {list(a.weight_clip)}: {len(cr)} knobs have events with a weight outside it (or NaN); '
              f'events of the {nrow} in the region, of those after all cuts, largest raw weight:')
        for knob, (n_reg, n_fin, mx) in sorted(cr.items(), key=lambda kv: -kv[1][0])[:10]:
            print(f'    {knob:<60} {n_reg:6d} {n_fin:4d}   {mx:.3g}')
    else:
        print(f'  weight clamp {list(a.weight_clip)}: no weight outside it (the clamp changes nothing)')
    print('  z stat: Poisson for the data with the MC-statistics uncertainty on the prediction; z all: also the '
          'flux + xsec band (both Gaussian on the prediction, marginalised)')
    pd.DataFrame(rows).to_csv(os.path.join(a.out, 'sysband.csv'), index=False)


KE_SPEC = dict(var=KE, range=(0.0, 2000.0), nbins=20, xlabel='Shower KE [MeV]', lines=[], name='shower KE')


def plot_specs(cl):
    """Every figure: each cut variable (and the shower KE) before and after all cuts, and each cut variable
    with all the other cuts (N-1), with the bin edges the cut values fall on."""
    specs = []
    for c in cl + [KE_SPEC]:
        edges = aligned_edges(*c['range'], c['nbins'], [x for x, _ in c.get('lines', [])])
        specs.append(dict(sub='before', c=c, kind='before', j=None, edges=edges, stage='preselection'))
        specs.append(dict(sub='after', c=c, kind='after', j=None, edges=edges, stage='all cuts'))
    for j, c in enumerate(cl):
        edges = aligned_edges(*c['range'], c['nbins'], [x for x, _ in c['lines']])
        specs.append(dict(sub='nm1', c=c, kind='nm1', j=j, edges=edges, stage=f'all cuts but {c["name"]}'))
    return specs


def plots(mc, on, p_mc, p_on, specs, out, pot_on, region, band=None):
    """The stacked truth groups with the data; a hatched band on the MC total: MC statistics, plus the
    flux + xsec band of the neutrino MC when `band` (sys_band) is given."""
    try:
        import matplotlib
    except ImportError:
        print('no matplotlib: plots skipped')
        return
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.transforms as mtransforms
    import warnings
    warnings.filterwarnings('ignore', message='The figure layout has changed to tight')
    here = os.path.dirname(os.path.abspath(__file__))
    sheet = os.path.join(here, '..', 'spineplot', 'styles', 'muon2024.mplstyle')
    if os.path.exists(sheet):
        plt.style.use(sheet)
    w_all = mc['w'].to_numpy(dtype=float)
    for sp in specs:
        c, edges = sp['c'], sp['edges']
        m_mc, m_on = job_mask(sp['kind'], sp['j'], p_mc), job_mask(sp['kind'], sp['j'], p_on)
        lines = c.get('lines', [])
        d = os.path.join(out, sp['sub'])
        os.makedirs(d, exist_ok=True)
        fig, ax = plt.subplots(figsize=(10, 7.5))
        lo, hi = edges[0], edges[-1]
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
        # uncertainty band on the MC total
        bi = bin_index(sel[c['var']].to_numpy(dtype=float), edges)
        ok = bi >= 0
        nb = len(edges) - 1
        tot = np.bincount(bi[ok], weights=w_all[m_mc][ok], minlength=nb)
        var = np.bincount(bi[ok], weights=w_all[m_mc][ok] ** 2, minlength=nb)
        what = 'MC stat.'
        key = (sp['sub'], c['var'])
        if band is not None and key in band['var']:
            var = var + band['var'][key]
            what = 'MC stat. + flux/xsec syst.'
        err = np.sqrt(var)
        if tot.sum() > 0:
            ax.bar(0.5 * (edges[1:] + edges[:-1]), 2 * err, width=np.diff(edges), bottom=tot - err, color='none',
                   edgecolor='dimgray', hatch='////', linewidth=0, label=what)
        x_on = np.clip(dsel[c['var']].to_numpy(dtype=float), lo, hi - 1e-9 * (hi - lo))
        k = np.histogram(x_on, edges)[0]
        ctr = 0.5 * (edges[1:] + edges[:-1])
        ax.errorbar(ctr, k, yerr=np.sqrt(k), fmt='o', color='k', ms=4, label=f'Data (BNB) ({k.sum()})')
        ax.set_xlim(lo, min(hi, c['range'][1]) if c['range'][1] > lo else hi)  # don't show past e.g. softmax = 1
        ax.set_ylim(0, 1.5 * max(ax.get_ylim()[1], 1.0))
        trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        span = 0.06 * (hi - lo)
        for n, (x, side) in enumerate(lines):
            ax.axvline(x, color='r', ls='--', label='cut (arrow: kept)' if n == 0 else None)
            x0, x1 = ax.get_xlim()
            room = (x1 - x) if side > 0 else (x - x0)
            ax.annotate('', xy=(x + side * min(span, 0.9 * room), 0.6), xytext=(x, 0.6), xycoords=trans,
                        textcoords=trans,
                        arrowprops=dict(arrowstyle='->', color='r', lw=2))
        ax.set_xlabel(c['xlabel'])
        ax.set_ylabel('Candidates')
        ax.text(0.0, 1.01, PRELIM, transform=ax.transAxes, ha='left', va='bottom', color='red',
                fontweight='bold', fontsize='medium')
        ax.text(1.0, 1.01, f'{pot_on:.2e} POT', transform=ax.transAxes, ha='right', va='bottom', fontsize='medium')
        ax.legend(loc='upper right', fontsize='x-small', frameon=False, title=f'{region}\n{sp["stage"]}',
                  title_fontsize='x-small')
        fig.tight_layout()
        for ext in ('png', 'pdf'):
            fig.savefig(os.path.join(d, f'{c["var"]}.{ext}'))
        plt.close(fig)


if __name__ == '__main__':
    sys.exit(main())
