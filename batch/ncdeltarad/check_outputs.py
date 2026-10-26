#!/usr/bin/env python3
"""
End-to-end check of a finished ncdeltarad campaign, from what it left on /pnfs.
Run on a gpvm with uproot + awkward (e.g. the medulla conda env):

    python batch/ncdeltarad/check_outputs.py /pnfs/icarus/persistent/users/$USER/ncdeltarad/<name>

FAIL lines set a non-zero exit; WARN lines need a look; the rest is a census.

  1. manifest: every stage good == expected; every published list and the
     files it names exist.
  2. CAF directory census, per job: which files are flat CAFs (recTree with
     plain rec.hdr.evt leaves), nested CAFs (recTree of caf::StandardRecord),
     TFileService files (no recTree) -- and that every job has a CAF.
  3. CAF truth (flat CAF preferred): every true neutrino NC and RES; the CAF
     events join to the gen sidecar on (run, subrun, event); the sidecar's
     radiative fraction is enhanced (~0.4) and the weighted one back at the
     nominal value; a true primary photon in every radiative event.
     Also: events with >= 1 slice; reweight groups per
     neutrino the same everywhere (see batch/diag_caf_weights.py).
  4. larcv (stage1, the SPINE input): files open and hold as many entries as
     the CAFs hold events.
"""
import argparse
import csv
import glob
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np

KEYS = ['rec.hdr.run', 'rec.hdr.subrun', 'rec.hdr.evt']
TRUTH = ['rec.mc.nu.iscc', 'rec.mc.nu.genie_mode', 'rec.mc.nu.pdg', 'rec.mc.nu.prim.pdg',
         'rec.mc.nu.prim.startE', 'rec.mc.nu.wgt..length']
GENIE_RES = 1                     # simb / caf genie_mode: kQE 0, kRes 1, kDIS 2, kCoh 3
fails, warns = [], []


def fail(msg):
    fails.append(msg)
    print(f'FAIL  {msg}')


def warn(msg):
    warns.append(msg)
    print(f'WARN  {msg}')


def info(msg):
    print(f'      {msg}')


def check_manifest(base):
    print('\n== 1. manifest and published lists')
    path = os.path.join(base, 'manifest.json')
    if not os.path.exists(path):
        fail(f'{path} missing (campaign not finished?)')
        return {}
    m = json.load(open(path))
    info(f"{m['name']}: {m['release']} {m['qual']}, BR(Ngamma)={m['br_gamma']}, finished {m.get('finished')}")
    for name, s in m['stages'].items():
        good, exp = s.get('good'), s.get('expected')
        line = f'{name:8s} good {good}/{exp}  {s["outdir"]}'
        if good is None or exp is None:
            fail(f'{line}: no counts')
        elif good > exp:
            warn(f'{line}: more good outputs than jobs -- the stage counts files, not jobs')
        elif good < exp:
            warn(f'{line}: {exp - good} jobs lost')
        else:
            info(line)
        for lname, entry in (s.get('lists') or {}).items():
            lpath = entry['list']                # orchestrator writes {list: path, n: count}
            if not os.path.exists(lpath):
                fail(f'{name}/{lname}: list {lpath} missing')
                continue
            files = [x.strip() for x in open(lpath) if x.strip()]
            missing = [f for f in files if not os.path.exists(f)]
            (fail if missing or not files else info)(
                f'{name}/{lname}: {len(files)} files listed, {len(missing)} missing'
                + (f' e.g. {missing[0]}' if missing else '') + ('' if files else ' -- EMPTY list'))
    return m


def classify(path, uproot):
    """'flat' | 'nested' | 'hist' | 'unreadable: ...' for one .root file."""
    try:
        f = uproot.open(path)
    except Exception as e:                       # noqa: BLE001
        return f'unreadable: {e}', 0
    if 'recTree' not in f:
        return 'hist', 0
    t = f['recTree']
    keys = set(t.keys())
    if 'rec.hdr.evt' in keys and 'rec' not in keys:
        return 'flat', t.num_entries
    return 'nested', t.num_entries


def census(base, uproot):
    print('\n== 2. CAF directory census')
    jobs = sorted(d for d in glob.glob(os.path.join(base, 'caf', 'out', '*')) if os.path.isdir(d))
    if not jobs:
        fail(f'no job directories under {base}/caf/out')
    by_kind = defaultdict(list)
    for d in jobs:
        kinds = []
        for f in sorted(glob.glob(os.path.join(d, '*.root'))):
            k, n = classify(f, uproot)
            by_kind[k.split(':')[0]].append((f, n))
            kinds.append(k)
            info(f'{os.path.basename(d)}  {k:7s} {n:6d} ev  {os.path.getsize(f) / 1e6:7.1f} MB  {os.path.basename(f)}')
        if not any(k in ('flat', 'nested') for k in kinds):
            fail(f'{os.path.basename(d)}: no CAF in the job directory')
    for k, v in sorted(by_kind.items()):
        info(f'{k}: {len(v)} files')
    if not by_kind['flat']:
        warn('no flat CAF found (CAFMaker CreateFlatCAF defaults to true): check the job log for the '
             'flat file and how it was renamed')
    for f, _ in by_kind['unreadable']:
        fail(f'unreadable: {f}')
    return by_kind


def load_caf(files, uproot, ak):
    """Event keys, truth and reco arrays from CAFs (flat preferred)."""
    parts = []
    for f in files:
        t = uproot.open(f)['recTree']
        have = set(t.keys())
        want = [b for b in KEYS + TRUTH + ['rec.nslc'] if b in have]
        parts.append(t.arrays(want, library='ak'))
    return ak.concatenate(parts) if parts else None


def check_truth(m, by_kind, uproot, ak, base):
    print('\n== 3. CAF truth, reco and weights vs the gen sidecar')
    files = [f for f, _ in by_kind['flat']] or [f for f, _ in by_kind['nested']]
    if not files:
        return None
    try:
        a = load_caf(files, uproot, ak)
    except Exception as e:                       # noqa: BLE001
        fail(f'reading {"flat" if by_kind["flat"] else "nested"} CAF branches: {e}')
        return None
    if a is None or not all(k in a.fields for k in KEYS):
        fail(f'CAF lacks {KEYS} -- nested CAF not readable by uproot; need the flat CAF')
        return None
    nev = len(a)
    info(f'{nev} CAF events in {len(files)} files')
    if 'rec.mc.nu.iscc' in a.fields:
        nnu = ak.num(a['rec.mc.nu.iscc'])
        info(f'true neutrinos per event: {dict(sorted(Counter(nnu.tolist()).items()))}')
        cc = ak.flatten(a['rec.mc.nu.iscc'])
        (fail if ak.sum(cc) else info)(f'{int(ak.sum(cc))}/{len(cc)} true neutrinos are CC (NCRES only -> 0)')
        if 'rec.mc.nu.genie_mode' in a.fields:
            mode = Counter(ak.flatten(a['rec.mc.nu.genie_mode']).tolist())
            (info if set(mode) <= {GENIE_RES} else fail)(f'genie_mode counts {dict(mode)} (RES = {GENIE_RES})')
        if 'rec.mc.nu.wgt..length' in a.fields:
            ng = Counter(ak.flatten(a['rec.mc.nu.wgt..length']).tolist())
            (info if len(ng) == 1 else warn)(f'reweight groups per neutrino: {dict(ng)}')
    if 'rec.nslc' in a.fields:
        ns = np.asarray(a['rec.nslc'])
        info(f'events with >= 1 slice: {np.mean(ns > 0):.3f}; mean slices/event {ns.mean():.2f}')
        if np.mean(ns > 0) < 0.5:
            warn('fewer than half the events have a reco slice')

    sc = (m.get('sidecar') or {}).get('weights') or os.path.join(base, 'gen_sidecar', 'ncdeltarad_weights.csv')
    if not os.path.exists(sc):
        fail(f'sidecar weights {sc} missing')
        return a
    rows = list(csv.DictReader(open(sc)))
    side = {(int(r['run']), int(r['subrun']), int(r['event'])): r for r in rows}
    if len(side) != len(rows):
        fail(f'sidecar keys not unique: {len(rows)} rows, {len(side)} keys')
    keys = list(zip(*(np.asarray(a[k]).tolist() for k in KEYS)))
    if len(set(keys)) != len(keys):
        fail(f'CAF (run, subrun, evt) not unique: {len(keys)} events, {len(set(keys))} keys')
    hit = [k for k in keys if k in side]
    info(f'sidecar: {len(side)} NC Delta events; {len(hit)}/{len(keys)} CAF events join to it, '
         f'{len(set(side) - set(keys))} sidecar events not in the CAFs (filtered or lost jobs)')
    miss = [k for k in keys if k not in side]
    if miss:
        # An NCRES event without a Δ+/Δ0 (a heavier resonance decaying
        # straight to Nπ) is legitimately absent; a whole (run, subrun) absent
        # from the sidecar points at a gen file the sidecar never read.
        grp = lambda ks: dict(sorted(Counter((k[0], k[1]) for k in ks).items()))
        info(f'CAF events not in the sidecar, per (run, subrun): {grp(miss)}')
        info(f'sidecar events per (run, subrun): {grp(side)}')
        absent = set((k[0], k[1]) for k in miss) - set((k[0], k[1]) for k in side)
        if absent:
            warn(f'(run, subrun) {sorted(absent)} has CAF events but no sidecar row at all: '
                 f'was its gen file in the sidecar dump? (check gen_sidecar/pilot_check.txt)')
    if not hit:
        fail('no CAF event joins the sidecar on (run, subrun, event)')
        return a
    rad = np.array([int(side[k]['radiative']) for k in hit])
    w = np.array([float(side[k]['weight']) for k in hit])
    info(f'joined events: radiative fraction {rad.mean():.3f} raw, {np.sum(w * rad) / np.sum(w):.4f} weighted; '
         f'mean weight {w.mean():.3f}')
    if not 0.2 < rad.mean() < 0.6:
        warn('raw radiative fraction outside 0.2-0.6: was the patched table used?')
    if np.sum(w * rad) / np.sum(w) > 0.02:
        warn('weighted radiative fraction above 2%: weights do not restore the nominal BR')
    if 'rec.mc.nu.prim.pdg' in a.fields:
        idx = {k: i for i, k in enumerate(keys)}
        prim = a['rec.mc.nu.prim.pdg']
        has_g = np.array([bool(ak.any(ak.flatten(prim[idx[k]], axis=None) == 22)) for k in hit])
        r = rad.astype(bool)
        # Every radiative Δ decay leaves a primary photon. The converse need not
        # hold: nuclear de-excitation and heavier-resonance photons are primaries
        # too (pilot: 4 of 50).
        (fail if np.any(r & ~has_g) else info)(
            f'radiative events without a true primary photon: {int(np.sum(r & ~has_g))}/{int(r.sum())}; '
            f'non-radiative events with one (de-excitation etc.): {int(np.sum(has_g & ~r))}/{int((~r).sum())}')
    return a


def check_larcv(m, caf_events, uproot):
    print('\n== 4. larcv (SPINE input)')
    lst = (((m.get('stages') or {}).get('stage1', {}).get('lists') or {}).get('larcv') or {}).get('list')
    if not lst or not os.path.exists(lst):
        fail(f'stage1 larcv list missing ({lst})')
        return
    files = [x.strip() for x in open(lst) if x.strip()]
    tot = 0
    for f in files:
        try:
            u = uproot.open(f)
            trees = sorted(k.split(';')[0] for k, c in u.classnames().items() if c == 'TTree')
            n = {u[t].num_entries for t in trees}
            if len(n) != 1:
                warn(f'{os.path.basename(f)}: trees disagree on entries {sorted(n)}')
            tot += max(n) if n else 0
            info(f'{os.path.basename(f)}: {len(trees)} trees, {sorted(n)} entries')
        except Exception as e:                   # noqa: BLE001
            fail(f'{f}: {e}')
    if caf_events is not None:
        (info if tot == caf_events else warn)(f'larcv entries {tot} vs CAF events {caf_events}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('base', help='persistent_base of the campaign (holds manifest.json, caf/, stage1/, gen_sidecar/)')
    a = ap.parse_args()
    try:
        import awkward as ak
        import uproot
    except ImportError as e:
        raise SystemExit(f'needs uproot + awkward ({e})')
    m = check_manifest(a.base)
    by_kind = census(a.base, uproot)
    arr = check_truth(m, by_kind, uproot, ak, a.base) if m else None
    check_larcv(m, None if arr is None else len(arr), uproot)
    print(f'\n{len(fails)} FAIL, {len(warns)} WARN')
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main())
