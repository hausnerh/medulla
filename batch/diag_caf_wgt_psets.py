#!/usr/bin/env python3
"""What does each reweight group of a flat CAF actually hold, and does a systematics toml read it right?

run_systematics takes knob `index` i from rec.mc.nu.wgt[i] and stores every universe it finds there; the
knob `name`, `type` and `nsigma` in the toml are labels only. This reads the CAF's own description of its
weight groups (globalTree: global.wgts.*: name, type, number of universes, the GENIE parameters each group
varies and, for a multisigma group, the sigma of each entry) and, with --toml, compares it with the
toml's knob list:
  * name mismatch at an index      -> the toml's indices do not fit this production (wrong labels)
  * universes vs the toml's nsigma -> a "multisigma" group with N entries is a sigma grid (see the
                                      printed sigmas), not N random universes
  * parameters shared by a multisim group and multisigma knobs that are both in the toml
                                   -> summing every knob counts those parameters twice

Read-only. Usage (one CAF from the weights list is enough):
  diag_caf_wgt_psets.py "$(head -1 build/weights_gOre_1g1p_sidebands_osc_run2.txt)" \\
      --toml systematics/toml/gOre_1g1p_sidebands_osc_run2_origin.toml \\
      [--sys build/output_gOre_1g1p_sidebands_osc_run2_origin_sys.root [--knob NAME ...]]
Reads flat ('global.wgts.name') and split ('global/wgts/wgts.name') CAF layouts.
"""
import argparse
import sys

import numpy as np


def split(flat, lengths):
    out, k = [], 0
    for n in lengths:
        out.append(flat[k:k + n])
        k += n
    return out


def canon(key):
    """Flat CAFs name branches 'global.wgts.name'; split (structured) ones 'global/wgts/wgts.name'."""
    leaf = key.split('/')[-1]
    return 'global.' + leaf if leaf.startswith('wgts') else leaf


def field(rec, sub):
    """The first record field whose last dotted component contains `sub` (split objects prefix their fields)."""
    for f in rec.fields:
        if sub in f.split('.')[-1]:
            return f
    return None


def parse_map(m):
    """global.wgts.map read as one object branch: per group, a list of {param: {name, ...}, vals: [...]}."""
    pf, vf = field(m, 'param'), field(m, 'vals')
    if pf is None:
        raise ValueError(f'no param field in {m.fields}')
    prm = m[pf]
    nf = field(prm, 'name') if prm.fields else None
    names = prm[nf] if nf else prm
    out = []
    for g in range(len(m)):
        out.append(([str(x) for x in names[g].tolist()],
                    [list(map(float, v)) for v in m[vf][g].tolist()] if vf else None))
    return out


def read_psets(path):
    import uproot
    f = uproot.open(path)
    gt = next((f[k] for k in f.keys() if k.split(';')[0].lower() == 'globaltree'), None)
    if gt is None:
        raise SystemExit(f'{path}: no globalTree; top-level keys: {[k.split(";")[0] for k in f.keys()]}')
    keys = {canon(k): k for k in gt.keys()}
    wk = [k for k in gt.keys() if 'wgts' in k]
    if not all(k in keys for k in ('global.wgts.name', 'global.wgts.nuniv')):
        print('globalTree weight branches found:\n  ' + '\n  '.join(wk))
        raise SystemExit('unexpected globalTree layout: paste the list above')

    def arr(k):
        return gt[keys[k]].array(library='np')[0] if k in keys else None
    names = [str(x) for x in arr('global.wgts.name')]
    psets = [dict(i=i, name=n) for i, n in enumerate(names)]
    for k, fld in (('global.wgts.type', 'type'), ('global.wgts.nuniv', 'nuniv')):
        v = arr(k)
        if v is not None:
            for p, x in zip(psets, v):
                p[fld] = int(x)
    nmap = arr('global.wgts.map..length')
    pnames = arr('global.wgts.map.param.name')
    vlen = arr('global.wgts.map.vals..length')
    vals = arr('global.wgts.map.vals')
    if nmap is not None and pnames is not None:                       # flat CAF
        per = split([str(x) for x in pnames], nmap)
        pv = split(list(split(vals, vlen)), nmap) if vals is not None and vlen is not None else [None] * len(per)
        for p, ps, vv in zip(psets, per, pv):
            p['params'] = ps
            p['vals'] = vv
    elif 'global.wgts.map' in keys:                                   # structured CAF: one object branch
        br = gt[keys['global.wgts.map']]
        try:
            m = br.array(library='ak')[0]
            for p, (ps, vv) in zip(psets, parse_map(m)):
                p['params'] = ps
                p['vals'] = vv
        except Exception as e:                                        # report it, keep the name/type/nuniv table
            print(f'NOTE: could not read {keys["global.wgts.map"]} ({type(e).__name__}: {e});\n'
                  f'  typename {br.typename!r}, interpretation {br.interpretation!r};\n'
                  f'  sub-branches {list(br.keys())}: parameters/sigmas unavailable')
    else:
        print('NOTE: no global.wgts.map branches; parameters/sigmas unavailable. Weight branches:\n  ' + '\n  '.join(wk))
    return psets


def observed_nuniv(path, n_nu=200):
    """Universes per group as stored for the first neutrinos (what run_systematics reads)."""
    import uproot
    t = uproot.open(path)['recTree']
    br = ['rec.mc.nu.wgt.univ..length']
    if not all(b in t.keys() for b in br):
        return None
    a = t.arrays(br, entry_stop=50, library='np')['rec.mc.nu.wgt.univ..length']
    rows = [np.asarray(x) for x in a if len(x)]
    return rows[0] if rows else None


def dump_sys(path, knobs, n=6, tree='events/full/selected_nonfid_NgNp_stage1_multisigmaTree'):
    """The stored entries of each knob for the first events it moves. A sigma grid gives the same monotone,
    symmetric shape in every event; random universes give an irregular order (the same one in every event).
    Also counts the NaN / negative entries."""
    import uproot
    f = uproot.open(path)
    t = f[tree] if tree in f else next(f[k] for k in f.keys() if k.split(';')[0].endswith('multisigmaTree'))
    for knob in knobs:
        if knob not in t.keys():
            print(f'{knob}: not in {t.object_path}')
            continue
        W = t[knob].array(entry_stop=20000, library='ak')
        flat = np.asarray(__import__('awkward').flatten(W), float)
        print(f'{knob} ({t.object_path}, first {len(W)} rows): {np.isnan(flat).sum()} NaN and '
              f'{(flat < 0).sum()} negative of {len(flat)} entries')
        if knob + '_sigma' in t.keys():
            print(f'  {knob}_sigma (toml labels): {np.asarray(t[knob + "_sigma"].array(entry_stop=1, library="np")[0])}')
        shown = 0
        for w in W:
            w = np.asarray(w, float)
            if len(w) and np.nanstd(w) > 0.01:
                print(f'  {len(w)} entries: ' + np.array2string(w, precision=3, max_line_width=250))
                shown += 1
                if shown == n:
                    break


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('caf')
    ap.add_argument('--toml', help='systematics toml to check against the CAF')
    ap.add_argument('--sys', help='run_systematics output: also print the stored entries of --knob for a few events')
    ap.add_argument('--knob', nargs='+', default=['GENIEReWeight_SBN_v1_multisigma_MaCCRES',
                                                  'GENIEReWeight_SBN_v1_multisigma_FrCEx_N'])
    a = ap.parse_args(argv)
    if a.sys:
        dump_sys(a.sys, a.knob)
    ps = read_psets(a.caf)
    obs = observed_nuniv(a.caf)
    print(f'{a.caf}: {len(ps)} weight groups')
    print(f'  {"idx":>4} {"type":>4} {"nuniv":>6} {"stored":>6}  name  [parameters]  (sigmas, multisigma)')
    for p in ps:
        st = int(obs[p['i']]) if obs is not None and p['i'] < len(obs) else -1
        prm = p.get('params') or []
        sig = ''
        vv = p.get('vals')
        if vv is not None and len(vv) == 1 and 'multisigma' in p['name']:
            sig = '  sigmas ' + np.array2string(np.asarray(vv[0], float), precision=2, max_line_width=200)
        print(f'  {p["i"]:4d} {p.get("type", -1):4d} {p.get("nuniv", -1):6d} {st:6d}  {p["name"]}  '
              f'[{", ".join(prm[:6])}{", ..." if len(prm) > 6 else ""}]{sig}')
    if not a.toml:
        return 0
    try:
        import tomllib                                   # python >= 3.11: binary file
        with open(a.toml, 'rb') as fh:
            knobs = tomllib.load(fh)['sys']
    except ImportError:
        import toml                                      # older pythons (gpvm 3.9): text
        with open(a.toml) as fh:
            knobs = toml.load(fh)['sys']
    print(f'\n{a.toml}: {len(knobs)} knobs')
    bad = 0
    for k in knobs:
        i = k['index']
        if i >= len(ps):
            print(f'  index {i:3d} {k["name"]}: beyond the CAF\'s {len(ps)} groups')
            bad += 1
        elif ps[i]['name'] != k['name']:
            print(f'  index {i:3d}: toml "{k["name"]}" but the CAF has "{ps[i]["name"]}"')
            bad += 1
    print(f'  {bad} index/name mismatches' if bad else '  every index carries the toml\'s name')
    grid = [(k, ps[k['index']]) for k in knobs if k['type'] == 'multisigma' and k['index'] < len(ps)
            and 'nsigma' in k and ps[k['index']].get('nuniv', len(k['nsigma'])) != len(k['nsigma'])]
    if grid:
        k, p = grid[0]
        print(f'  {len(grid)} multisigma knobs store a different number of entries than their nsigma '
              f'({len(k["nsigma"])}), e.g. {k["name"]}: {p.get("nuniv")} (spineplot then treats them as universes)')
    used = {k['index'] for k in knobs if k['index'] < len(ps)}
    single = {}
    for k in knobs:
        if k['type'] == 'multisigma' and k['index'] < len(ps):
            for prm in ps[k['index']].get('params') or []:
                single.setdefault(prm, []).append(k['name'])
    over = []
    for k in knobs:
        if k['type'] == 'multisim' and k['index'] in used:
            shared = [prm for prm in (ps[k['index']].get('params') or []) if prm in single]
            if shared:
                over.append((k['name'], shared))
    if over:
        print('  parameters varied by a multisim group AND by multisigma knobs in the same toml (counted twice '
              'when every knob is summed):')
        for n, sh in over:
            print(f'    {n}: {", ".join(sh)}')
    elif any(p.get('params') for p in ps):
        print('  no parameter is in both a multisim group and a multisigma knob')
    return 0


if __name__ == '__main__':
    sys.exit(main())
