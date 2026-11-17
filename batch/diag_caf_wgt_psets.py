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
      --toml systematics/toml/gOre_1g1p_sidebands_osc_run2_origin.toml
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


def read_psets(path):
    import uproot
    f = uproot.open(path)
    gt = next((f[k] for k in f.keys() if k.split(';')[0].lower() == 'globaltree'), None)
    if gt is None:
        raise SystemExit(f'{path}: no globalTree; top-level keys: {[k.split(";")[0] for k in f.keys()]}')
    keys = list(gt.keys())
    wk = [k for k in keys if 'wgts' in k]
    need = ['global.wgts.name', 'global.wgts.nuniv']
    if not all(k in keys for k in need):
        print('globalTree weight branches found:\n  ' + '\n  '.join(wk))
        raise SystemExit('unexpected globalTree layout: paste the list above')

    def arr(k):
        return gt[k].array(library='np')[0] if k in keys else None
    names = [str(x) for x in arr('global.wgts.name')]
    psets = [dict(i=i, name=n) for i, n in enumerate(names)]
    for k, field in (('global.wgts.type', 'type'), ('global.wgts.nuniv', 'nuniv')):
        v = arr(k)
        if v is not None:
            for p, x in zip(psets, v):
                p[field] = int(x)
    nmap = arr('global.wgts.map..length')
    pnames = arr('global.wgts.map.param.name')
    vlen = arr('global.wgts.map.vals..length')
    vals = arr('global.wgts.map.vals')
    if nmap is not None and pnames is not None:
        per = split([str(x) for x in pnames], nmap)
        pv = split(list(split(vals, vlen)), nmap) if vals is not None and vlen is not None else [None] * len(per)
        for p, ps, vv in zip(psets, per, pv):
            p['params'] = ps
            p['vals'] = vv
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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('caf')
    ap.add_argument('--toml', help='systematics toml to check against the CAF')
    a = ap.parse_args(argv)
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
        import tomllib
    except ImportError:
        import toml as tomllib  # noqa
    with open(a.toml, 'rb') as fh:
        knobs = tomllib.load(fh)['sys']
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
