#!/usr/bin/env python3
"""
Δ(1232) decay-channel probabilities as GENIE R-3_04_02 computes them, the
patched PDG table for an enhanced-radiative NC Δ sample, and the per-event
weight that takes that sample back to nominal.

GENIE (src/Physics/Decay/BaryonResonanceDecayer.cxx, EvolveDeltaBR /
EvolveDeltaDecayWidth) re-evaluates the Δ+ / Δ0 branching ratios for every
decay at the Δ invariant mass W:
    Γ_i(W) = BR_i · Γ · s_i(W),   P_i(W) = Γ_i(W) / Σ_j Γ_j(W)
    s_π(W) = (p_π(W) / p_π(m))³
    s_γ(W) = (E_γ(W) / E_γ(m))³ · (F_γ(W) / F_γ(m))²,  F_γ(E) = 1/(1 + E²/0.706)²
with Γ_i = 0 below the channel's final-state mass. m is the Δ0 table mass,
mN is kProtonMass for a proton daughter and kNucleonMass otherwise, and the
threshold uses the PDG-table daughter masses. Γ cancels in P_i.

Enhancement keeps the two Nπ BRs at their nominal values and raises only the
Nγ BR. Both pion channels then scale by the same factor, so every
non-radiative decay gets the same weight, which needs only (radiative?, W,
Δ charge):
    radiative:      w = P_γ^nom(W) / P_γ^enh(W)
    non-radiative:  w = (1 − P_γ^nom(W)) / (1 − P_γ^enh(W))
W is the Δ invariant mass GENIE decayed with. CAF truth does not have it
(SRTrueInteraction::w is nugen's free-nucleon-at-rest W, off by Fermi motion
and binding), so the inputs come from the gen-level sidecar written by
dump_delta_decays.C, joined on (run, subrun, event).

Usage
  delta_br.py table   --in genie_pdg_table.txt --out ncdeltarad_pdg_table.txt [--br-gamma 0.663]
  delta_br.py split   [--br-gamma 0.663]     # channel split and weights vs W
  delta_br.py weights --in delta_decays.csv --out ncdeltarad_weights.csv [--br-gamma 0.663]
"""
import argparse
import re

import numpy as np

# GENIE R-3_04_02 constants (src/Framework/Conventions/Constants.h) and
# genie_pdg_table.txt masses.
K_PION, K_PI0 = 1.3957018e-01, 1.349766e-01
K_PROTON, K_NEUTRON = 9.38272081e-01, 9.39565413e-01
K_NUCLEON = 0.5 * (K_PROTON + K_NEUTRON)
M_DELTA = 1.233          # utils::res::Mass(kP33_1232) = Δ0 table mass
FF_SCALING = 0.706       # BaryonResonanceDecayer.xml FFScaling
PDG_MASS = {2212: K_PROTON, 2112: K_NEUTRON, 111: K_PI0, 211: K_PION, -211: K_PION, 22: 0.0}

# Nominal channels: (daughters, BR) in table order.
NOMINAL = {
    2214: [((2212, 111), 0.663), ((2112, 211), 0.331), ((2212, 22), 0.006)],
    2114: [((2212, -211), 0.331), ((2112, 111), 0.663), ((2112, 22), 0.006)],
}
# Default enhancement: BR(Nγ) = BR(Nπ0), i.e. radiative : π0 = 1 : 1 at the
# Δ pole, with the charged-pion channel kept at its nominal ratio.
BR_GAMMA_DEFAULT = 0.663


def channel_scaling(daughters, W):
    """s_i(W) for one channel (0 below threshold), as EvolveDeltaDecayWidth."""
    W = np.asarray(W, dtype=float)
    nucleon = next(d for d in daughters if abs(d) in (2212, 2112))
    other = next(d for d in daughters if d != nucleon)
    mN = K_PROTON if nucleon == 2212 else K_NUCLEON
    thresh = sum(PDG_MASS[d] for d in daughters)
    m = M_DELTA
    with np.errstate(invalid='ignore', divide='ignore'):
        if other == 22:
            e_w, e_m = (W**2 - mN**2) / (2 * W), (m**2 - mN**2) / (2 * m)
            f_w, f_m = 1 / (1 + e_w**2 / FF_SCALING)**2, 1 / (1 + e_m**2 / FF_SCALING)**2
            s = (e_w / e_m)**3 * (f_w / f_m)**2
        else:
            mpi = K_PION if abs(other) == 211 else K_PI0
            a1, a2 = (mN + mpi)**2, (mN - mpi)**2
            p_w = np.sqrt((W**2 - a1) * (W**2 - a2)) / (2 * W)
            p_m = np.sqrt((m**2 - a1) * (m**2 - a2)) / (2 * m)
            s = (p_w / p_m)**3
    return np.where(W < thresh, 0.0, np.nan_to_num(s))


def channel_probs(pdg, W, br_gamma=None):
    """P_i(W) for each channel of Δ `pdg` (2214 or 2114). br_gamma replaces the
    Nγ BR (None = nominal). Returns array (nchannel, len(W))."""
    chans = NOMINAL[pdg]
    widths = []
    for daughters, br in chans:
        if 22 in daughters and br_gamma is not None:
            br = br_gamma
        widths.append(br * channel_scaling(daughters, W))
    widths = np.array(widths)
    tot = widths.sum(0)
    return np.divide(widths, tot, out=np.zeros_like(widths), where=tot > 0)


def p_gamma(pdg, W, br_gamma=None):
    return channel_probs(pdg, W, br_gamma)[2]


def weight(radiative, W, delta_pdg, br_gamma=BR_GAMMA_DEFAULT):
    """Per-event weight back to nominal for the enhanced sample.
    radiative: bool array; W: GeV; delta_pdg: 2214 / 2114 per event."""
    radiative = np.asarray(radiative, bool)
    W = np.asarray(W, float)
    delta_pdg = np.asarray(delta_pdg)
    w = np.ones_like(W)
    for pdg in (2214, 2114):
        sel = delta_pdg == pdg
        if not sel.any():
            continue
        pn, pe = p_gamma(pdg, W[sel], None), p_gamma(pdg, W[sel], br_gamma)
        rad = radiative[sel]
        ws = np.where(rad, np.divide(pn, pe, out=np.ones_like(pn), where=pe > 0),
                      np.divide(1 - pn, 1 - pe, out=np.ones_like(pn), where=pe < 1))
        w[sel] = ws
    return w


def sidecar_weights(path, br_gamma=BR_GAMMA_DEFAULT):
    """Per-event weights from the gen-level sidecar (dump_delta_decays.C CSV).
    The exact Δ mass and decay channel come from GENIE's record; CAF truth
    has neither (SRTrueInteraction::w is the free-nucleon-at-rest W).
    Returns a list of dicts keyed by (run, subrun, event); raises on a
    duplicated key, since the join to the CAF would then be ambiguous."""
    import csv
    rows = list(csv.DictReader(open(path)))
    seen, out = set(), []
    for r in rows:
        key = (int(r['run']), int(r['subrun']), int(r['event']))
        if key in seen:
            raise ValueError(f'{path}: duplicated (run, subrun, event) {key}')
        seen.add(key)
        rad = '22' in r['channel'].split(':')
        pdg, wd = int(r['delta_pdg']), float(r['w_delta'])
        out.append(dict(run=key[0], subrun=key[1], event=key[2], delta_pdg=pdg, w_delta=wd,
                        channel=r['channel'], radiative=int(rad),
                        weight=float(weight([rad], [wd], [pdg], br_gamma)[0])))
    return out


def delta_pdg_from_hitnuc(hitnuc):
    """NC conserves charge: struck proton -> Δ+, struck neutron -> Δ0."""
    hitnuc = np.asarray(hitnuc)
    return np.where(hitnuc == 2212, 2214, np.where(hitnuc == 2112, 2114, 0))


_CHAN_RE = re.compile(r'^(\s+\d+\s+\d+\s+)([0-9.eE+-]+)(\s+2\s+)(-?\d+)(\s+)(-?\d+)(.*)$')


def patch_table(text, br_gamma):
    """Replace the Nγ BR of Δ+ (2214) and Δ0 (2114) in genie_pdg_table.txt
    text; everything else is byte-identical. Returns (text, n_replaced)."""
    out, cur, n = [], None, 0
    for line in text.splitlines(keepends=True):
        tok = line.split()
        if len(tok) >= 3 and not line.lstrip().startswith('#') and _is_particle_line(tok):
            cur = int(tok[2])
        m = _CHAN_RE.match(line.rstrip('\n'))
        if m and cur in (2214, 2114) and 22 in (int(m.group(4)), int(m.group(6))):
            old = m.group(2)
            new = f'{br_gamma:.5e}'.rjust(len(old))
            line = line.replace(old, new, 1)
            n += 1
        out.append(line)
    return ''.join(out), n


def _is_particle_line(tok):
    # "<index> <name> <pdg> ..." (decay lines start with two integers then a float)
    try:
        int(tok[0]); int(tok[2])
    except ValueError:
        return False
    return not re.fullmatch(r'-?\d+', tok[1])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    t = sub.add_parser('table', help='write the enhanced PDG table')
    t.add_argument('--in', dest='inp', required=True)
    t.add_argument('--out', required=True)
    t.add_argument('--br-gamma', type=float, default=BR_GAMMA_DEFAULT)
    s = sub.add_parser('split', help='print channel split and weights vs W')
    s.add_argument('--br-gamma', type=float, default=BR_GAMMA_DEFAULT)
    w = sub.add_parser('weights', help='per-event weights CSV from the gen-level sidecar')
    w.add_argument('--in', dest='inp', required=True, help='dump_delta_decays.C output')
    w.add_argument('--out', required=True)
    w.add_argument('--br-gamma', type=float, default=BR_GAMMA_DEFAULT)
    a = ap.parse_args()

    if a.cmd == 'weights':
        import csv
        rows = sidecar_weights(a.inp, a.br_gamma)
        with open(a.out, 'w', newline='') as f:
            wr = csv.DictWriter(f, fieldnames=['run', 'subrun', 'event', 'delta_pdg', 'w_delta', 'channel',
                                               'radiative', 'weight'])
            wr.writeheader()
            wr.writerows(rows)
        print(f'wrote {a.out}: {len(rows)} NC Δ events, BR(Nγ)={a.br_gamma:g}')
        return

    if a.cmd == 'table':
        text, n = patch_table(open(a.inp).read(), a.br_gamma)
        if n != 2:
            raise SystemExit(f'expected to patch 2 Nγ channels (Δ+, Δ0), patched {n}')
        open(a.out, 'w').write(text)
        print(f'wrote {a.out}: BR(Δ+→pγ) = BR(Δ0→nγ) = {a.br_gamma:g} (Nπ BRs unchanged)')
        return

    W = np.array([1.08, 1.10, 1.15, 1.20, 1.233, 1.30, 1.40, 1.60])
    print('W [GeV]          ' + ''.join(f'{x:8.3f}' for x in W))
    for pdg, name, names in ((2214, 'Δ+', ('pπ0', 'nπ+', 'pγ')), (2114, 'Δ0', ('pπ-', 'nπ0', 'nγ'))):
        pn, pe = channel_probs(pdg, W), channel_probs(pdg, W, a.br_gamma)
        for i, cn in enumerate(names):
            print(f'{name} {cn:4s} nominal  ' + ''.join(f'{x:8.4f}' for x in pn[i]))
        for i, cn in enumerate(names):
            print(f'{name} {cn:4s} enhanced ' + ''.join(f'{x:8.4f}' for x in pe[i]))
        print(f'{name} w(radiative)    ' + ''.join(f'{x:8.4f}' for x in weight(np.ones(len(W)), W, np.full(len(W), pdg), a.br_gamma)))
        print(f'{name} w(Nπ)           ' + ''.join(f'{x:8.4f}' for x in weight(np.zeros(len(W)), W, np.full(len(W), pdg), a.br_gamma)))


if __name__ == '__main__':
    main()
