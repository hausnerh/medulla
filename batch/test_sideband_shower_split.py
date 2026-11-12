"""
Tests for batch/sideband_shower_split.py: the shower-origin truth groups, the
selection-output input (events/cv with overlay cosmics dropped, dirt from the
same file) and that the sidebands config saves the origin branches.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sideband_shower_split as S  # noqa: E402


def origin(rows):
    return pd.DataFrame(rows, columns=S.ORIGIN)


def test_groups_from_shower_origin():
    d = origin([[1, 0, 0],             # one photon from the nu
                [0, 1, 0],             # one electron from the nu (CC nue)
                [1, 1, 0],             # two from the nu
                [1, 0, 2],             # one from the nu plus pi0 photons
                [0, 0, 2],             # pi0 decay photons only
                [0, 0, 1],             # one reinteraction photon
                [0, 0, 0],             # no true shower
                [np.nan] * 3])         # unmatched
    assert list(S.group_of(d)) == [S.SINGLE_G, S.SINGLE_E, S.MULTI, S.MULTI,
                                   S.SECONDARY, S.SECONDARY, S.OTHER, S.OTHER]


def test_poisson_z_sign_and_zero():
    assert S.poisson_z(10, 10.0) == pytest.approx(0.0)
    assert S.poisson_z(20, 10.0) > 2.5
    assert S.poisson_z(0, 3.0) < -2.0


def write_selection(path, rng):
    import uproot

    def sample(n, cat, ng, ne, ns, nsh):
        return {S.GPS: rng.uniform(0.998, 1, n), S.PPS: rng.uniform(0.96, 1, n),
                S.DM: rng.normal(1232, 20, n), S.NSH: np.full(n, float(nsh)),
                'true_category': np.full(n, float(cat)),
                'reco_leading_primary_gOre_start_dedx': np.full(n, 5.0),
                'reco_leading_primary_gOre_axial_spread': np.zeros(n),
                'reco_leading_primary_gOre_directional_spread': np.full(n, 0.01),
                'reco_leading_primary_gOre_photon_softmax': np.full(n, 0.99),
                'reco_leading_primary_gOre_shower_ke': np.full(n, 200.0),
                S.ORIGIN[0]: np.full(n, float(ng)), S.ORIGIN[1]: np.full(n, float(ne)),
                S.ORIGIN[2]: np.full(n, float(ns))}

    def cat(*parts):
        return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    with uproot.recreate(path) as f:
        # cv: 100 single-gamma (1 shower), 100 secondary (2 showers), 50 overlay cosmics (dropped)
        cv = cat(sample(100, 0, 1, 0, 0, 1), sample(100, 3, 0, 0, 2, 2), sample(50, 9, 0, 0, 0, 1))
        for key, d, pot, lt in (('cv', cv, 1e21, 0.0), ('offbeam', sample(40, 9, 0, 0, 0, 1), 0.0, 200.0),
                                ('onbeam', sample(30, 9, 0, 0, 0, 1), 1e20, 50.0),
                                ('dirt', sample(20, 5, 0, 0, 0, 2), 5e20, 0.0)):
            f[f'events/{key}/{S.TREE}'] = d
            f[f'events/{key}/POT'] = (np.array([pot]), np.array([0., 1.]))
            f[f'events/{key}/Livetime'] = (np.array([lt]), np.array([0., 1.]))


def test_selection_output_input(tmp_path):
    pytest.importorskip('uproot')
    path = tmp_path / 'sel.root'
    write_selection(str(path), np.random.default_rng(0))
    want = [S.GPS, S.PPS, S.DM, S.NSH, 'true_category'] + S.ORIGIN
    a = type('A', (), dict(sys=str(path), dirt=None, exposure_from=None, overlay_cosmics=False))()
    mc, on, pot_on = S.build(a, want)
    assert pot_on == pytest.approx(1e20)
    w = mc.groupby('g')['w'].sum()
    # overlay cosmics dropped; nu rows scaled by POT, off-beam by livetime, dirt from the same file
    assert w[S.SINGLE_G] == pytest.approx(100 * 0.1)
    assert w[S.SECONDARY] == pytest.approx(100 * 0.1)
    assert w[S.COSMIC] == pytest.approx(40 * 50 / 200)
    assert w[S.DIRT] == pytest.approx(20 * 0.2)
    assert S.OVERLAY not in w.index
    assert len(on) == 30
    # --overlay-cosmics: the 50 unmatched cv rows come back as their own POT-scaled group
    a.overlay_cosmics = True
    w2 = S.build(a, want)[0].groupby('g')['w'].sum()
    assert w2[S.OVERLAY] == pytest.approx(50 * 0.1)
    assert w2.drop(S.OVERLAY).to_dict() == pytest.approx(w.to_dict())


def test_one_shower_region_cutflow_and_nm1_outputs(tmp_path):
    pytest.importorskip('uproot')
    path = tmp_path / 'sel.root'
    write_selection(str(path), np.random.default_rng(0))
    out = tmp_path / 'out'
    assert S.main([str(path), '--out', str(out)]) == 0
    flow = pd.read_csv(out / 'cutflow.csv')
    # region = 1 reco shower: the 2-shower secondary rows and 2-shower dirt are outside it
    pre = flow.iloc[0]
    assert pre['stage'] == 'preselection'
    assert pre['True secondary showers'] == pytest.approx(0.0)
    assert pre['True single gamma primary'] == pytest.approx(10.0)
    n_cuts = len(S.cuts(dict(gps=0.9985, pps=0.9513, dm=95, dedx=4.07, axial=-0.5, dirs=0.037, phs=0.9225)))
    assert len(flow) == 1 + 2 * n_cuts          # preselection, cumulative cuts, N-1
    for sub in ('before', 'after', 'nm1'):
        assert any((out / sub).glob('*.png')) or not pytest.importorskip('matplotlib')


def test_missing_origin_branches_explained(tmp_path):
    import uproot
    path = tmp_path / 'old.root'
    with uproot.recreate(str(path)) as f:
        f[f'events/cv/{S.TREE}'] = {S.GPS: np.ones(3)}
    with pytest.raises(SystemExit, match='shower-origin branches'):
        S.load(str(path), 'cv', [S.GPS] + S.ORIGIN)


def test_sidebands_config_saves_origin_branches():
    text = (HERE.parent / 'selection' / 'toml' / 'gOre_1g1p_sidebands_osc_run2.toml').read_text()
    blocks = re.split(r'\n\[\[tree\]\]', text)
    nonfid = [b for b in blocks if re.search(r'name = "selected_nonfid_', b)]
    assert nonfid
    for b in nonfid:
        for v in ('n_nu_photons', 'n_nu_electrons', 'n_secondary_true_showers'):
            assert re.search(rf'name = "{v}",\s*type = "true"', b), (v, b[:80])
