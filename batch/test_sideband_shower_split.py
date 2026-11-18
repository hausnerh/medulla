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
    d = origin([[1, 0, 0, 0],          # one photon from the nu
                [0, 1, 0, 0],          # one electron from the nu (CC nue)
                [1, 1, 0, 0],          # two from the nu
                [1, 0, 2, 2],          # one from the nu plus pi0 photons
                [0, 0, 2, 2],          # pi0 decay photons, both start in the argon
                [0, 0, 2, 1],          # pi0 decay photons, one starts outside the argon
                [0, 0, 1, 0],          # one secondary photon, starting outside
                [0, 0, 1, 1],          # one secondary photon, starting inside
                [0, 0, 0, 0],          # no true shower
                [np.nan] * 4])         # unmatched
    assert list(S.group_of(d)) == [S.SINGLE_G, S.SINGLE_E, S.MULTI, S.MULTI, S.SEC_IN, S.SEC_OUT,
                                   S.SEC_OUT, S.SEC_IN, S.OTHER, S.OTHER]


def test_groups_split_pi0_from_other_secondaries():
    d = origin([[0, 0, 2, 2], [0, 0, 2, 2], [0, 0, 1, 1], [0, 0, 1, 1], [1, 0, 0, 0]])
    pi0 = pd.DataFrame([[2, 0, 0],      # pi0, both photons contained in the argon
                        [2, 1, 0],      # pi0, one photon leaves the argon
                        [1, 0, 1],      # pi0 with its partner photon missing from the record
                        [0, 0, 0],      # one secondary shower, not from a pi0 (delta ray / Michel)
                        [0, 0, 0]], columns=S.PI0)
    assert list(S.group_of(pd.concat([d, pi0], axis=1))) == [S.SEC_IN, S.SEC_OUT, S.SEC_OUT, S.SEC_OTHER,
                                                             S.SINGLE_G]


@pytest.mark.parametrize('rng_, nb, cuts', [((0.99, 1.0), 20, [0.9985]), ((0.99, 1.0), 20, [0.9991]),
                                            ((0.0, 1.0), 50, [0.9513]), ((950.0, 1600.0), 26, [1137.0, 1327.0]),
                                            ((950.0, 1600.0), 26, [1172.3, 1291.7]), ((0.0, 0.2), 40, [0.037])])
def test_aligned_edges_put_cuts_on_edges(rng_, nb, cuts):
    e = S.aligned_edges(*rng_, nb, cuts)
    for c in cuts:
        assert np.min(np.abs(e - c)) < 1e-9 * max(1.0, abs(c))
    assert e[0] <= rng_[0] + 1e-12 and e[-1] >= rng_[1] - 1e-12
    w = np.diff(e)
    assert np.allclose(w, w[0]) and 0.5 < w[0] / ((rng_[1] - rng_[0]) / nb) < 2


def test_poisson_z_sign_and_zero():
    assert S.poisson_z(10, 10.0) == pytest.approx(0.0)
    assert S.poisson_z(20, 10.0) > 2.5
    assert S.poisson_z(0, 3.0) < -2.0


def write_selection(path, rng, pi0=False):
    import uproot

    def sample(n, cat, ng, ne, ns, nsh, ns_in=None, npi0=0):
        return {S.GPS: rng.uniform(0.998, 1, n), S.PPS: rng.uniform(0.96, 1, n),
                S.DM: rng.normal(1232, 20, n), S.NSH: np.full(n, float(nsh)),
                'true_category': np.full(n, float(cat)),
                'reco_leading_primary_gOre_start_dedx': np.full(n, 5.0),
                'reco_leading_primary_gOre_axial_spread': np.zeros(n),
                'reco_leading_primary_gOre_directional_spread': np.full(n, 0.01),
                'reco_leading_primary_gOre_photon_softmax': np.full(n, 0.99),
                'reco_leading_primary_gOre_shower_ke': np.full(n, 200.0),
                S.ORIGIN[0]: np.full(n, float(ng)), S.ORIGIN[1]: np.full(n, float(ne)),
                S.ORIGIN[2]: np.full(n, float(ns)), S.ORIGIN[3]: np.full(n, float(ns if ns_in is None else ns_in)),
                **({S.PI0[0]: np.full(n, float(npi0)), S.PI0[1]: np.zeros(n), S.PI0[2]: np.zeros(n),
                    S.TRIG[0]: np.where(np.arange(n) % 10 == 9, 0.0, 1.0)} if pi0 else {})}

    def cat(*parts):
        return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    with uproot.recreate(path) as f:
        # cv: 100 single-gamma (1 shower), 100 secondary (2 showers), 50 overlay cosmics
        cv = cat(sample(100, 0, 1, 0, 0, 1), sample(100, 3, 0, 0, 2, 2, npi0=2), sample(50, 9, 0, 0, 0, 1))
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
    a = type('A', (), dict(sys=str(path), dirt=None, exposure_from=None, no_overlay_cosmics=False,
                           no_trigger_emulation=False))()
    mc, on, pot_on = S.build(a, want)
    assert pot_on == pytest.approx(1e20)
    w = mc.groupby('g')['w'].sum()
    # nu rows scaled by POT, off-beam by livetime, dirt from the same file
    assert w[S.SINGLE_G] == pytest.approx(100 * 0.1)
    assert w[S.SEC_IN] == pytest.approx(100 * 0.1)
    assert w[S.DIRT] == pytest.approx(20 * 0.2)
    # one cosmic category: 40 off-beam rows scaled by livetime + 50 overlay rows (default on) by POT
    assert w[S.COSMIC] == pytest.approx(40 * 50 / 200 + 50 * 0.1)
    assert len(on) == 30
    # --no-overlay-cosmics: the old matched-only scheme
    a.no_overlay_cosmics = True
    w2 = S.build(a, want)[0].groupby('g')['w'].sum()
    assert w2[S.COSMIC] == pytest.approx(40 * 50 / 200)
    assert w2.drop(S.COSMIC).to_dict() == pytest.approx(w.drop(S.COSMIC).to_dict())


def test_untriggered_mc_dropped_and_pi0_groups(tmp_path):
    pytest.importorskip('uproot')
    path = tmp_path / 'sel.root'
    write_selection(str(path), np.random.default_rng(0), pi0=True)
    want = [S.GPS, S.NSH, 'true_category'] + S.ORIGIN
    a = type('A', (), dict(sys=str(path), dirt=None, exposure_from=None, no_overlay_cosmics=False,
                           no_trigger_emulation=False))()
    w = S.build(a, want)[0].groupby('g')['w'].sum()
    # every 10th MC row fails the trigger emulation and is dropped; data rows are untouched
    assert w[S.SINGLE_G] == pytest.approx(90 * 0.1)
    assert w[S.SEC_IN] == pytest.approx(90 * 0.1)            # pi0 photons, contained
    assert w[S.DIRT] == pytest.approx(18 * 0.2)
    assert w[S.COSMIC] == pytest.approx(40 * 50 / 200 + 45 * 0.1)   # off-beam (data: no trigger drop) + overlay
    assert S.GROUPS[S.SEC_IN][0] == S.PI0_LABELS[S.SEC_IN]
    a.no_trigger_emulation = True
    w2 = S.build(a, want)[0].groupby('g')['w'].sum()
    assert w2[S.SINGLE_G] == pytest.approx(100 * 0.1)


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
    assert pre[S.plain(S.GROUPS[S.SEC_IN][0])] == pytest.approx(0.0)
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
        for v in ('n_nu_photons', 'n_nu_electrons', 'n_secondary_true_showers', 'n_secondary_true_showers_in_av',
                  'in_beam_gate', 'n_pi0_true_showers', 'n_pi0_true_showers_escaping', 'n_pi0_lone_photons'):
            assert re.search(rf'name = "{v}",\s*type = "true"', b), (v, b[:80])
        assert re.search(r'name = "trigger_emulated",\s*type = "event"', b), b[:80]
        for v in ('flash', 'flash_time'):
            assert re.search(rf'name = "{v}",\s*type = "reco"', b), (v, b[:80])


def test_marginal_z_reduces_to_poisson_and_shrinks_with_sigma():
    assert S.poisson_tail(7, 2.9, True) == pytest.approx(0.02872, abs=1e-4)
    z0 = S.marginal_z(7, 2.9, 0.0)
    assert z0 == pytest.approx(1.90, abs=0.01)
    assert S.marginal_z(7, 2.9, 1.0) < S.marginal_z(7, 2.9, 0.3) < z0
    assert S.marginal_z(1, 5.0, 0.0) < 0


def test_knob_universe_sums_interpolates_multisigma_like_spineplot():
    S1 = np.ones((1, 4))
    W = np.tile(1 + 0.2 * S.SIGMA_LEVELS[6], (4, 1))          # w(s) = 1 + 0.2 s on every event
    N = S.knob_universe_sums(W, S1, 4000, np.random.default_rng(1))
    s = np.random.default_rng(1).normal(0, 1, 4000)
    assert np.allclose(N[0], 4 * (1 + 0.2 * np.clip(s, -3, 3)))
    Wu = np.array([[1.1, 0.9]] * 4)                            # 2 multisim universes, used as they are
    assert np.allclose(S.knob_universe_sums(Wu, S1, 10, np.random.default_rng(0)), [[4.4, 3.6]])


def write_sys(path, n_sel=40, n_rej=60):
    """A run_systematics-like file: events/full/TREE + row-aligned multisim/multisigma trees."""
    import uproot
    ak = pytest.importorskip('awkward')
    n = n_sel + n_rej
    gps = np.r_[np.full(n_sel, 0.9999), np.full(n_rej, 0.5)]   # only the first n_sel pass every cut
    base = {S.GPS: gps, S.PPS: np.full(n, 0.99), S.DM: np.full(n, 1232.0), S.NSH: np.ones(n),
            'reco_leading_primary_gOre_start_dedx': np.full(n, 5.0),
            'reco_leading_primary_gOre_axial_spread': np.zeros(n),
            'reco_leading_primary_gOre_directional_spread': np.full(n, 0.01),
            'reco_leading_primary_gOre_photon_softmax': np.full(n, 0.99)}
    ms = ak.Array([[1.1, 0.9]] * n)                             # +-10%: sigma = 0.1 N
    sg = ak.Array([list(1 + 0.2 * S.SIGMA_LEVELS[6])] * n)      # 20% per sigma
    with uproot.recreate(path) as f:
        f[f'events/full/{S.TREE}'] = base
        f[f'events/full/{S.TREE}_multisimTree'] = {'flux': ms, 'Run': np.zeros(n, np.int32)}
        f[f'events/full/{S.TREE}_multisigmaTree'] = {'xsec': sg, 'xsec_sigma': ak.Array([[1.0]] * n)}
    return n


def test_sys_band_matches_the_analytic_spread(tmp_path):
    pytest.importorskip('uproot')
    path = str(tmp_path / 'sys.root')
    n = write_sys(path)
    cl = S.cuts(dict(gps=0.9985, pps=0.9513, dm=95, dedx=4.07, axial=-0.5, dirs=0.037, phs=0.9225))
    specs = S.plot_specs(cl)
    after = next(s for s in specs if s['sub'] == 'after' and s['c']['var'] == S.GPS)
    jobs = S.stage_jobs(cl) + [dict(key='after_gps', kind='after', j=None, var=S.GPS, edges=after['edges'])]
    # step=7: the weights are read in chunks that do not divide the 100 rows (chunk edges inside the region)
    b = S.sys_band(path, cl, lambda d: d[S.NSH].to_numpy() == 1, jobs, 0.5, nuniv=20000, step=7)
    cv = [b['cv'][('stage', k)][0] for k in range(len(cl) + 1)]
    var = [b['var'][('stage', k)][0] for k in range(len(cl) + 1)]
    per = b['per_knob']
    assert cv[0] == pytest.approx(0.5 * n) and cv[-1] == pytest.approx(0.5 * 40)
    assert set(per) == {'flux', 'xsec'}
    assert np.sqrt(per['flux'][-1]) == pytest.approx(0.1 * cv[-1])
    assert np.sqrt(per['xsec'][-1]) == pytest.approx(0.2 * cv[-1], rel=0.03)
    assert np.sqrt(var[-1]) == pytest.approx(np.hypot(0.1, 0.2) * cv[-1], rel=0.03)
    # the after-cuts histogram holds the same events in one bin here (all at gps 0.9999): same band
    assert b['cv']['after_gps'].sum() == pytest.approx(cv[-1])
    assert b['var']['after_gps'].sum() == pytest.approx(var[-1])
    assert b['points']['flux'] == pytest.approx([1.1, 0.9])


def test_sys_band_clamps_pathological_weights(tmp_path):
    import uproot
    ak = pytest.importorskip('awkward')
    n = 10
    path = str(tmp_path / 'sys.root')
    base = {S.GPS: np.full(n, 0.9999), S.PPS: np.full(n, 0.99), S.DM: np.full(n, 1232.0), S.NSH: np.ones(n),
            'reco_leading_primary_gOre_start_dedx': np.full(n, 5.0),
            'reco_leading_primary_gOre_axial_spread': np.zeros(n),
            'reco_leading_primary_gOre_directional_spread': np.full(n, 0.01),
            'reco_leading_primary_gOre_photon_softmax': np.full(n, 0.99)}
    w = [[1.0, 1.0]] * n
    w[3] = [1e6, 1.0]                                          # one event with a GENIE-style runaway weight
    with uproot.recreate(path) as f:
        f[f'events/full/{S.TREE}'] = base
        f[f'events/full/{S.TREE}_multisimTree'] = {'k': ak.Array(w)}
    cl = S.cuts(dict(gps=0.9985, pps=0.9513, dm=95, dedx=4.07, axial=-0.5, dirs=0.037, phs=0.9225))
    b = S.sys_band(path, cl, lambda d: d[S.NSH].to_numpy() == 1, S.stage_jobs(cl), 1.0)
    # clamped to 10: universe 0 counts 9 + 10 = 19 vs CV 10, universe 1 counts 10 -> variance 81 / 2
    assert b['var'][('stage', 0)][0] == pytest.approx(40.5)
    # and the clamp is reported: 1 event touched in the region, 1 after all cuts, raw weight 1e6
    assert b['clipped'] == {'k': [1, 1, pytest.approx(1e6), pytest.approx(1.0), 0, 0, 1]}


def test_caf_pset_diagnostic_flags_labels_grids_and_overlaps(tmp_path, monkeypatch, capsys):
    import diag_caf_wgt_psets as D
    ps = [dict(i=0, name='GENIEReWeight_SBN_v1_multisigma_MaCCRES', type=1, nuniv=30, params=['MaCCRES'],
               vals=[np.linspace(-3, 3, 30)]),
          dict(i=1, name='GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse', type=0, nuniv=500,
               params=['MaCCRES', 'MvCCRES'], vals=None),
          dict(i=2, name='some_flux_knob', type=0, nuniv=1000, params=['flux'], vals=None)]
    monkeypatch.setattr(D, 'read_psets', lambda p: ps)
    monkeypatch.setattr(D, 'observed_nuniv', lambda p: np.array([30, 500, 1000]))
    t = tmp_path / 'sys.toml'
    t.write_text('[[sys]]\nname = "GENIEReWeight_SBN_v1_multisigma_MaCCRES"\ntype = "multisigma"\nindex = 0\n'
                 'nsigma = [-1.0, 1.0, -2.0, 2.0, -3.0, 3.0]\n\n'
                 '[[sys]]\nname = "GENIEReWeight_SBN_v1_multisim_CCRESVariationResponse"\ntype = "multisim"\n'
                 'index = 1\n\n[[sys]]\nname = "wrong_label"\ntype = "multisim"\nindex = 2\n')
    assert D.main(['x.root', '--toml', str(t)]) == 0
    out = capsys.readouterr().out
    assert 'index   2: toml "wrong_label" but the CAF has "some_flux_knob"' in out
    assert '1 multisigma knobs store a different number of entries' in out
    assert 'CCRESVariationResponse: MaCCRES' in out
    assert 'sigmas [-3.' in out


def test_uncertainty_table_end_to_end(tmp_path):
    pytest.importorskip('uproot')
    sel = tmp_path / 'sel.root'
    write_selection(str(sel), np.random.default_rng(0))
    sysf = tmp_path / 'sys.root'
    write_sys(str(sysf))
    out = tmp_path / 'out'
    assert S.main([str(sel), '--out', str(out), '--weights', str(sysf)]) == 0
    band = pd.read_csv(out / 'sysband.csv')
    assert len(band) == len(pd.read_csv(out / 'cutflow.csv').query('stage == "preselection" or '
                                                                    'stage.str.startswith("+")'))
    assert (band['syst_frac'].dropna() > 0).all()
    assert (band['z_stat_syst'].abs() <= band['z_stat'].abs() + 1e-9).all()


def test_clamp_report_separates_nan_negative_and_high(tmp_path):
    import uproot
    ak = pytest.importorskip('awkward')
    n = 10
    path = str(tmp_path / 'sys.root')
    base = {S.GPS: np.full(n, 0.9999), S.PPS: np.full(n, 0.99), S.DM: np.full(n, 1232.0), S.NSH: np.ones(n),
            'reco_leading_primary_gOre_start_dedx': np.full(n, 5.0),
            'reco_leading_primary_gOre_axial_spread': np.zeros(n),
            'reco_leading_primary_gOre_directional_spread': np.full(n, 0.01),
            'reco_leading_primary_gOre_photon_softmax': np.full(n, 0.99)}
    w = [[1.0, 1.0]] * n
    w[1] = [np.nan, 1.0]
    w[2] = [-0.5, 1.0]
    w[3] = [50.0, 1.0]
    with uproot.recreate(path) as f:
        f[f'events/full/{S.TREE}'] = base
        f[f'events/full/{S.TREE}_multisimTree'] = {'k': ak.Array(w)}
    cl = S.cuts(dict(gps=0.9985, pps=0.9513, dm=95, dedx=4.07, axial=-0.5, dirs=0.037, phs=0.9225))
    b = S.sys_band(path, cl, lambda d: d[S.NSH].to_numpy() == 1, S.stage_jobs(cl), 1.0)
    assert b['clipped'] == {'k': [3, 3, pytest.approx(50.0), pytest.approx(-0.5), 1, 1, 1]}


def test_caf_pset_diagnostic_reads_split_layout():
    import diag_caf_wgt_psets as D
    ak = pytest.importorskip('awkward')
    assert D.canon('global/wgts/wgts.name') == 'global.wgts.name'
    assert D.canon('global/wgts/wgts.map') == 'global.wgts.map'
    assert D.canon('global.wgts.nuniv') == 'global.wgts.nuniv'
    m = ak.Array([[{'wgts.map.param': {'name': 'MaCCRES', 'id': 3}, 'wgts.map.vals': [-1.0, 1.0, -2.0, 2.0]}],
                  [{'wgts.map.param': {'name': 'MaCCRES', 'id': 3}, 'wgts.map.vals': []},
                   {'wgts.map.param': {'name': 'MvCCRES', 'id': 4}, 'wgts.map.vals': []}]])
    out = D.parse_map(m)
    assert out[0] == (['MaCCRES'], [[-1.0, 1.0, -2.0, 2.0]])
    assert out[1][0] == ['MaCCRES', 'MvCCRES']


def test_sys_band_uses_one_copy_of_repeated_weight_vectors(tmp_path):
    """run_systematics with 5 [[sysvar]]s stored each knob 5 times over (a 6-point multisigma knob as 30
    entries). The band must equal the one from a single copy, not treat the 30 entries as universes."""
    import uproot
    ak = pytest.importorskip('awkward')
    n = 12
    rng = np.random.default_rng(3)
    base = {S.GPS: np.full(n, 0.9999), S.PPS: np.full(n, 0.99), S.DM: np.full(n, 1232.0), S.NSH: np.ones(n),
            'reco_leading_primary_gOre_start_dedx': np.full(n, 5.0),
            'reco_leading_primary_gOre_axial_spread': np.zeros(n),
            'reco_leading_primary_gOre_directional_spread': np.full(n, 0.01),
            'reco_leading_primary_gOre_photon_softmax': np.full(n, 0.99)}
    resp = 1 + 0.2 * rng.random(n)[:, None] * np.array([-1, 1, -2, 2, -3, 3])   # 6 sigma points per event
    univ = 1 + 0.1 * rng.normal(size=(n, 4))                                     # 4 multisim universes
    flat = np.ones((n, 6))                                                       # a knob that never moves
    cl = S.cuts(dict(gps=0.9985, pps=0.9513, dm=95, dedx=4.07, axial=-0.5, dirs=0.037, phs=0.9225))
    bands = []
    for copies in (1, 5):
        path = str(tmp_path / f'sys{copies}.root')
        with uproot.recreate(path) as f:
            f[f'events/full/{S.TREE}'] = base
            f[f'events/full/{S.TREE}_multisimTree'] = {'u': ak.Array(np.tile(univ, copies).tolist())}
            f[f'events/full/{S.TREE}_multisigmaTree'] = {'m': ak.Array(np.tile(resp, copies).tolist()),
                                                         'f': ak.Array(np.tile(flat, copies).tolist())}
        bands.append(S.sys_band(path, cl, lambda d: d[S.NSH].to_numpy() == 1, S.stage_jobs(cl), 1.0))
    assert bands[0]['copies'] == 1 and bands[1]['copies'] == 5
    for k in bands[0]['var']:
        assert bands[1]['var'][k] == pytest.approx(bands[0]['var'][k])
    assert len(bands[1]['points']['m']) == 6
