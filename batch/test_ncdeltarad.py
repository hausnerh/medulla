"""
Tests for batch/ncdeltarad: GENIE Δ branching-ratio emulation, PDG-table
patching, the back-to-nominal weight, the pilot check, and the orchestrator
state machine against a fake project.py / jobsub_q.
"""
import csv
import json
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest

NCDR = Path(__file__).resolve().parent / 'ncdeltarad'
sys.path.insert(0, str(NCDR))
import delta_br as D  # noqa: E402
import orchestrate as O  # noqa: E402

# The Δ blocks of GENIE R-3_04_02 genie_pdg_table.txt, verbatim.
TABLE = textwrap.dedent("""\
      186 Delta0             2114  1 100 Baryon       0 1.23300e+00 1.20000e-01 -100 -1 -100 -1   -1   3
    #----------------------------------------------------------------------
    #    decay  type(PY6)    BR     Nd         daughters(codes, then names)
    #----------------------------------------------------------------------
             1     0   3.31000e-01   2     2212     -211                                     p+         pi-
             2     0   6.63000e-01   2     2112      111                                     n0         pi0
             3     0   6.00000e-03   2     2112       22                                     n0         gamma
      187 Delta0_bar      -2114     186     0
      194 Delta+           2214    1 100 Baryon      3 1.23200e+00 1.20000e-01 -100 -1 -100 -1   -1   3
    #----------------------------------------------------------------------
    #    decay  type(PY6)    BR     Nd         daughters(codes, then names)
    #----------------------------------------------------------------------
             1     0   6.63000e-01   2     2212      111                                     p+         pi0
             2     0   3.31000e-01   2     2112      211                                     n0         pi+
             3     0   6.00000e-03   2     2212       22                                     p+         gamma
      195 Delta-          -2214     194     0
      196 Delta++            2224  1 100 Baryon      6 1.23100e+00 1.20000e-01 -100 -1 -100 -1   -1   1
    #----------------------------------------------------------------------
    #    decay  type(PY6)    BR     Nd         daughters(codes, then names)
    #----------------------------------------------------------------------
             1     0   1.00000e+00   2     2212      211                                     p+         pi+
    """)


# --------------------------------------------------------------------------- #
# delta_br
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize('pdg,brs', [(2214, (0.663, 0.331, 0.006)), (2114, (0.331, 0.663, 0.006))])
def test_pole_reproduces_table_brs(pdg, brs):
    np.testing.assert_allclose(D.channel_probs(pdg, np.array([D.M_DELTA]))[:, 0], brs, rtol=1e-9)


def test_enhanced_pole_is_one_to_one():
    p = D.channel_probs(2214, np.array([D.M_DELTA]), 0.663)[:, 0]
    assert p[0] == pytest.approx(p[2])            # pπ0 : pγ = 1 : 1
    assert p[1] / p[0] == pytest.approx(0.331 / 0.663)


def test_below_pion_threshold_is_purely_radiative():
    W = np.array([1.05])                          # below p + π0
    assert D.p_gamma(2214, W)[0] == pytest.approx(1.0)
    assert D.weight([True], W, [2214])[0] == pytest.approx(1.0)


def test_weight_closure():
    """Events drawn with the enhanced BRs, weighted, give the nominal rates."""
    rng = np.random.default_rng(1)
    W = rng.uniform(1.09, 1.6, 400000)
    for pdg in (2214, 2114):
        pe = D.p_gamma(pdg, W, 0.663)
        rad = rng.random(W.size) < pe
        w = D.weight(rad, W, np.full(W.size, pdg), 0.663)
        bins = np.linspace(1.09, 1.6, 6)
        idx = np.digitize(W, bins) - 1
        for b in range(5):
            s = idx == b
            got = (w * rad)[s].sum() / w[s].sum()
            assert got == pytest.approx(D.p_gamma(pdg, W[s]).mean(), rel=0.05)
            assert w[s].sum() / s.sum() == pytest.approx(1.0, rel=0.02)


def test_patch_table_only_touches_gamma_lines():
    out, n = D.patch_table(TABLE, 0.663)
    assert n == 2
    diff = [(a, b) for a, b in zip(TABLE.splitlines(), out.splitlines()) if a != b]
    assert len(diff) == 2
    for a, b in diff:
        assert ' 22 ' in a and '6.63000e-01' in b and len(a) == len(b)


def test_fcl_files_are_ascii():
    """The FHiCL parser rejects non-ASCII bytes, even in comments."""
    for f in NCDR.glob('*.fcl'):
        f.read_bytes().decode('ascii')


def test_hitnuc_to_delta():
    assert list(D.delta_pdg_from_hitnuc([2212, 2112, 0])) == [2214, 2114, 0]


# --------------------------------------------------------------------------- #
# pilot_check
# --------------------------------------------------------------------------- #
def _pilot_csv(path, br_gamma, n=20000, seed=2):
    rng = np.random.default_rng(seed)
    with open(path, 'w', newline='') as f:
        wr = csv.writer(f)
        wr.writerow(['run', 'subrun', 'event', 'delta_pdg', 'w_delta', 'w_kine', 'hitnuc', 'channel',
                     'n_daughters', 'n_photon_status1_gt30mev'])
        for i in range(n):
            pdg = 2214 if rng.random() < 0.5 else 2114
            W = min(max(rng.normal(1.232, 0.08), 1.08), 1.8)
            p = D.channel_probs(pdg, np.array([W]), br_gamma)[:, 0]
            ch = rng.choice(3, p=p)
            chans = D.NOMINAL[pdg][ch][0]
            wr.writerow([1, 1, i, pdg, W, W, 2212 if pdg == 2214 else 2112, ':'.join(map(str, chans)), 2,
                         int(22 in chans)])


def test_pilot_check_passes_on_enhanced_sample(tmp_path):
    _pilot_csv(tmp_path / 'ok.csv', 0.663)
    r = subprocess.run([sys.executable, str(NCDR / 'pilot_check.py'), str(tmp_path / 'ok.csv')],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr


def test_pilot_check_catches_nominal_table(tmp_path):
    _pilot_csv(tmp_path / 'nominal.csv', None)
    r = subprocess.run([sys.executable, str(NCDR / 'pilot_check.py'), str(tmp_path / 'nominal.csv')],
                       capture_output=True, text=True)
    assert r.returncode != 0 and 'table not applied' in r.stdout


# --------------------------------------------------------------------------- #
# orchestrator
# --------------------------------------------------------------------------- #
class FakeGrid:
    """Stands in for project.py / jobsub_q / htgettoken. `good` maps stage ->
    list of good-output counts returned by successive checks."""

    def __init__(self, cfg, good, queue_ticks=1, token_ok=True, held=False):
        self.cfg, self.good, self.queue_ticks = cfg, {k: list(v) for k, v in good.items()}, queue_ticks
        self.token_ok, self.held = token_ok, held
        self.calls, self.next_id, self.live = [], 100, {}
        self.dry_run = False

    def log(self, msg):
        pass

    def __call__(self, cmd, mutating=True, timeout=0):
        self.calls.append(cmd)
        prog = cmd[0] if isinstance(cmd, list) else cmd.split()[0]
        if prog == 'htgettoken':
            return (0 if self.token_ok else 1), ''
        if prog == 'project.py':
            stage, action = cmd[4], cmd[5]
            if action in ('--submit', '--makeup'):
                self.next_id += 1
                cl = f'{self.next_id}@jobsub01.fnal.gov'
                self.live[cl] = self.queue_ticks
                return 0, f'Use job id {self.next_id}.0@jobsub01.fnal.gov to retrieve output'
            if action in ('--check', '--checkana'):
                st = next(s for s in self.cfg['stages'] if s['name'] == stage)
                dirs = O.stage_dirs(self.cfg, st)
                book = Path(dirs['bookdir'])
                book.mkdir(parents=True, exist_ok=True)
                n = self.good[stage].pop(0)
                good = []
                for i in range(n):          # one output directory per good job
                    jd = Path(dirs['outdir']) / f'job_{i}'
                    jd.mkdir(parents=True, exist_ok=True)
                    main = jd / (f'{stage}_{i}.flat.caf.root' if st['ana'] else f'{stage}_{i}.root')
                    main.write_text('x')
                    good.append(str(main))
                    if st['ana']:
                        (jd / f'{stage}_{i}.caf.root').write_text('x')
                    if stage == 'stage1':
                        (jd / f'larcv_{i}.root').write_text('x')
                # leftovers of a failed job must never be published
                bad = Path(dirs['outdir']) / 'job_failed'
                bad.mkdir(parents=True, exist_ok=True)
                (bad / 'larcv_failed.root').write_text('x')
                name = 'filesana.list' if st['ana'] else 'files.list'
                (book / name).write_text(''.join(f + '\n' for f in good))
                return 0, ''
        if prog == 'jobsub_q':
            cl = cmd[cmd.index('--jobid') + 1]
            if self.live.get(cl, 0) > 0:
                if not self.held:
                    self.live[cl] -= 1
                st = 'H' if self.held else 'R'
                return 0, (f'JOBSUBJOBID OWNER SUBMITTED RUNTIME ST PRI SIZE CMD\n'
                           f'{cl.split("@")[0]}.0@jobsub01.fnal.gov me 10/20 12:00 0+01:00:00 {st} 0 0.0 condor_lar.sh\n')
            return 0, 'JOBSUBJOBID OWNER SUBMITTED RUNTIME ST PRI SIZE CMD\n'
        if prog == 'jobsub_rm':
            self.live[cmd[cmd.index('--jobid') + 1]] = 0
            return 0, ''
        return 0, ''


@pytest.fixture
def campaign(tmp_path):
    def make(stages=('gen', 'g4', 'caf'), gen_numjobs=10, **over):
        st = '\n'.join(f'[[stage]]\nname = "{s}"\nfcl = "{"simulation_genie_icarus_bnb_overlays_volDetEnclosure_fiducial_ncdeltarad.fcl" if s == "gen" else s + ".fcl"}"\n'
                       + ('ana = true\n' if s == 'caf' else '') for s in stages)
        extra = '\n'.join(f'{k} = {json.dumps(v)}' for k, v in over.items())
        (tmp_path / 'table.txt').write_text(TABLE)
        cfg = tmp_path / 'c.toml'
        cfg.write_text(f'[campaign]\nname = "t"\nworkdir = "{tmp_path / "wd"}"\npnfs_base = "{tmp_path / "pnfs"}"\n'
                       f'overlay_inputdef = "def"\ngen_numjobs = {gen_numjobs}\n'
                       f'genie_pdg_table = "{tmp_path / "table.txt"}"\n{extra}\n{st}')
        O.cmd_init(type('A', (), dict(config=str(cfg), force=False, dry_run=False))())
        return tmp_path / 'wd'
    return make


def _drive(wd, grid, max_ticks=50):
    results = []
    for _ in range(max_ticks):
        o = O.Orchestrator(wd, shell=grid, clock=grid_clock)
        r = o.tick()
        o.save()
        results.append(r)
        if r in ('done', 'failed', 'auth'):
            break
    return results, json.loads((wd / 'state.json').read_text())


_t = [0.0]


def grid_clock():
    _t[0] += 3600.0          # each tick = one hour
    return _t[0]


def test_init_stages_files_and_renders_xml(campaign):
    wd = campaign()
    pnfs_fcl = wd.parent / 'pnfs' / 'fcl'
    assert '6.63000e-01' in (pnfs_fcl / 'ncdeltarad_pdg_table.txt').read_text()
    env = (pnfs_fcl / 'ncdeltarad_genie_env.sh').read_text()
    assert '@NCDR' not in env and 'NCDR_BR_GAMMA="0.663"' in env
    xml = (wd / 't.xml').read_text()
    assert '<inputdef>def</inputdef>' in xml and '<initsource>' in xml
    assert xml.count('<initsource>') == 1 and '-f ' in xml.split('<stage name="g4">')[0]


def test_full_chain_with_makeup(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[8, 10], g4=[10], caf=[10]))
    results, state = _drive(wd, grid)
    assert results[-1] == 'done'
    assert state['stages']['gen']['attempts'] == 2          # one makeup round
    assert state['stages']['g4']['expected'] == 10          # follows gen's good count
    assert all(s['status'] == 'done' for s in state['stages'].values())
    assert any(c[5] == '--checkana' for c in grid.calls if isinstance(c, list) and c[0] == 'project.py')


def test_accepts_below_target_then_propagates(campaign):
    wd = campaign(max_makeup=1)
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[8, 9], g4=[9], caf=[9]))   # 90% < 95%, >= 80%
    _, state = _drive(wd, grid)
    assert state['stages']['gen']['status'] == 'done'
    assert state['stages']['g4']['expected'] == 9


def test_fails_below_accept_fraction(campaign):
    wd = campaign(max_makeup=1)
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[5, 6]))
    results, state = _drive(wd, grid)
    assert results[-1] == 'failed' and state['stages']['gen']['status'] == 'failed'
    assert state['stages']['g4']['status'] == 'pending'


def test_waits_for_token(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[10]), token_ok=False)
    results, state = _drive(wd, grid)
    assert results[-1] == 'auth' and state['stages']['gen']['status'] == 'pending'
    assert not any(isinstance(c, list) and c[0] == 'project.py' for c in grid.calls)


def test_publishes_lists_and_manifest(tmp_path):
    (tmp_path / 'table.txt').write_text(TABLE)
    cfgp = tmp_path / 'c.toml'
    cfgp.write_text(textwrap.dedent(f"""\
        [campaign]
        name = "t"
        workdir = "{tmp_path / 'wd'}"
        pnfs_base = "{tmp_path / 'scratch'}"
        persistent_base = "{tmp_path / 'persistent'}"
        overlay_inputdef = "def"
        gen_numjobs = 4
        genie_pdg_table = "{tmp_path / 'table.txt'}"
        [[stage]]
        name = "gen"
        fcl = "simulation_genie_icarus_bnb_overlays_volDetEnclosure_fiducial_ncdeltarad.fcl"
        [[stage]]
        name = "stage1"
        fcl = "stage1.fcl"
        persistent = true
        publish = [ {{ name = "larcv", glob = "*larcv*.root" }} ]
        [[stage]]
        name = "caf"
        fcl = "caf.fcl"
        ana = true
        persistent = true
        publish = [ {{ name = "flatcaf", glob = "*.flat.caf.root" }},
                    {{ name = "caf", glob = "*.caf.root", exclude = "*.flat.caf.root" }} ]
        """))
    O.cmd_init(type('A', (), dict(config=str(cfgp), force=False, dry_run=False))())
    wd = tmp_path / 'wd'
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[4], stage1=[4], caf=[4]))
    results, state = _drive(wd, grid)
    assert results[-1] == 'done'

    s1 = tmp_path / 'persistent' / 'stage1'
    larcv = (s1 / 'larcv_files.list').read_text().split()
    assert len(larcv) == 4 and not any('failed' in f for f in larcv)
    assert len((s1 / 'files.list').read_text().split()) == 4
    caf = tmp_path / 'persistent' / 'caf'
    flat = (caf / 'flatcaf_files.list').read_text().split()
    plain = (caf / 'caf_files.list').read_text().split()
    assert len(flat) == 4 and len(plain) == 4 and not any(f.endswith('.flat.caf.root') for f in plain)
    assert (tmp_path / 'scratch' / 'gen' / 'files.list').exists()     # non-persistent stage stays on scratch

    m = json.loads((tmp_path / 'persistent' / 'manifest.json').read_text())
    assert m['stages']['stage1']['lists']['larcv']['n'] == 4
    assert m['stages']['caf']['lists']['flatcaf']['list'] == str(caf / 'flatcaf_files.list')
    assert m['br_gamma'] == pytest.approx(0.663) and m['event_generator_list'] == 'NCRES'
    assert '6.63000e-01' in (tmp_path / 'persistent' / 'ncdeltarad_pdg_table.txt').read_text()
    assert (tmp_path / 'persistent' / 'campaign_config.json').exists()


def test_removes_long_held_jobs(campaign):
    wd = campaign(held_timeout_hours=2.0)
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[0, 10], g4=[10], caf=[10]), held=True)
    o = O.Orchestrator(wd, shell=grid, clock=grid_clock)
    for _ in range(6):
        o.tick()
        if any(isinstance(c, list) and c[0] == 'jobsub_rm' for c in grid.calls):
            break
    assert any(isinstance(c, list) and c[0] == 'jobsub_rm' for c in grid.calls)
