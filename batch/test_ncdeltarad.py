"""
Tests for batch/ncdeltarad: GENIE Δ branching-ratio emulation, PDG-table
patching, the back-to-nominal weight, the pilot check, and the orchestrator
state machine against a fake project.py / jobsub_q.
"""
import csv
import os
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

    def __init__(self, cfg, good, queue_ticks=1, token_ok=True, held=False, sidecar_br=0.663):
        self.cfg, self.good, self.queue_ticks = cfg, {k: list(v) for k, v in good.items()}, queue_ticks
        self.token_ok, self.held = token_ok, held
        self.sidecar_br = sidecar_br      # BR(Nγ) the fake gen "used" (None = nominal table)
        self.sl7_scripts = []
        self.calls, self.next_id, self.live = [], 100, {}
        self.dry_run = False

    def log(self, msg):
        pass

    def __call__(self, cmd, mutating=True, timeout=0):
        self.calls.append(cmd)
        prog = cmd[0] if isinstance(cmd, list) else cmd.split()[0]
        if prog == 'apptainer':           # sl7_exec: the root call is the script's last line
            script = Path(cmd.split()[-1].strip("'"))
            self.sl7_scripts.append(script.read_text())
            cmd = script.read_text().strip().splitlines()[-1].split(' ')
            prog = 'root'
        if prog == 'root':                # dump_delta_decays.C("<list>","<csv>")
            out_csv = ' '.join(cmd[4:]).strip("'").rsplit('"', 2)[-2]
            _pilot_csv(out_csv, self.sidecar_br, n=4000, seed=5)
            return 0, ''
        if prog == sys.executable:        # pilot_check.py / delta_br.py run for real
            r = subprocess.run(cmd, capture_output=True, text=True)
            return r.returncode, r.stdout + r.stderr
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


@pytest.fixture(autouse=True)
def no_real_project_py(monkeypatch, request):
    """There is no project.py on a laptop: stub the startup check, except in
    the tests of the check itself."""
    if 'preflight' not in request.node.name:
        monkeypatch.setattr(O, 'preflight', lambda sh: None)


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


def test_reinit_replaces_staged_files(campaign, tmp_path):
    wd = campaign()
    staged = wd.parent / 'pnfs' / 'fcl' / 'ncdeltarad_pdg_table.txt'
    staged.write_text('stale')
    cfg = next(wd.parent.glob('c.toml'))
    O.cmd_init(type('A', (), dict(config=str(cfg), force=True, dry_run=False))())
    assert '6.63000e-01' in staged.read_text()


def test_copy_to_pnfs_removes_destination_first(tmp_path, monkeypatch):
    """dCache refuses to overwrite: the destination must be deleted, not opened for write."""
    src, dst = tmp_path / 'a.txt', tmp_path / 'out' / 'a.txt'
    src.write_text('new')
    dst.parent.mkdir()
    dst.write_text('old')
    real_open, removed = open, []
    real_unlink = Path.unlink

    def guarded_open(path, mode='r', *a, **k):
        if Path(path) == dst and 'w' in mode and dst.exists():
            raise PermissionError('dCache: file exists')
        return real_open(path, mode, *a, **k)

    def tracking_unlink(self, missing_ok=False):
        removed.append(self)
        return real_unlink(self, missing_ok=missing_ok)

    monkeypatch.setattr('builtins.open', guarded_open)
    monkeypatch.setattr(Path, 'unlink', tracking_unlink)
    sh = type('S', (), dict(dry_run=False, log=lambda self, m: None,
                            __call__=lambda self, c, **k: (1, 'ifdh should not be needed')))()
    O.copy_to_pnfs(sh, src, dst)
    assert dst in removed and dst.read_text() == 'new'


def test_job_env_prepends_paths(tmp_path):
    cfg = dict(path_prepend=['/cvmfs/larbatch/bin'], pythonpath_prepend=['/cvmfs/a', '/cvmfs/b'])
    env = O.job_env(cfg, base={'PATH': '/usr/bin', 'PYTHONPATH': '/x'})
    assert env['PATH'] == os.pathsep.join(['/cvmfs/larbatch/bin', '/usr/bin'])
    assert env['PYTHONPATH'] == os.pathsep.join(['/cvmfs/a', '/cvmfs/b', '/x'])
    assert O.job_env({}, base={'PATH': '/usr/bin'}) == {'PATH': '/usr/bin'}
    sh = O.Shell(lambda m: None, env=O.job_env(dict(pythonpath_prepend=[str(tmp_path)])))
    rc, out = sh([sys.executable, '-c', 'import os; print(os.environ["PYTHONPATH"].split(os.pathsep)[0])'])
    assert rc == 0 and out.strip() == str(tmp_path)


def test_dry_run_submits_nothing(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[10]))
    o = O.Orchestrator(wd, shell=O.Shell(lambda m: None, dry_run=True), dry_run=True)
    o.sh = lambda cmd, mutating=True, timeout=0: grid(cmd) if not mutating else (0, '')
    assert o.tick() == 'dry-run'
    o.save()
    assert json.loads((wd / 'state.json').read_text())['stages']['gen']['status'] == 'pending'
    assert not any(isinstance(c, list) and c[0] == 'project.py' for c in grid.calls)


def test_fresh_submit_cleans_stage_first_makeup_does_not(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[8, 10], g4=[10], caf=[10]))
    _drive(wd, grid)
    acts = [(c[4], c[5]) for c in grid.calls if isinstance(c, list) and c[0] == 'project.py']
    gen = [a for s, a in acts if s == 'gen']
    assert gen[:2] == ['--clean_one', '--submit']             # clean right before the fresh submit
    assert '--makeup' in gen and gen[gen.index('--makeup') - 1] != '--clean_one'
    assert gen.count('--clean_one') == 1
    assert ('g4', '--clean_one') in acts and ('g4', '--clean') not in acts


LOCATE = {  # samweb locate-file output, as seen for the run2 off-beam min-bias files
    'a.root': 'dcache:/pnfs/sbn/data_add/sbn_fd/poms_production/x/fc/42\n'
              'enstore:/pnfs/icarus/archive/sbn/sbn_fd/data/compressed_raw/offbeambnbminbias/00/00/95/32(nearline)\n',
    'b.root': 'enstore:/pnfs/icarus/archive/sbn/y/00/00/93/80(nearline)\n'
              'dcache:/pnfs/sbn/data_add/z/3d/5e\n',
    'c.root': 'dcache:/pnfs/sbn/data_add/only/here\n',
}


def _fake_samweb(names):
    def sh(cmd, mutating=True, timeout=0):
        # on AL9 nothing sets SAM_EXPERIMENT: samweb must get -e explicitly
        assert cmd[:3] == ['samweb', '-e', 'icarus'], cmd
        if cmd[3] == 'list-definition-files':
            return 0, ''.join(n + '\n' for n in names)
        if cmd[3] == 'locate-file':
            return 0, LOCATE[cmd[4]]
        raise AssertionError(cmd)
    return sh


def test_gen_must_take_one_overlay_file_per_job(tmp_path):
    base = textwrap.dedent("""\
        [campaign]
        name = "t"
        workdir = "/w"
        pnfs_base = "/p"
        overlay_inputdef = "d"
        gen_numjobs = 2
        {extra}
        [[stage]]
        name = "gen"
        fcl = "gen.fcl"
        maxfilesperjob = {n}
        """)
    p = tmp_path / 'c.toml'
    p.write_text(base.format(extra='', n=5))
    with pytest.raises(SystemExit, match='maxfilesperjob must be 1'):
        O.load_config(p)
    p.write_text(base.format(extra='allow_multifile_gen = true', n=5))
    assert O.load_config(p)['stages'][0]['maxfilesperjob'] == 5
    p.write_text(base.format(extra='', n=1))
    assert O.load_config(p)['stages'][0]['maxfilesperjob'] == 1


def test_parse_location():
    assert O.parse_location('enstore:/pnfs/icarus/archive/a/b(nearline)') == '/pnfs/icarus/archive/a/b'
    assert O.parse_location('dcache:/pnfs/sbn/data_add/x') == '/pnfs/sbn/data_add/x'


def test_resolve_input_list_prefers_archive(tmp_path):
    c = dict(overlay_inputdef='def', input_location_prefer='/pnfs/icarus/archive')
    n = O.resolve_input_list(_fake_samweb(['a.root', 'b.root']), c, tmp_path / 'in.list')
    lines = (tmp_path / 'in.list').read_text().split()
    assert n == 2 and all(p.startswith('/pnfs/icarus/archive/') for p in lines)
    assert lines[0].endswith('/00/00/95/32/a.root') and lines[1].endswith('/00/00/93/80/b.root')


def test_resolve_input_list_fails_without_preferred_replica(tmp_path):
    c = dict(overlay_inputdef='def', input_location_prefer='/pnfs/icarus/archive')
    with pytest.raises(SystemExit, match='1/3 files'):
        O.resolve_input_list(_fake_samweb(['a.root', 'b.root', 'c.root']), c, tmp_path / 'in.list')


def test_locality_and_online_split(tmp_path):
    d = tmp_path / 'archive'
    d.mkdir()
    for name, loc in (('a.root', 'ONLINE_AND_NEARLINE'), ('b.root', 'NEARLINE'), ('c.root', 'ONLINE')):
        (d / name).write_text('x')
        (d / f'.(get)({name})(locality)').write_text(loc + '\n')
    (d / 'd.root').write_text('x')                       # no locality info -> treated as not on disk
    on, off = O.online_split([str(d / n) for n in ('a.root', 'b.root', 'c.root', 'd.root')])
    assert [Path(p).name for p in on] == ['a.root', 'c.root']
    assert [Path(p).name for p in off] == ['b.root', 'd.root']
    assert O.locality(d / 'd.root') == 'UNKNOWN'


def test_gen_stage_streams_with_schema_root(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    cfg['stages'][0]['schema'] = 'root'
    xml = O.render_xml(cfg, {})
    gen = xml.split('<stage name="gen">')[1].split('</stage>')[0]
    g4 = xml.split('<stage name="g4">')[1].split('</stage>')[0]
    assert '<schema>root</schema>' in gen and '<schema>' not in g4


def test_xml_uses_inputlist_when_resolved(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    cfg['overlay_inputlist'] = '/exp/x/overlay_input.list'
    xml = O.render_xml(cfg, {})
    gen = xml.split('<stage name="gen">')[1].split('</stage>')[0]
    assert '<inputlist>/exp/x/overlay_input.list</inputlist>' in gen and '<inputdef>' not in gen


def test_xml_sets_project_os_for_start_stop_jobs(campaign):
    xml = (campaign() / 't.xml').read_text()
    project_part = xml.split('<stage ')[0]
    assert '<os>SL7</os>' in project_part
    assert '--singularity-image' not in xml


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


def test_status_shows_why_a_stage_failed(campaign, capsys):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[0, 0, 0, 0]))   # each tick re-checks before retrying the makeup
    real = grid.__call__

    def failing_makeup(cmd, mutating=True, timeout=0):
        if isinstance(cmd, list) and cmd[0] == 'project.py' and cmd[5] == '--makeup':
            return 1, 'Traceback ...\nRuntimeError: something about the makeup'
        return real(cmd, mutating, timeout)
    _drive(wd, failing_makeup)
    state = json.loads((wd / 'state.json').read_text())
    assert state['stages']['gen']['status'] == 'failed' and state['stages']['gen']['attempts'] == 1
    O.cmd_status(type('A', (), dict(workdir=str(wd)))())
    out = capsys.readouterr().out
    assert '--makeup rc 1' in out and 'something about the makeup' in out


def test_preflight_reports_missing_python_paths():
    def sh(cmd, mutating=True, timeout=0):
        assert cmd == ['project.py', '-h']
        return 1, 'Traceback (most recent call last):\nModuleNotFoundError: No module named \'larbatch_posix\''
    with pytest.raises(SystemExit, match='pythonpath_prepend'):
        O.preflight(sh)
    O.preflight(lambda cmd, mutating=True, timeout=0: (0, 'usage'))     # passes silently


def test_recheck_recovers_a_stage_without_resubmitting(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[10], g4=[10], caf=[10]))
    o = O.Orchestrator(wd, shell=grid, clock=grid_clock)
    o.tick(); o.save()                                    # gen submitted
    st = json.loads((wd / 'state.json').read_text())
    st['stages']['gen'].update(status='failed', last_error='--makeup rc 1: ModuleNotFoundError')
    (wd / 'state.json').write_text(json.dumps(st))       # what the broken run loop left behind
    O.cmd_recheck(type('A', (), dict(workdir=str(wd), stage='gen'))())
    n_before = len(grid.calls)
    results, state = _drive(wd, grid)
    gen_actions = [c[5] for c in grid.calls[n_before:] if isinstance(c, list) and c[0] == 'project.py'
                   and c[4] == 'gen']
    assert gen_actions == ['--check']                     # no clean, no resubmit
    assert state['stages']['gen']['status'] == 'done' and results[-1] == 'done'


def test_recheck_refuses_a_never_submitted_stage(campaign):
    wd = campaign()
    with pytest.raises(SystemExit, match='never submitted'):
        O.cmd_recheck(type('A', (), dict(workdir=str(wd), stage='g4'))())


def test_reconfig_keeps_state_and_input_list(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    cfg['overlay_inputlist'] = str(wd / 'overlay_input.list')
    (wd / 'config.json').write_text(json.dumps(cfg))
    st = json.loads((wd / 'state.json').read_text())
    st['stages']['gen'].update(status='failed', attempts=1)
    (wd / 'state.json').write_text(json.dumps(st))
    toml_path = next(wd.parent.glob('c.toml'))
    toml_path.write_text(toml_path.read_text().replace(
        '[campaign]\n', '[campaign]\npythonpath_prepend = ["/cvmfs/larbatch/python"]\n'))
    O.cmd_reconfig(type('A', (), dict(config=str(toml_path)))())
    new = json.loads((wd / 'config.json').read_text())
    assert new['pythonpath_prepend'] == ['/cvmfs/larbatch/python']
    assert new['overlay_inputlist'] == str(wd / 'overlay_input.list')
    assert json.loads((wd / 'state.json').read_text())['stages']['gen'] == {'status': 'failed', 'attempts': 1}


def test_nonzero_exit_with_job_id_counts_as_submitted(campaign):
    """project.py --makeup printed 'Use job id ...' then exited 1: no duplicate makeup."""
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[1, 10], g4=[10], caf=[10]))
    real = grid.__call__

    def makeup_rc1(cmd, mutating=True, timeout=0):
        rc, out = real(cmd, mutating, timeout)
        if isinstance(cmd, list) and cmd[0] == 'project.py' and cmd[5] == '--makeup':
            return 1, out + '\n\nError:\n\n'
        return rc, out
    makeup_rc1.dry_run, makeup_rc1.log = False, grid.log      # publish uses the shell's attributes
    results, state = _drive(wd, makeup_rc1)
    gen = [c[5] for c in grid.calls if isinstance(c, list) and c[0] == 'project.py' and c[4] == 'gen']
    assert gen.count('--makeup') == 1                       # exactly one makeup, no duplicate
    assert state['stages']['gen']['attempts'] == 2 and state['stages']['gen']['status'] == 'done'
    assert 'rc 1 but submitted' in state['stages']['gen'].get('warning', '')


def test_adopt_tracks_an_external_cluster(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[1, 10], g4=[10], caf=[10]))
    o = O.Orchestrator(wd, shell=grid, clock=grid_clock)
    o.tick(); o.save()                                    # gen submitted (cluster 101)
    st = json.loads((wd / 'state.json').read_text())
    st['stages']['gen'].update(status='checking')         # as left by a "failed" makeup submission
    (wd / 'state.json').write_text(json.dumps(st))
    grid.live['72592754@jobsub03.fnal.gov'] = 1            # the makeup that really is in the queue
    O.cmd_adopt(type('A', (), dict(workdir=str(wd), stage='gen', jobid='72592754.0@jobsub03.fnal.gov'))())
    st = json.loads((wd / 'state.json').read_text())['stages']['gen']
    assert st['status'] == 'running' and st['clusters'] == ['72592754@jobsub03.fnal.gov'] and st['attempts'] == 2
    n = len(grid.calls)
    o = O.Orchestrator(wd, shell=grid, clock=grid_clock)
    assert o.tick() == 'running'                          # polls the adopted cluster
    assert not any(isinstance(c, list) and c[0] == 'project.py' for c in grid.calls[n:])
    with pytest.raises(SystemExit, match='not a jobsub id'):
        O.cmd_adopt(type('A', (), dict(workdir=str(wd), stage='gen', jobid='garbage'))())


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


def test_gen_sidecar_published_with_weights(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    _, state = _drive(wd, FakeGrid(cfg, good=dict(gen=[10], g4=[10], caf=[10])))
    side = Path(state['sidecar']['weights'])
    rows = list(csv.DictReader(open(side)))
    assert len(rows) == 4000 and {'run', 'subrun', 'event', 'weight', 'radiative'} <= set(rows[0])
    rad = [float(r['weight']) for r in rows if r['radiative'] == '1']
    assert 0.005 < np.median(rad) < 0.05            # ~ P_nom / P_enh
    assert 'PASS' in Path(state['sidecar']['pilot_check']).read_text()


def test_sidecar_dump_runs_in_sl7_container_when_configured(campaign):
    wd = campaign(sl7_exec='apptainer exec -B /cvmfs,/exp IMAGE bash {script}')
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[10], g4=[10], caf=[10]))
    _, state = _drive(wd, grid)
    assert state['stages']['gen']['status'] == 'done'
    assert len(grid.sl7_scripts) == 1
    script = grid.sl7_scripts[0]
    assert 'setup icaruscode v10_06_00_13 -q e26:prof' in script and 'dump_delta_decays.C' in script
    assert not any(isinstance(c, list) and c[0] == 'root' for c in grid.calls)   # never run on the host


def test_sl7_script_checks_each_step(tmp_path):
    fakebin = tmp_path / 'bin'
    fakebin.mkdir()
    (fakebin / 'root').write_text('#!/bin/bash\necho "fake root $@"\n')
    (fakebin / 'root').chmod(0o755)
    env = {'PATH': f'{fakebin}:/usr/bin:/bin', 'HOME': str(tmp_path)}

    def run(setup, cmd, path=None):
        script = tmp_path / 's.sh'
        script.write_text(O.sl7_script(setup, cmd))
        e = dict(env, PATH=path or env['PATH'])
        return subprocess.run(['bash', str(script)], capture_output=True, text=True, env=e)

    text = O.sl7_script('source /x/setup_icarus.sh\nsetup icaruscode v1 -q e26:prof', ['root', '-b', 'm.C'])
    assert 'set -e' not in text and 'BEARER_TOKEN_FILE' in text
    # a sourced setup that returns non-zero internally but succeeds overall is fine
    r = run('f() { false; true; }; f\ntrue', ['root', '-b', 'm.C("a","b")'])
    assert r.returncode == 0 and 'fake root -b m.C("a","b")' in r.stdout
    r = run('true\nfalse', ['root', '-b', 'm.C'])
    assert r.returncode == 2 and 'FAILED: false' in r.stderr
    r = run('true', ['root', '-b', 'm.C'], path='/usr/bin:/bin')
    assert r.returncode == 3 and 'no root after setup' in r.stderr


def test_nominal_table_stops_campaign_before_g4(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[10]), sidecar_br=None)
    results, state = _drive(wd, grid)
    assert results[-1] == 'failed' and state['stages']['gen']['status'] == 'failed'
    assert state['stages']['g4']['status'] == 'pending'
    assert not any(isinstance(c, list) and c[0] == 'project.py' and c[4] == 'g4' for c in grid.calls)


def test_sidecar_rejects_duplicate_events(tmp_path):
    p = tmp_path / 'd.csv'
    _pilot_csv(p, 0.663, n=10)
    lines = p.read_text().splitlines()
    p.write_text('\n'.join(lines + [lines[1]]) + '\n')
    with pytest.raises(ValueError, match='duplicated'):
        D.sidecar_weights(p)


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
