"""
Tests for batch/ncdeltarad: GENIE Δ branching-ratio emulation, PDG-table
patching, the back-to-nominal weight, the pilot check, and the orchestrator
state machine against a fake project.py / jobsub_q.
"""
import csv
import os
import json
import re
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

    def __init__(self, cfg, good, queue_ticks=1, token_ok=True, held=False, sidecar_br=0.663, krb_ok=True):
        self.cfg, self.good, self.queue_ticks = cfg, {k: list(v) for k, v in good.items()}, queue_ticks
        self.token_ok, self.held = token_ok, held
        self.sidecar_br = sidecar_br      # BR(Nγ) the fake gen "used" (None = nominal table)
        self.sl7_scripts = []
        self.krb_ok = krb_ok
        self.calls, self.next_id, self.live, self.removed = [], 100, {}, set()
        self.hold_ad = None               # ClassAd text of the held jobs, for jobsub_q --long
        self.summary = '0 total; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended\n'
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
        if prog in ('klist', 'kinit'):
            return (0 if self.krb_ok else 1), ''
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
                    if st['ana']:         # checkana lists every ROOT file of the job
                        (jd / f'{stage}_{i}.caf.root').write_text('x')
                        good.append(str(jd / f'{stage}_{i}.caf.root'))
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
            cid = cl.split('@')[0]
            if '--long' in cmd:
                if self.hold_ad is None:
                    return 1, 'unknown option --long'
                return 0, '\n\n'.join(f'ClusterId = {cid}\nProcId = {i}\n{self.hold_ad}' for i in (0, 1)) + '\n'
            if '-long' in cmd:
                return 1, 'unknown option -long'
            head = 'JOBSUBJOBID OWNER SUBMITTED RUNTIME ST PRI SIZE CMD\n'
            row = '{}.{}@jobsub01.fnal.gov me 10/20 12:00 0+01:00:00 {} 0 0.0 condor_lar.sh\n'
            if self.held == 'mixed':      # job 0 runs to completion, job 1 stays held until removed
                rows = ''
                if self.live.get(cl, 0) > 0:
                    self.live[cl] -= 1
                    rows += row.format(cid, 0, 'R')
                if f'{cid}.1@jobsub01.fnal.gov' not in self.removed:
                    rows += row.format(cid, 1, 'H')
                return 0, head + rows + self.summary
            if self.live.get(cl, 0) > 0:
                if not self.held:
                    self.live[cl] -= 1
                return 0, head + row.format(cid, 0, 'H' if self.held else 'R') + self.summary
            return 0, head + self.summary
        if prog == 'jobsub_rm':
            jid = cmd[cmd.index('--jobid') + 1]
            self.removed.add(jid)
            self.live[re.sub(r'\.\d+@', '@', jid)] = 0
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
        over.setdefault('empty_polls_to_finish', 1)
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


def test_resolve_input_list_drops_files_without_preferred_replica(tmp_path, capsys):
    c = dict(overlay_inputdef='def', input_location_prefer='/pnfs/icarus/archive')
    n = O.resolve_input_list(_fake_samweb(['a.root', 'b.root', 'c.root']), c, tmp_path / 'in.list')
    assert n == 2 and len((tmp_path / 'in.list').read_text().split()) == 2
    assert '1/3 files have no replica' in capsys.readouterr().out
    with pytest.raises(SystemExit, match='no file of def'):
        O.resolve_input_list(_fake_samweb(['c.root']), c, tmp_path / 'in2.list')


def test_init_uses_first_on_disk_files_of_a_larger_definition(campaign, monkeypatch):
    """10x: tape recalls of the archive copies never came (prestage-dataset
    stages SAM's dead replica); draw more files and keep the on-disk ones."""
    def resolve(sh, c, out_path, threads=8):
        Path(out_path).write_text(''.join(f'/pnfs/icarus/archive/f{i}.root\n' for i in range(10)))
        return 10
    monkeypatch.setattr(O, 'resolve_input_list', resolve)
    monkeypatch.setattr(O, 'online_split', lambda ps: ([p for p in ps if int(p[-6]) % 2 == 0],
                                                       [p for p in ps if int(p[-6]) % 2]))
    wd = campaign(gen_numjobs=3, input_location_prefer='/pnfs/icarus/archive')
    used = (wd / 'overlay_input.list').read_text().split()
    assert used == [f'/pnfs/icarus/archive/f{i}.root' for i in (0, 2, 4)]
    assert len((wd / 'overlay_input_offline.list').read_text().split()) == 5


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


def _outdir(cfg, name):
    return Path(O.stage_dirs(cfg, next(s for s in cfg['stages'] if s['name'] == name))['outdir'])


def test_cleanup_intermediate_removes_inputs_of_done_stages(campaign):
    """Only the final outputs matter: a scratch stage's out/ goes once the
    stage reading it is done; the last stage, logs and lists stay."""
    wd = campaign(cleanup_intermediate=True)
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[8, 10], g4=[10], caf=[10]))
    results, state = _drive(wd, grid)
    assert results[-1] == 'done'
    assert not _outdir(cfg, 'gen').exists() and not _outdir(cfg, 'g4').exists()
    assert any(_outdir(cfg, 'caf').iterdir())
    assert state['stages']['gen']['cleaned'] and state['stages']['g4']['cleaned']
    assert 'cleaned' not in state['stages']['caf']
    assert (wd / 'published' / 'gen' / 'files.list').exists()        # the record stays
    m = json.loads((wd / 'manifest.json').read_text())
    assert m['stages']['gen']['cleaned'] and m['stages']['caf']['cleaned'] is None
    log = (wd / 'orchestrator.log').read_text()
    assert log.index('gen: removed') > log.index('g4: check rc')       # only after g4 consumed gen


def test_cleanup_is_off_by_default_and_command_applies_it_later(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[10], g4=[10], caf=[10]))
    results, _ = _drive(wd, grid)
    assert results[-1] == 'done' and _outdir(cfg, 'gen').exists() and _outdir(cfg, 'g4').exists()
    O.Shell.__call__, real = (lambda self, cmd, mutating=True, timeout=0: (0, '')), O.Shell.__call__
    try:
        O.cmd_cleanup(type('A', (), dict(workdir=str(wd)))())
    finally:
        O.Shell.__call__ = real
    state = json.loads((wd / 'state.json').read_text())
    assert not _outdir(cfg, 'gen').exists() and not _outdir(cfg, 'g4').exists() and _outdir(cfg, 'caf').exists()
    assert json.loads((wd / 'manifest.json').read_text())['stages']['g4']['cleaned'] == state['stages']['g4']['cleaned']


def test_cleanup_never_touches_persistent_keep_or_unfinished(campaign):
    wd = campaign()
    o = O.Orchestrator(wd, shell=FakeGrid(json.loads((wd / 'config.json').read_text()), good={}), clock=grid_clock)
    gen = o.cfg['stages'][0]
    for st, s in ((dict(gen, persistent=True), dict(status='done')),      # larcv / CAF stages
                  (dict(gen, keep=True), dict(status='done')),
                  (dict(gen, publish=[dict(name='larcv', glob='*larcv*.root')]), dict(status='done')),
                  (gen, dict(status='running'))):                          # a makeup may still need it
        out = Path(O.stage_dirs(o.cfg, st)['outdir'])
        out.mkdir(parents=True, exist_ok=True)
        assert o.clean_stage(st, s) == '' and out.exists()


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
    # filesana.list names 2 ROOT files per job: good counts jobs (pilot: "6/2 good")
    assert state['stages']['caf']['good'] == 10


def test_republish_with_exclude_list_finds_renamed_flat_caf(campaign):
    """Pilot: the grid renamed the flat CAF to <input>_<uuid>.root, so the
    *.flat.caf.root glob published nothing. Fix the glob in the TOML,
    reconfig, then `republish` the done stage -- no jobs."""
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[10], g4=[10], caf=[10]))
    results, _ = _drive(wd, grid)
    assert results[-1] == 'done'
    st = next(s for s in cfg['stages'] if s['name'] == 'caf')
    for d in Path(O.stage_dirs(cfg, st)['outdir']).glob('job_[0-9]*'):
        (d / f'in-Overlay_2026_{d.name}.root').write_text('x')      # the renamed flat CAF
        (d / f'Supplemental-in_{d.name}.root').write_text('x')
    st['publish'] = [dict(name='flatcaf', glob='*.root', exclude=['*.caf.root', 'Supplemental-*'])]
    (wd / 'config.json').write_text(json.dumps(cfg))
    O.Shell.__call__, real = (lambda self, cmd, mutating=True, timeout=0: (0, '')), O.Shell.__call__
    try:
        O.cmd_republish(type('A', (), dict(workdir=str(wd), stage='caf'))())
    finally:
        O.Shell.__call__ = real
    caf = json.loads((wd / 'state.json').read_text())['stages']['caf']
    assert caf['good'] == 10 and caf['published']['flatcaf']['n'] == 10
    flat = (wd / 'published' / 'caf' / 'flatcaf_files.list').read_text().split()
    assert all('Overlay_2026' in f for f in flat)
    m = json.loads((wd / 'manifest.json').read_text())
    assert m['stages']['caf']['lists']['flatcaf']['n'] == 10


def _write_flat_caf(path, keys, rad_prim, nslc, ndlp=None, pot=1e18):
    import awkward as ak
    import uproot
    with uproot.recreate(path) as f:
        tree = {
            'rec.hdr.run': np.array([k[0] for k in keys], dtype=np.int32),
            'rec.hdr.subrun': np.array([k[1] for k in keys], dtype=np.int32),
            'rec.hdr.evt': np.array([k[2] for k in keys], dtype=np.int32),
            'rec.nslc': np.array(nslc, dtype=np.int32),
            'rec.mc.nu.iscc': ak.Array([[0] for _ in keys]),
            'rec.mc.nu.genie_mode': ak.Array([[1] for _ in keys]),
            'rec.mc.nu.wgt..length': ak.Array([[189] for _ in keys]),
        }
        if ndlp is not None:              # merged ML CAF: SPINE interactions per event
            tree['rec.dlp..length'] = np.array(ndlp, dtype=np.int32)
            tree['rec.dlp_true..length'] = np.array(ndlp, dtype=np.int32)
        f['recTree'] = tree
        f['TotalPOT'] = (np.array([pot]), np.array([0.0, 1.0]))


def test_check_outputs_ml_flat_cafs(tmp_path):
    pytest.importorskip('uproot')
    import check_outputs as C
    base = tmp_path / 'p'
    keys = [(1, 1, e) for e in range(1, 5)]
    d = base / 'caf' / 'out' / '7_0'
    d.mkdir(parents=True)
    _write_flat_caf(str(d / 'x-Overlay_u.root'), keys, None, [1] * 4)
    (base / 'spine' / 'flatcaf').mkdir(parents=True)
    ml = base / 'spine' / 'flatcaf' / 'x.spine.flat.caf.root'
    _write_flat_caf(str(ml), keys, None, [1] * 4, ndlp=[2, 1, 0, 3])
    (base / 'spine' / 'flatcaf_files.list').write_text(f'{ml}\n')
    by_kind = C.census(str(base / 'caf' / 'out'), __import__('uproot'))
    import awkward as ak
    a = C.load_caf([f for f, _ in by_kind['flat']], __import__('uproot'), ak)
    C.fails.clear(), C.warns.clear()
    C.check_spine(str(base), a, [f for f, _ in by_kind['flat']], __import__('uproot'), ak)
    assert not C.fails, C.fails
    assert C.warns == ['1 events with no SPINE products at all (no HDF5 match in the merge?)']


def test_check_outputs_on_a_synthetic_campaign(tmp_path):
    pytest.importorskip('uproot')
    import check_outputs as C
    base = tmp_path / 'persistent'
    base.mkdir()
    scratch = tmp_path / 'scratch'                # the CAFs are found through the manifest
    keys = [(1, 1, e) for e in range(1, 9)]
    rad = [1, 0, 1, 0, 1, 0, 0, 1]
    for j, part in enumerate((keys[:4], keys[4:])):
        d = scratch / 'caf' / 'out' / f'7_{j}'
        d.mkdir(parents=True)
        _write_flat_caf(str(d / f'x_{j}.flat.caf.root'), part, None, [1] * len(part))
        import uproot
        with uproot.recreate(str(d / f'Supplemental-x_{j}.root')) as f:
            f['h'] = np.histogram([1, 2])
    (base / 'gen_sidecar').mkdir()
    with open(base / 'gen_sidecar' / 'ncdeltarad_weights.csv', 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['run', 'subrun', 'event', 'delta_pdg', 'w_delta', 'channel', 'radiative', 'weight'])
        for k, r in zip(keys, rad):
            w.writerow([*k, 2214, 1.23, '2212:22' if r else '2212:111', r, 0.01 if r else 1.6])
    caflist = scratch / 'caf' / 'flatcaf_files.list'
    caflist.write_text(''.join(f'{p}\n' for p in sorted(scratch.glob('caf/out/*/*.flat.caf.root'))))
    m = dict(name='t', release='r', qual='q', br_gamma=0.663, sidecar={},
             stages=dict(caf=dict(outdir=str(scratch / 'caf' / 'out'), good=2, expected=2,
                                  lists=dict(flatcaf=dict(list=str(caflist), n=2)))))   # as publish() writes it
    (base / 'manifest.json').write_text(json.dumps(m))
    C.fails.clear(), C.warns.clear()
    sys.argv = ['check_outputs.py', str(base)]
    rc = C.main()
    # only the larcv list is missing from this campaign
    assert C.fails == ['stage1 larcv list missing (None)'] and rc == 1, C.fails
    assert not any('flat' in w_ for w_ in C.warns), C.warns


def test_spine_caf_pairs_jobs_with_hdf5_by_event_keys():
    import spine_caf as S
    jobs = {'7_0': {(1, 1, e) for e in range(9)}, '7_1': {(2, 1, e) for e in range(50)}, '7_2': {(3, 1, 0)},
            '7_3': set()}
    h5s = {'b.h5': {(2, 1, e) for e in range(50)}, 'a.h5': {(1, 1, e) for e in range(9)}, 'x.h5': {(9, 9, 9)}}
    pairs, unpaired, empty = S.pair_jobs(jobs, h5s)
    assert pairs == [('7_0', 'a.h5', 9, 9), ('7_1', 'b.h5', 50, 50)] and unpaired == ['7_2'] and empty == ['7_3']


def test_spine_caf_pairing_problems():
    import spine_caf as S
    a, b = {(1, 1, e) for e in range(3)}, {(2, 1, e) for e in range(3)}
    assert S.pairing_problems([('j0', 'a.h5', 3, 3), ('j1', 'b.h5', 3, 3)], [], {'a.h5': a, 'b.h5': b}) == []
    # the same event in two HDF5 files is double processing, whatever the pairing says
    dup = S.pairing_problems([('j0', 'a.h5', 3, 3), ('j1', 'b.h5', 3, 3)], [], {'a.h5': a, 'b.h5': b | {(1, 1, 0)}})
    assert dup == ['event (1, 1, 0) is in both a.h5 and b.h5']
    assert S.pairing_problems([('j0', 'a.h5', 5, 3)], [], {'a.h5': a}) == ['j0: only 3/5 events in a.h5']
    assert S.pairing_problems([], ['j2'], {}) == ['j2: NO SPINE output holds its events']
    assert S.pairing_problems([('j0', 'a.h5', 3, 3)], [], {'a.h5': a, 'b.h5': b}) == ['unused HDF5: b.h5']
    shared = S.pairing_problems([('j0', 'a.h5', 1, 1), ('j1', 'a.h5', 2, 2)], [], {'a.h5': a})
    assert shared == ["a.h5 holds events of several jobs: ['j0', 'j1']"]


def test_spine_caf_h5_destination_never_overwrites_a_source():
    import spine_caf as S
    base = '/pnfs/icarus/scratch/x'
    # already on scratch: no copy, since the merge removes its destination first
    assert S.h5_destination(base, f'{base}/spine/hdf5_incoming/pass1/a.h5', '7_0') is None
    assert S.h5_destination(base, f'{base}/spine/hdf5/a.h5', '7_0') is None
    # larcv names repeat across jobs, so the job goes in the name
    assert S.h5_destination(base, '/exp/d/pass1/a.h5', '7_0') == f'{base}/spine/hdf5/7_0_a.h5'
    assert S.h5_destination(base, '/exp/d/pass2/a.h5', '7_1') == f'{base}/spine/hdf5/7_1_a.h5'
    sc = S.merge_script('setup sbnana', '/b/merge', '/c.caf.root', f'{base}/spine/hdf5_incoming/pass1/a.h5', '/w/t.root',
                        '/w/c.flat.root', '/pnfs/flat/c.flat.root', None)
    assert 'ifdh cp /w/c.flat.root /pnfs/flat/c.flat.root' in sc
    assert 'hdf5_incoming/pass1/a.h5' in sc and 'ifdh rm ' + base not in sc and 'ifdh cp ' + base not in sc


def test_spine_caf_setup_from_build_cache_and_script_order(tmp_path):
    import spine_caf as S
    exe = tmp_path / 'build' / 'merge_sources_simulation'
    exe.parent.mkdir()
    assert S.setup_from_cmakecache(str(exe)) is None
    (exe.parent / 'CMakeCache.txt').write_text(
        'SBNANAOBJ_DIR:PATH=/cvmfs/sbn.opensciencegrid.org/products/sbn/sbnanaobj/v10_01_02/slf7.x86_64.e26.prof\n'
        'CMAKE_PREFIX_PATH=/cvmfs/sbn.opensciencegrid.org/products/sbn/sbnana/v10_01_04/slf7.x86_64.e26.prof/lib\n')
    setup = S.setup_from_cmakecache(str(exe))
    assert setup.endswith('setup sbnana v10_01_04 -q e26:prof') and 'setup_icarus.sh' in setup
    sc = S.merge_script(setup, str(exe), '/pnfs/c.caf.root', '/exp/x.h5', '/w/tmp.root', '/w/c.spine.flat.caf.root',
                        '/pnfs/spine/flatcaf/c.spine.flat.caf.root', '/pnfs/spine/hdf5/x.h5')
    order = [sc.index(s) for s in ('setup sbnana', 'merge_sources_simulation /w/tmp.root /pnfs/c.caf.root /exp/x.h5',
                                   'flatten_caf /w/tmp.root /w/c.spine.flat.caf.root',
                                   'ifdh rm /pnfs/spine/flatcaf/c.spine.flat.caf.root',
                                   'ifdh cp /w/c.spine.flat.caf.root /pnfs/spine/flatcaf/c.spine.flat.caf.root',
                                   'ifdh cp /exp/x.h5 /pnfs/spine/hdf5/x.h5')]
    assert order == sorted(order)


def test_spine_caf_dry_run_pairs_campaign_jobs(campaign, tmp_path, capsys):
    h5py = pytest.importorskip('h5py')
    pytest.importorskip('uproot')
    import spine_caf as S
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    out = Path(O.stage_dirs(cfg, next(s for s in cfg['stages'] if s['name'] == 'caf'))['outdir'])
    h5dir = tmp_path / 'h5'
    h5dir.mkdir()
    good = []
    for j, run in enumerate((9532, 9534)):
        d = out / f'7_{j}'
        d.mkdir(parents=True)
        keys = [(run, 1, e) for e in range(1, 4)]
        _write_flat_caf(str(d / f'in{j}-Overlay_2026_u{j}.root'), keys, None, [1, 1, 1])
        good.append(str(d / f'in{j}-Overlay_2026_v{j}.caf.root'))
        (d / f'in{j}-Overlay_2026_v{j}.caf.root').write_text('x')
        (d / f'Supplemental-in{j}.root').write_text('x')
        ri = np.array(keys, dtype=[('run', 'i4'), ('subrun', 'i4'), ('event', 'i4')])
        with h5py.File(h5dir / f'larcv_mc_{1 - j}_spine.h5', 'w') as f:   # names say nothing about the job
            f['run_info'] = ri
    stale = out / '6_0'                      # a failed attempt, not in files.list, with the same events
    stale.mkdir()
    _write_flat_caf(str(stale / 'in0-Overlay_2026_u9.root'), [(9532, 1, e) for e in range(1, 4)], None, [1, 1, 1])
    (stale / 'in0-Overlay_2026_v9.caf.root').write_text('x')
    book = Path(O.stage_dirs(cfg, next(s for s in cfg['stages'] if s['name'] == 'caf'))['bookdir'])
    book.mkdir(parents=True, exist_ok=True)
    (book / 'filesana.list').write_text(''.join(g + '\n' for g in good))
    cmk = tmp_path / 'ml_cafmaker'           # as build_ml_cafmaker.sh leaves it
    cmk.mkdir()
    (cmk / 'setup_ml_cafmaker.sh').write_text('source /x/setup_icarus.sh\nsetup sbnana v10_06_00_99 -q e26:prof\n'
                                             'export LD_LIBRARY_PATH=/h5/lib:$LD_LIBRARY_PATH\n')
    sys.argv = ['spine_caf.py', str(wd), '--h5', str(h5dir / '*.h5'), '--cafmaker-dir', str(cmk),
                '--dry-run']
    assert S.main() == 0
    log = capsys.readouterr().out
    assert '7_0: 3/3 events in larcv_mc_1_spine.h5' in log and '7_1: 3/3 events in larcv_mc_0_spine.h5' in log
    sc = (wd / 'spine_caf' / 'merge_7_1.sh').read_text()
    assert 'in1-Overlay_2026_v1.caf.root' in sc and 'larcv_mc_0_spine.h5' in sc
    assert 'in1-Overlay_2026_v1.spine.flat.caf.root' in sc
    assert f'{cmk}/sbn_ml_cafmaker/build/merge_sources_simulation' in sc
    assert 'setup sbnana v10_06_00_99 -q e26:prof ||' in sc and 'export LD_LIBRARY_PATH=/h5/lib' in sc
    assert '2 jobs with 6 CAF events, 2 HDF5 files with 6 events, 0 jobs without events' in log
    assert not (wd / 'spine_caf' / 'merge_6_0.sh').exists()
    # the same events in another HDF5: nothing is merged
    for old in (wd / 'spine_caf').glob('merge_*.sh'):
        old.unlink()
    with h5py.File(h5dir / 'again.h5', 'w') as f:
        f['run_info'] = np.array([(9532, 1, 1)], dtype=[('run', 'i4'), ('subrun', 'i4'), ('event', 'i4')])
    assert S.main() == 1
    log = capsys.readouterr().out
    assert 'is in both' in log and 'nothing merged' in log and not list((wd / 'spine_caf').glob('merge_*.sh'))


def test_decode_before_gen_production_order(campaign):
    """Overlay production order: decoder_run2 -> gen -> g4 -> ... The first
    stage reads the overlay data; gen keeps the GENIE table and the sidecar."""
    wd = campaign(stages=('decode', 'gen', 'g4', 'caf'))
    cfg = json.loads((wd / 'config.json').read_text())
    xml = (wd / 't.xml').read_text()
    blocks = {m.group(1): m.group(0) for m in re.finditer(r'<stage name="(\w+)">.*?</stage>', xml, re.S)}
    assert '<inputdef>def</inputdef>' in blocks['decode'] and '<numjobs>10</numjobs>' in blocks['decode']
    assert '<input' not in blocks['gen'] and '<numjobs>1</numjobs>' in blocks['gen']
    assert 'initsource' in blocks['gen'] and 'initsource' not in blocks['decode']
    assert 'ncdeltarad_pdg_table.txt' in blocks['gen'] and 'ncdeltarad_pdg_table.txt' not in blocks['decode']
    grid = FakeGrid(cfg, good=dict(decode=[10], gen=[10], g4=[10], caf=[10]))
    results, state = _drive(wd, grid)
    assert results[-1] == 'done'
    assert state['stages']['gen']['expected'] == 10                 # follows decode's good count
    assert any(isinstance(c, list) and any('pilot_check.py' in str(x) for x in c) for c in grid.calls)  # gen's gate ran


def test_multifile_refused_for_stages_before_gen(tmp_path):
    p = tmp_path / 'c.toml'
    p.write_text('[campaign]\nname = "t"\nworkdir = "/w"\npnfs_base = "/p"\noverlay_inputdef = "d"\ngen_numjobs = 2\n'
                 '[[stage]]\nname = "decode"\nfcl = "d.fcl"\nmaxfilesperjob = 3\n'
                 '[[stage]]\nname = "gen"\nfcl = "g.fcl"\nmaxfilesperjob = 1\n')
    with pytest.raises(SystemExit, match="maxfilesperjob must be 1 for \\['decode', 'gen'\\]"):
        O.load_config(p)
    p.write_text(p.read_text().replace('name = "gen"', 'name = "g4"'))
    with pytest.raises(SystemExit, match="named 'gen' is required"):
        O.load_config(p)


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


def test_one_empty_queue_answer_does_not_finish_the_stage(campaign):
    """10x stage1: jobsub_q came back empty with the jobs still running; the
    check found 2/971 and the makeup processed every input a second time."""
    wd = campaign(empty_polls_to_finish=3)
    o = O.Orchestrator(wd)
    s = o.state['stages']['gen']
    st = next(x for x in o.cfg['stages'] if x['name'] == 'gen')
    s.update(status='running', clusters=['7@jobsub01.fnal.gov'], submitted=0, expected=10, numjobs=10, attempts=1)
    answers = iter([[], [], [('7.0@jobsub01.fnal.gov', 'R')]])
    o.ensure_token = lambda: True
    o.cluster_states = lambda clusters: ({clusters[0]: next(answers)}, '')
    assert o.poll(st, s) == 'waiting' and o.poll(st, s) == 'waiting'
    assert o.poll(st, s) == 'running' and s['empty_polls'] == 0       # the jobs were there after all
    assert s['status'] == 'running' and s['attempts'] == 1


def test_queue_answer_without_summary_line_is_a_failed_query(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[10]))
    grid.summary = ''
    O.Shell.__call__, real = (lambda self, cmd, mutating=True, timeout=0: grid(cmd)), O.Shell.__call__
    try:
        o = O.Orchestrator(wd)
        q, err = o.cluster_states(['7@jobsub01.fnal.gov'])
    finally:
        O.Shell.__call__ = real
    assert q is None and 'no summary line' in err


def test_dedupe_keeps_one_output_per_input():
    u1, u2 = 'ac977456-ac9d-444d-ba2e-08f1246f0e71', '2fe109a2-4d31-492b-8ecd-ae3c1e478d37'
    both_a, both_b = f'/o/30_0/a-Overlay_2026_{u1}.root', f'/o/72_5/a-Overlay_2026_{u2}.root'
    only_old, only_new = f'/o/72_6/b-Overlay_2026_{u1}.root', f'/o/30_1/c-Overlay_2026_{u2}.root'
    files = [both_a, both_b, only_old, only_new]
    kept, dropped = O.dedupe_by_input(files, prefer='72', clean=lambda f: True)
    assert kept == sorted([both_b, only_old, only_new]) and dropped == [both_a]
    # a restarted job's directory is passed over; an input with nothing else is dropped
    kept, dropped = O.dedupe_by_input(files, prefer='72', clean=lambda f: f not in (both_b, only_old))
    assert kept == sorted([both_a, only_new]) and dropped == sorted([both_b, only_old])


def test_single_attempt_spots_a_restarted_job(tmp_path):
    u1, u2 = 'ac977456-ac9d-444d-ba2e-08f1246f0e71', '2fe109a2-4d31-492b-8ecd-ae3c1e478d37'
    one, two = tmp_path / '72_0', tmp_path / '72_1'
    one.mkdir(), two.mkdir()
    for d, uuids in ((one, [u1]), (two, [u1, u2])):
        for i, u in enumerate(uuids):
            (d / f'in-Overlay_2026_{u}.root').write_text('x')
            (d / f'larcv_mc_{i}.root').write_text('x')
            (d / f'Supplemental-in_{u}.root').write_text('x')
    assert O.single_attempt(str(one / f'in-Overlay_2026_{u1}.root'))
    assert not O.single_attempt(str(two / f'in-Overlay_2026_{u1}.root'))


def test_dedupe_command_republishes_the_reduced_list(campaign, capsys):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    O.Shell.__call__, real = (lambda self, cmd, mutating=True, timeout=0: (0, '')), O.Shell.__call__
    try:
        o = O.Orchestrator(wd)
        st = next(x for x in o.cfg['stages'] if x['name'] == 'gen')
        o.state['stages']['gen'].update(status='done', good=4, expected=3)
        o.save()
        book = Path(O.stage_dirs(cfg, st)['bookdir'])
        book.mkdir(parents=True, exist_ok=True)
        u1, u2 = 'ac977456-ac9d-444d-ba2e-08f1246f0e71', '2fe109a2-4d31-492b-8ecd-ae3c1e478d37'
        out = wd.parent / 'o'
        names = [str(out / f'72_{i}' / f'in{i}_{u1}.root') for i in range(3)] + [str(out / '30_0' / f'in0_{u2}.root')]
        for n in names:
            Path(n).parent.mkdir(parents=True)
            Path(n).write_text('x')
            (Path(n).parent / 'larcv_mc.root').write_text('x')
        (book / 'files.list').write_text('\n'.join(names) + '\n')
        args = dict(workdir=str(wd), stage='gen', prefer_cluster='72')
        O.cmd_dedupe(type('A', (), dict(args, dry_run=True))())
        assert len((book / 'files.list').read_text().split()) == 4       # dry run changes nothing
        O.cmd_dedupe(type('A', (), dict(args, dry_run=False))())
    finally:
        O.Shell.__call__ = real
    assert (book / 'files.list').read_text().split() == sorted(names[:3])
    assert len((book / 'files.list.with_duplicates').read_text().split()) == 4
    assert O.Orchestrator(wd).state['stages']['gen']['good'] == 3


def test_more_good_outputs_than_expected_stops_the_chain(campaign):
    """10x stage1: 1912/971 passed the >= 95% test and caf then got 1912 jobs."""
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[11]))
    results, state = _drive(wd, grid)
    assert results[-1] == 'failed' and state['stages']['gen']['status'] == 'failed'
    assert 'duplicate events' in state['stages']['gen']['note']
    assert 'published' not in state['stages']['gen']
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


def test_reconfig_inserts_stage_and_resets_downstream(campaign):
    """Pilot: overlay_waveforms was missing between g4/detsim and stage0. Inserting
    it keeps the finished upstream stages and resets the ones after it."""
    wd = campaign(stages=('gen', 'g4', 'stage0', 'caf'))
    st = json.loads((wd / 'state.json').read_text())
    st['stages']['gen'] = dict(status='done', attempts=1, good=2, expected=2)
    st['stages']['g4'] = dict(status='done', attempts=1, good=2, expected=2)
    st['stages']['stage0'] = dict(status='failed', attempts=3, good=0, expected=2, note='0/2 good')
    (wd / 'state.json').write_text(json.dumps(st))
    toml_path = next(wd.parent.glob('c.toml'))
    toml_path.write_text(toml_path.read_text().replace(
        '[[stage]]\nname = "stage0"', '[[stage]]\nname = "overlay"\nfcl = "overlay_waveforms.fcl"\n\n[[stage]]\nname = "stage0"'))
    O.cmd_reconfig(type('A', (), dict(config=str(toml_path)))())
    new = json.loads((wd / 'state.json').read_text())['stages']
    assert list(new) == ['gen', 'g4', 'overlay', 'stage0', 'caf']
    assert new['g4']['status'] == 'done' and new['g4']['good'] == 2
    assert new['overlay'] == new['stage0'] == new['caf'] == {'status': 'pending'}
    xml = (wd / 't.xml').read_text()
    assert xml.index('<stage name="g4">') < xml.index('<stage name="overlay">') < xml.index('<stage name="stage0">')
    # removing (or reordering) still needs init --force
    toml_path.write_text(re.sub(r'\[\[stage\]\]\nname = "g4"\n[^\[]*', '', toml_path.read_text()))
    with pytest.raises(SystemExit):
        O.cmd_reconfig(type('A', (), dict(config=str(toml_path)))())


def test_reconfig_refuses_insert_before_running_stage(campaign):
    wd = campaign(stages=('gen', 'g4', 'caf'))
    st = json.loads((wd / 'state.json').read_text())
    st['stages']['g4'] = dict(status='running', attempts=1, clusters=['1@jobsub01.fnal.gov'])
    (wd / 'state.json').write_text(json.dumps(st))
    toml_path = next(wd.parent.glob('c.toml'))
    toml_path.write_text(toml_path.read_text().replace(
        '[[stage]]\nname = "g4"', '[[stage]]\nname = "overlay"\nfcl = "o.fcl"\n\n[[stage]]\nname = "g4"'))
    with pytest.raises(SystemExit, match='running'):
        O.cmd_reconfig(type('A', (), dict(config=str(toml_path)))())
    assert 'overlay' not in json.loads((wd / 'state.json').read_text())['stages']


def test_reconfig_while_run_loop_is_alive_updates_xml_memory(campaign):
    """The run loop holds only .runlock between ticks: a reconfig made while it
    sleeps must reach the XML that the next makeup reads (pilot: detsim makeup
    went out at 6000 MB after a reconfig to 20000)."""
    wd = campaign()
    toml_path = next(wd.parent.glob('c.toml'))
    toml_path.write_text(toml_path.read_text().replace('name = "g4"\n', 'name = "g4"\nmemory = 20000\n'))
    o = O.Orchestrator(wd)

    def reconfig_during_run():
        O.cmd_reconfig(type('A', (), dict(config=str(toml_path)))())
        return o.xml_value('g4', 'memory')
    assert O.with_lock(wd, reconfig_during_run, name='.runlock', wait=0) == '20000'
    assert o.xml_value('gen', 'memory') == str(O.STAGE_DEFAULTS['memory'])
    with pytest.raises(SystemExit, match='runlock'):     # still only one run loop
        O.with_lock(wd, lambda: O.with_lock(wd, lambda: None, name='.runlock', wait=0), name='.runlock', wait=0)


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


def test_waits_for_kerberos_ticket(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[10]), krb_ok=False)
    o = O.Orchestrator(wd, shell=grid, clock=grid_clock)
    assert o.tick() == 'auth'
    assert ['kinit', '-R'] in grid.calls                       # tried the non-interactive renewal
    assert not any(isinstance(c, list) and c[0] == 'project.py' for c in grid.calls)
    assert o.state['auth_alerted'] == 'kerberos ticket needed'


def test_auth_failure_in_submit_output_is_not_a_stage_error(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[10]))
    real = grid.__call__

    def vault_fails(cmd, mutating=True, timeout=0):
        if isinstance(cmd, list) and cmd[0] == 'project.py' and cmd[5] == '--submit':
            return 1, ("Failed to process job credential requests (1): ... condor_vault_storer failed\n"
                       "Error: Authentication needed for icarus\nhtgettoken: Failure getting token")
        return real(cmd, mutating, timeout)
    vault_fails.dry_run, vault_fails.log = False, grid.log
    for _ in range(4):
        o = O.Orchestrator(wd, shell=vault_fails, clock=grid_clock)
        assert o.tick() == 'auth'
        o.save()
    gen = json.loads((wd / 'state.json').read_text())['stages']['gen']
    assert gen['status'] == 'pending' and gen.get('errors', 0) == 0      # never counted, never failed


def test_waits_for_token(campaign):
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[10]), token_ok=False)
    results, state = _drive(wd, grid)
    assert results[-1] == 'auth' and state['stages']['gen']['status'] == 'pending'
    assert not any(isinstance(c, list) and c[0] == 'project.py' for c in grid.calls)


def test_failing_queue_query_is_logged_and_notified(campaign):
    """Pilot: jobsub_q timed out on every poll for 2.6 days (stage0) while the
    log showed only the commands -- a failed query must say so, and alert."""
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[10]))
    real = grid.__call__

    def q_hangs(cmd, mutating=True, timeout=0):
        if isinstance(cmd, list) and cmd[0] == 'jobsub_q':
            return 124, 'timeout after 300s'
        return real(cmd, mutating, timeout)
    q_hangs.dry_run, q_hangs.log = False, grid.log
    o = O.Orchestrator(wd, shell=q_hangs, clock=grid_clock)
    assert o.tick() == 'submitted'
    o.save()
    for _ in range(5):
        o = O.Orchestrator(wd, shell=q_hangs, clock=grid_clock)
        assert o.tick() == 'waiting'
        o.save()
    log = (wd / 'orchestrator.log').read_text()
    assert log.count('gen: jobsub_q failed') == 5 and 'timeout after 300s' in log
    assert log.count('NOTIFY') == 1 and 'queue query failing' in log
    gen = json.loads((wd / 'state.json').read_text())['stages']['gen']
    assert gen['status'] == 'running' and gen['query_failures'] == 5
    o = O.Orchestrator(wd, shell=grid, clock=grid_clock)         # query works again
    o.tick()
    o.save()
    assert json.loads((wd / 'state.json').read_text())['stages']['gen'].get('query_failures', 0) == 0


def test_poll_refreshes_token_and_waits_for_auth(campaign):
    """The bearer token outlives no long stage: polling must refresh it too,
    and stop querying (instead of hanging) when that needs a human."""
    wd = campaign()
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[10]), queue_ticks=10)
    o = O.Orchestrator(wd, shell=grid, clock=grid_clock)
    assert o.tick() == 'submitted'
    o.save()
    grid.calls.clear()
    o = O.Orchestrator(wd, shell=grid, clock=grid_clock)
    assert o.tick() == 'running'
    o.save()
    progs = [c[0] for c in grid.calls if isinstance(c, list)]
    assert progs.index('htgettoken') < progs.index('jobsub_q')
    grid.token_ok, grid.krb_ok = False, False
    grid.calls.clear()
    o = O.Orchestrator(wd, shell=grid, clock=grid_clock)
    assert o.tick() == 'auth'
    o.save()
    assert not any(isinstance(c, list) and c[0] == 'jobsub_q' for c in grid.calls)
    assert json.loads((wd / 'state.json').read_text())['stages']['gen']['status'] == 'running'


def test_shell_closes_stdin(tmp_path):
    """An auth prompt from jobsub/htgettoken must fail fast, not wait on the
    tmux terminal until the timeout."""
    sh = O.Shell(lambda m: None)
    rc, out = sh([sys.executable, '-c', 'import sys; print(repr(sys.stdin.read()))'], timeout=10)
    assert rc == 0 and out.strip() == "''"


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


def test_held_job_beside_running_jobs_is_removed_and_made_up(campaign):
    """A job held for memory never runs: it must be removed after the timeout
    even while its siblings run, so the stage reaches its check and makeup
    (pilot: detsim 1 R + 1 H at 6 GB)."""
    wd = campaign(held_timeout_hours=0.5)
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, good=dict(gen=[5, 10], g4=[10], caf=[10]), held='mixed', queue_ticks=20)
    o = O.Orchestrator(wd, shell=grid, clock=grid_clock)
    for _ in range(40):
        o.tick()
        if o.state['stages']['gen']['attempts'] >= 2:
            break
    rm = [c for c in grid.calls if isinstance(c, list) and c[0] == 'jobsub_rm']
    assert rm and all(c[c.index('--jobid') + 1].endswith('.1@jobsub01.fnal.gov') for c in rm)
    assert o.state['stages']['gen']['attempts'] == 2       # removed job redone by the makeup
    assert any(isinstance(c, list) and c[:1] == ['project.py'] and '--makeup' in c for c in grid.calls)


@pytest.mark.parametrize('code,subcode,kinds', [
    (1, 0, ['manual']), (26, 1, ['memory']), (26, 2, ['disk']), (26, 4, ['restarts']), (26, 8, ['lifetime']),
    (26, 3, ['memory', 'disk']), (26, 0, ['other']), (34, 0, ['memory']), (34, 102, ['memory']), (34, 104, ['disk']),
    (46, 0, ['lifetime']), (47, 0, ['lifetime']), (12, 256, ['other']), (3, 0, ['other'])])
def test_classify_hold(code, subcode, kinds):
    assert O.classify_hold(code, subcode) == kinds


def test_parse_hold_ads():
    txt = ('ClusterId = 7\nProcId = 3\nHoldReasonCode = 26\nHoldReasonSubCode = 1\n'
           'HoldReason = "Memory usage exceeded"\nResidentSetSize = 10240000\n'
           'MemoryUsage = ((ResidentSetSize + 1023) / 1024)\nDiskUsage = 3145728\n\n'
           'ClusterId = 7\nProcId = 4\nJobStatus = 2\n')
    ads = O.parse_hold_ads(txt)
    assert list(ads) == ['7.3']
    assert ads['7.3']['kinds'] == ['memory'] and ads['7.3']['reason'] == 'Memory usage exceeded'
    assert ads['7.3']['peak']['memory'] == pytest.approx(10000) and ads['7.3']['peak']['disk'] == pytest.approx(3)


@pytest.mark.parametrize('kind,text,value,raised', [
    ('memory', 6000, 6000, 9000), ('disk', '50GB', 50, '75GB'), ('disk', '500MB', 0.5, '1GB'),
    ('lifetime', '8h', 8, '12h'), ('lifetime', '2d', 48, '72h')])
def test_request_value_and_text(kind, text, value, raised):
    assert O.request_value(kind, text) == pytest.approx(value)
    assert O.request_text(kind, value * 1.5) == raised


def test_named_lifetime_cannot_be_scaled():
    assert O.request_value('lifetime', 'long') is None


def _held_grid(wd, ad, **kw):
    cfg = json.loads((wd / 'config.json').read_text())
    grid = FakeGrid(cfg, **kw)
    grid.hold_ad = ad
    return grid, O.Orchestrator(wd, shell=grid, clock=grid_clock)


MEMORY_AD = 'HoldReasonCode = 26\nHoldReasonSubCode = 1\nHoldReason = "memory usage exceeded"'


def _rm_calls(grid):
    return [c for c in grid.calls if isinstance(c, list) and c[0] == 'jobsub_rm']


def test_memory_hold_raises_the_request_for_the_makeup(campaign):
    wd = campaign(held_timeout_hours=1000.0)
    grid, o = _held_grid(wd, MEMORY_AD, good=dict(gen=[0, 10], g4=[10], caf=[10]), held=True)
    o.tick()                  # submit
    o.tick()                  # held: removed at once, not after the 1000 h
    assert _rm_calls(grid)
    assert o.state['stages']['gen']['raised'] == dict(memory=6000)
    assert '<memory>6000</memory>' in o.xml.read_text()
    grid.held = False
    o.tick()                  # queue empty -> check -> makeup
    s = o.state['stages']['gen']
    assert s['attempts'] == 2 and s['raise_rounds'] == 1
    assert 'code 26/1' in (wd / 'orchestrator.log').read_text()
    assert 'memory usage exceeded' in (wd / 'orchestrator.log').read_text()


def test_raise_rounds_do_not_use_up_makeups(campaign):
    """max_makeup = 0: a stage with no makeups left still gets the round that
    follows a raise, and a later failure with no hold ends it."""
    wd = campaign(held_timeout_hours=1000.0, max_makeup=0, accept_fraction=1.0)
    grid, o = _held_grid(wd, MEMORY_AD, good=dict(gen=[0, 0]), held=True)
    o.tick()
    o.tick()
    grid.held = False
    o.tick()
    assert o.state['stages']['gen']['status'] == 'running' and o.state['stages']['gen']['attempts'] == 2
    for _ in range(4):
        o.tick()
    assert o.state['stages']['gen']['status'] == 'failed'


def test_request_raised_once_per_round(campaign):
    """Jobs of one round held on different polls all ran at the same request."""
    wd = campaign(held_timeout_hours=1000.0)
    grid, o = _held_grid(wd, MEMORY_AD, good=dict(gen=[10], g4=[10], caf=[10]), held=True)
    o.tick()
    s = o.state['stages']['gen']
    st = o.cfg['stages'][0]
    new = O.parse_hold_ads(f'ClusterId = 1\nProcId = 0\n{MEMORY_AD}\n')
    o.act_on_new_holds(st, s, new)
    o.act_on_new_holds(st, s, new)
    assert st['memory'] == 6000 and s['raise_rounds'] == 1


def test_raise_follows_the_peak_when_it_is_larger(campaign):
    wd = campaign()
    grid, o = _held_grid(wd, None, good=dict(gen=[10]))
    st, s = o.cfg['stages'][0], o.state['stages']['gen']
    assert o.raise_request(st, s, 'memory', 10000.0)
    assert st['memory'] == 12000


def test_raise_stops_at_the_cap(campaign):
    wd = campaign(max_memory_mb=5000)
    grid, o = _held_grid(wd, None, good=dict(gen=[10]))
    st, s = o.cfg['stages'][0], o.state['stages']['gen']
    assert o.raise_request(st, s, 'memory', 0.0)
    assert st['memory'] == 5000
    assert not o.raise_request(st, s, 'memory', 0.0)
    assert st['memory'] == 5000
    assert 'at the cap' in (wd / 'orchestrator.log').read_text()


def test_raised_request_survives_a_restart_unless_the_config_asks_for_more(campaign):
    wd = campaign()
    grid, o = _held_grid(wd, None, good=dict(gen=[10]))
    o.state['stages']['gen']['raised'] = dict(memory=9000)
    o.save()
    assert O.Orchestrator(wd, shell=grid, clock=grid_clock).cfg['stages'][0]['memory'] == 9000
    cfg = json.loads((wd / 'config.json').read_text())
    cfg['stages'][0]['memory'] = 12000
    (wd / 'config.json').write_text(json.dumps(cfg))
    assert O.Orchestrator(wd, shell=grid, clock=grid_clock).cfg['stages'][0]['memory'] == 12000


def test_other_hold_waits_for_the_timeout_and_notifies(campaign):
    wd = campaign(held_timeout_hours=1000.0)
    ad = 'HoldReasonCode = 12\nHoldReasonSubCode = 256\nHoldReason = "output transfer failed"'
    grid, o = _held_grid(wd, ad, good=dict(gen=[10]), held=True)
    o.tick()
    o.tick()
    assert not _rm_calls(grid) and 'raised' not in o.state['stages']['gen']
    assert 'jobs held (other) output transfer failed' in (wd / 'orchestrator.log').read_text()


def test_manual_hold_is_left_alone(campaign):
    wd = campaign(held_timeout_hours=0.001)
    grid, o = _held_grid(wd, 'HoldReasonCode = 1\nHoldReason = "via condor_hold"', good=dict(gen=[10]), held=True)
    for _ in range(4):
        o.tick()
    assert not _rm_calls(grid)


def test_hold_reason_is_read_once_per_job(campaign):
    wd = campaign(held_timeout_hours=1000.0)
    ad = 'HoldReasonCode = 12\nHoldReason = "output transfer failed"'
    grid, o = _held_grid(wd, ad, good=dict(gen=[10]), held=True)
    for _ in range(4):
        o.tick()
    assert sum(isinstance(c, list) and '--long' in c for c in grid.calls) == 1


def test_whole_cluster_removed_at_one_go(campaign):
    wd = campaign()
    grid, o = _held_grid(wd, None, good=dict(gen=[10]))
    cl = '7@jobsub01.fnal.gov'
    jobs = [f'7.{i}@jobsub01.fnal.gov' for i in range(3)]
    o.remove_jobs({cl: [(j, 'H') for j in jobs]}, set(jobs))
    assert [c[c.index('--jobid') + 1] for c in _rm_calls(grid)] == [cl]
    grid.calls.clear()
    o.remove_jobs({cl: [(jobs[0], 'H'), (jobs[1], 'H'), (jobs[2], 'R')]}, set(jobs[:2]))
    assert [c[c.index('--jobid') + 1] for c in _rm_calls(grid)] == jobs[:2]


def test_template_keeps_stage_outputs_on_scratch_and_their_products_through_cleanup():
    cfg = O.load_config(str(Path(O.__file__).parent / 'campaign_run2.toml'))
    assert cfg['cleanup_intermediate']
    assert not [s['name'] for s in cfg['stages'] if s['persistent']]
    assert [s['name'] for s in cfg['stages'] if s['publish']] == ['stage1', 'caf']
    stage1 = next(s for s in cfg['stages'] if s['name'] == 'stage1')
    assert O.stage_root(cfg, stage1).startswith(cfg['pnfs_base'])
