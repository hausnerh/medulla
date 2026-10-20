#!/usr/bin/env python3
"""
Orchestrator for the enhanced-radiative NC Δ grid campaign.

Drives the ICARUS production chain stage by stage with larbatch's project.py
(gen -> g4 -> detsim -> stage0 -> stage1 -> caf). The SPINE step needs a GPU
node and is not run here: the campaign's product is the stage outputs on
/pnfs plus file lists and a manifest the SPINE step reads. Optional
"external" stages (shell commands) can still be chained after caf. For each
larbatch stage it:

  1. renders the project.py XML (numjobs follows the previous stage's good
     output, so no input file is silently dropped),
  2. refreshes the token and submits (`project.py --submit`),
  3. watches the submitted clusters with `jobsub_q`, removing jobs held
     longer than `held_timeout_hours`,
  4. runs `project.py --check` (or `--checkana`) once the clusters are gone,
  5. resubmits failures with `project.py --makeup` up to `max_makeup` times,
     then accepts the stage (>= accept_fraction good) or marks it failed,
  6. publishes the stage: `files.list` (good outputs per project.py) and one
     `<name>_files.list` per [[stage.publish]] glob (e.g. larcv, flat CAF),
     written next to the stage's outdir on /pnfs.

When every larbatch stage is done it writes <persistent_base>/manifest.json
(stage outdirs, counts, list paths) and copies the patched PDG table and the
campaign config next to it.

State lives in <workdir>/state.json (written atomically, guarded by a lock),
so `tick` can run from cron or `run` can loop in tmux, and either can be
killed and restarted at any point.

    orchestrate.py init   campaign.toml          # render, stage files, write state
    orchestrate.py status <workdir>
    orchestrate.py tick   <workdir> [--dry-run]  # one pass of the state machine
    orchestrate.py run    <workdir> [--interval 900]
    orchestrate.py retry  <workdir> <stage>      # reset a failed stage
"""
import argparse
import fcntl
import fnmatch
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

try:
    import tomllib

    def _read_toml(path):
        with open(path, 'rb') as f:
            return tomllib.load(f)
except ImportError:  # python < 3.11 (the gpvm default): the toml package
    import toml

    def _read_toml(path):
        return toml.load(path)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import delta_br  # noqa: E402

DEFAULTS = dict(experiment='icarus', release='v10_06_00_13', qual='e26:prof', file_type='overlay',
                run_type='physics', resource='DEDICATED,OPPORTUNISTIC', jobsub_extra='',
                br_gamma=delta_br.BR_GAMMA_DEFAULT, genie_pdg_table='', min_good_fraction=0.95,
                max_makeup=2, accept_fraction=0.80, held_timeout_hours=6.0, fallback_wait_hours=12.0,
                notify_cmd='', persistent_base='', numevents=1000000000, delta_sidecar=True,
                xrootd_prefix='root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/')
STAGE_DEFAULTS = dict(datatier='', memory=4000, disk='20GB', lifetime='8h', ana=False,
                      maxfilesperjob=1, persistent=False, jobsub='', publish=[])
JOBID_RE = re.compile(r'(\d+)\.\d+@(jobsub\d+\.fnal\.gov)')
STATE_TOKENS = {'I', 'R', 'H', 'X', 'C', 'S'}


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


# --------------------------------------------------------------------------- #
# configuration and rendering
# --------------------------------------------------------------------------- #
def load_config(path):
    cfg = _read_toml(path)
    c = {**DEFAULTS, **cfg['campaign']}
    c['stages'] = [{**STAGE_DEFAULTS, **s} for s in cfg.get('stage', [])]
    c['external'] = cfg.get('external', [])
    for k in ('name', 'workdir', 'pnfs_base', 'overlay_inputdef', 'gen_numjobs'):
        if not c.get(k):
            raise SystemExit(f'config: campaign.{k} is required')
    if not c['stages'] or c['stages'][0]['name'] != 'gen':
        raise SystemExit("config: the first [[stage]] must be 'gen'")
    return c


def stage_dirs(c, st):
    root = stage_root(c, st)
    return dict(outdir=f'{root}/out', logdir=f'{root}/log', workdir=f"{c['pnfs_base']}/{st['name']}/work",
                bookdir=f"{c['workdir']}/book/{st['name']}")


def stage_root(c, st):
    """/pnfs directory holding the stage's out/ and log/ and its published lists."""
    base = c['persistent_base'] if (st['persistent'] and c['persistent_base']) else c['pnfs_base']
    return f"{base}/{st['name']}"


def manifest_dir(c):
    return c['persistent_base'] or c['pnfs_base']


def render_xml(c, numjobs):
    """project.py XML; numjobs: stage name -> number of jobs."""
    fcldir = f"{c['pnfs_base']}/fcl"
    lines = ['<?xml version="1.0"?>', '<job>', f'<project name="{c["name"]}">',
             f'  <numevents>{c["numevents"]}</numevents>',
             f'  <resource>{c["resource"]}</resource>',
             '  <larsoft>', f'    <tag>{c["release"]}</tag>', f'    <qual>{c["qual"]}</qual>', '  </larsoft>',
             f'  <version>{c["release"]}</version>', f'  <filetype>{c["file_type"]}</filetype>',
             f'  <runtype>{c["run_type"]}</runtype>', f'  <fcldir>{fcldir}</fcldir>']
    for st in c['stages']:
        d = stage_dirs(c, st)
        jobsub = f"--expected-lifetime={st['lifetime']} {c['jobsub_extra']} {st['jobsub']}"
        if st['name'] == 'gen':
            jobsub += f" -f {fcldir}/ncdeltarad_pdg_table.txt"
        lines += [f'  <stage name="{st["name"]}">', f'    <fcl>{st["fcl"]}</fcl>']
        if st['name'] == 'gen':
            lines += [f'    <inputdef>{c["overlay_inputdef"]}</inputdef>',
                      f'    <initsource>{fcldir}/ncdeltarad_genie_env.sh</initsource>']
        lines += [f'    <{k}>{v}</{k}>' for k, v in d.items()]
        default = c['gen_numjobs'] if st['name'] == 'gen' else 1
        lines += [f'    <numjobs>{numjobs.get(st["name"], default)}</numjobs>',
                  f'    <maxfilesperjob>{st["maxfilesperjob"]}</maxfilesperjob>',
                  f'    <memory>{st["memory"]}</memory>', f'    <disk>{st["disk"]}</disk>',
                  f'    <jobsub>{" ".join(jobsub.split())}</jobsub>']
        if st['datatier']:
            lines.append(f'    <datatier>{st["datatier"]}</datatier>')
        if st['ana']:
            lines.append('    <ana>1</ana>')
        lines.append('  </stage>')
    lines += ['</project>', '</job>', '']
    return '\n'.join(lines)


# --------------------------------------------------------------------------- #
# shell backend (one replaceable entry point, so the logic is testable)
# --------------------------------------------------------------------------- #
class Shell:
    def __init__(self, log, dry_run=False):
        self.log, self.dry_run = log, dry_run

    def __call__(self, cmd, mutating=True, timeout=1800):
        """Run `cmd` (list or shell string). In dry-run, mutating commands are
        only logged. Returns (returncode, stdout+stderr)."""
        text = cmd if isinstance(cmd, str) else ' '.join(shlex.quote(x) for x in cmd)
        if self.dry_run and mutating:
            self.log(f'[dry-run] {text}')
            return 0, ''
        self.log(f'$ {text}')
        try:
            p = subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            return 124, f'timeout after {timeout}s'
        return p.returncode, (p.stdout or '') + (p.stderr or '')


class Orchestrator:
    def __init__(self, workdir, shell=None, dry_run=False, clock=time.time):
        self.dir = Path(workdir)
        self.cfg = json.loads((self.dir / 'config.json').read_text())
        self.state_path = self.dir / 'state.json'
        self.state = json.loads(self.state_path.read_text())
        self.xml = self.dir / f"{self.cfg['name']}.xml"
        self.sh = shell or Shell(self.log, dry_run)
        self.dry_run = dry_run
        self.clock = clock

    # ---- persistence / logging -------------------------------------------- #
    def log(self, msg):
        line = f'{now()} {msg}'
        print(line, flush=True)
        with open(self.dir / 'orchestrator.log', 'a') as f:
            f.write(line + '\n')

    def save(self):
        if self.dry_run:
            return
        tmp = self.state_path.with_suffix('.tmp')
        tmp.write_text(json.dumps(self.state, indent=1))
        os.replace(tmp, self.state_path)

    def notify(self, subject, body=''):
        self.log(f'NOTIFY {subject} {body}')
        if self.cfg.get('notify_cmd'):
            self.sh(self.cfg['notify_cmd'].format(subject=shlex.quote(subject), body=shlex.quote(body)))

    # ---- helpers ------------------------------------------------------------ #
    def project(self, stage, action, mutating=True):
        return self.sh(['project.py', '--xml', str(self.xml), '--stage', stage, f'--{action}'], mutating=mutating)

    def write_xml(self):
        numjobs = {n: s['numjobs'] for n, s in self.state['stages'].items() if s.get('numjobs')}
        self.xml.write_text(render_xml(self.cfg, numjobs))

    def ensure_token(self):
        rc, out = self.sh(['htgettoken', '-a', 'htvaultprod.fnal.gov', '-i', self.cfg['experiment'], '--nooidc'],
                          mutating=False, timeout=120)
        if rc != 0:
            if not self.state.get('auth_alerted'):
                self.notify(f"[{self.cfg['name']}] token refresh needed",
                            f"run: htgettoken -a htvaultprod.fnal.gov -i {self.cfg['experiment']}")
                self.state['auth_alerted'] = True
            return False
        self.state['auth_alerted'] = False
        return True

    def cluster_states(self, clusters):
        """{cluster: list of job state letters still in the queue}."""
        out = {}
        for cl in clusters:
            rc, txt = self.sh(['jobsub_q', '-G', self.cfg['experiment'], '--jobid', cl], mutating=False, timeout=300)
            if rc != 0 and 'no matching' not in txt.lower() and txt.strip():
                return None  # query failed; try again next tick
            cid = cl.split('@')[0]
            states = []
            for line in txt.splitlines():
                tok = line.split()
                if tok and tok[0].startswith(cid + '.'):
                    st = next((t for t in tok[3:] if t in STATE_TOKENS), '?')
                    if st not in ('C', 'X'):
                        states.append(st)
            out[cl] = states
        return out

    def count_good(self, st):
        book = Path(stage_dirs(self.cfg, st)['bookdir'])
        lst = book / ('filesana.list' if st['ana'] else 'files.list')
        if not lst.exists():
            return 0
        return sum(1 for line in lst.read_text().splitlines() if line.strip())

    # ---- state machine -------------------------------------------------------- #
    def tick(self):
        stages = self.cfg['stages']
        for i, st in enumerate(stages):
            s = self.state['stages'][st['name']]
            if s['status'] == 'done':
                continue
            if s['status'] == 'failed':
                return 'failed'
            if s['status'] == 'pending':
                expected = self.cfg['gen_numjobs'] if i == 0 else self.state['stages'][stages[i - 1]['name']]['good']
                s['numjobs'] = s['expected'] = int(expected)
                self.write_xml()
                if not self.ensure_token():
                    return 'auth'
                return self.submit(st, s, 'submit')
            if s['status'] == 'running':
                return self.poll(st, s)
            if s['status'] == 'checking':
                return self.check(st, s)
        if not self.state.get('manifest'):
            self.write_manifest()
        for ext in self.cfg['external']:
            s = self.state['stages'][ext['name']]
            if s['status'] in ('pending', 'running'):
                return self.external(ext, s)
            if s['status'] == 'failed':
                return 'failed'
        if not self.state.get('finished'):
            self.state['finished'] = now()
            self.notify(f"[{self.cfg['name']}] campaign complete")
        return 'done'

    def submit(self, st, s, action):
        rc, out = self.project(st['name'], action)
        if rc != 0:
            s['errors'] = s.get('errors', 0) + 1
            self.log(f"{st['name']}: project.py --{action} failed (rc {rc}):\n{out[-2000:]}")
            if s['errors'] >= 3:
                s['status'] = 'failed'
                self.notify(f"[{self.cfg['name']}] {st['name']} submission failing", out[-500:])
            return 'error'
        clusters = sorted({f'{a}@{b}' for a, b in JOBID_RE.findall(out)})
        if not clusters and not self.dry_run:
            self.log(f"{st['name']}: no job id in project.py output; will check after "
                     f"{self.cfg['fallback_wait_hours']} h")
        s.update(status='running', clusters=clusters, submitted=self.clock(), held_since=None, errors=0)
        s['attempts'] = s.get('attempts', 0) + 1
        self.log(f"{st['name']}: --{action} attempt {s['attempts']} submitted {s['numjobs']} jobs {clusters}")
        return 'submitted'

    def poll(self, st, s):
        if not s['clusters']:
            if self.clock() - s['submitted'] < 3600 * self.cfg['fallback_wait_hours']:
                return 'waiting'
            s['status'] = 'checking'
            return self.check(st, s)
        q = self.cluster_states(s['clusters'])
        if q is None:
            return 'waiting'
        left = [x for v in q.values() for x in v]
        s['queue'] = {k: ''.join(sorted(v)) for k, v in q.items()}
        if not left:
            s['status'] = 'checking'
            return self.check(st, s)
        if all(x == 'H' for x in left):
            s['held_since'] = s.get('held_since') or self.clock()
            if self.clock() - s['held_since'] > 3600 * self.cfg['held_timeout_hours']:
                for cl in s['clusters']:
                    self.sh(['jobsub_rm', '-G', self.cfg['experiment'], '--jobid', cl])
                self.log(f"{st['name']}: removed {len(left)} jobs held > {self.cfg['held_timeout_hours']} h")
        else:
            s['held_since'] = None
        self.log(f"{st['name']}: {len(left)} jobs in queue "
                 f"({sum(x == 'R' for x in left)} R, {sum(x == 'I' for x in left)} I, {sum(x == 'H' for x in left)} H)")
        return 'running'

    def check(self, st, s):
        rc, out = self.project(st['name'], 'checkana' if st['ana'] else 'check', mutating=False)
        good = self.count_good(st)
        s['good'] = good
        exp = max(s['expected'], 1)
        self.log(f"{st['name']}: check rc {rc}: {good}/{s['expected']} good outputs")
        if good >= math.ceil(self.cfg['min_good_fraction'] * exp):
            return self.finish(st, s, 'done')
        if s['attempts'] <= self.cfg['max_makeup']:
            if not self.ensure_token():
                return 'auth'
            return self.submit(st, s, 'makeup')
        if good >= math.ceil(self.cfg['accept_fraction'] * exp):
            return self.finish(st, s, 'done', note=f'accepted below target ({good}/{s["expected"]})')
        return self.finish(st, s, 'failed', note=f'{good}/{s["expected"]} good after {s["attempts"]} attempts')

    def finish(self, st, s, status, note=''):
        s['status'], s['finished'] = status, now()
        if status == 'done':
            note = f'{note} {self.publish(st, s)}'.strip()
            if st['name'] == 'gen' and self.cfg.get('delta_sidecar', True):
                ok, msg = self.gen_sidecar(s)
                note = f'{note} {msg}'
                if not ok:
                    s['status'] = status = 'failed'
        self.notify(f"[{self.cfg['name']}] {st['name']} {status}", f"{s['good']}/{s['expected']} good. {note}".strip())
        # 'done' is reserved for the whole campaign; a finished stage advances.
        return 'advanced' if status == 'done' else status

    def publish(self, st, s):
        """Write the stage's file lists next to its outdir on /pnfs. Only the
        job directories of project.py's good outputs are searched, so files
        left by failed attempts are never published. Returns a short summary."""
        root = Path(stage_root(self.cfg, st))
        book = Path(stage_dirs(self.cfg, st)['bookdir'])
        src = book / ('filesana.list' if st['ana'] else 'files.list')
        good = [x.strip() for x in src.read_text().splitlines() if x.strip()] if src.exists() else []
        dirs = sorted({str(Path(f).parent) for f in good}) or [str(root / 'out')]
        lists = {'files': good}
        for p in st['publish']:
            found = set()
            for d in dirs:
                for f in Path(d).rglob(p['glob']):
                    if not (p.get('exclude') and fnmatch.fnmatch(f.name, p['exclude'])):
                        found.add(str(f))
            lists[p['name']] = sorted(found)
        published, warn = {}, []
        for name, files in lists.items():
            fname = 'files.list' if name == 'files' else f'{name}_files.list'
            local = self.dir / 'published' / st['name'] / fname
            local.parent.mkdir(parents=True, exist_ok=True)
            local.write_text(''.join(f + '\n' for f in files))
            copy_to_pnfs(self.sh, local, str(root / fname))
            published[name] = dict(list=str(root / fname), n=len(files))
            if name != 'files' and len(files) < len(good):
                warn.append(f'{name}: {len(files)} files for {len(good)} good jobs')
        s['published'] = published
        msg = 'published ' + ', '.join(f"{k}={v['n']}" for k, v in published.items()) + f' to {root}'
        if warn:
            msg += ' WARNING ' + '; '.join(warn)
        self.log(f"{st['name']}: {msg}")
        return msg

    def gen_sidecar(self, s):
        """Gen-level truth for the weights, and the gate on the patched table.
        Dumps every NC Δ decay of the good gen files (dump_delta_decays.C, read
        over xrootd), requires pilot_check.py to pass (radiative fraction vs W
        matches the enhanced table; join key unique), then writes the per-event
        weights and publishes both CSVs to <persistent_base>/gen_sidecar/.
        Returns (ok, message); a failure stops the campaign before g4."""
        pub = self.dir / 'published' / 'gen'
        good = [x for x in (pub / 'files.list').read_text().split() if x] if (pub / 'files.list').exists() else []
        prefix = self.cfg['xrootd_prefix']
        urls = [prefix + f[len('/pnfs/'):] if f.startswith('/pnfs/') else f for f in good]
        (pub / 'gen_xrootd.list').write_text(''.join(u + '\n' for u in urls))
        csv_path, w_path = pub / 'delta_decays.csv', pub / 'ncdeltarad_weights.csv'
        br = f"{self.cfg['br_gamma']:g}"
        rc, out = self.sh(['root', '-l', '-b', '-q',
                           f'{HERE / "dump_delta_decays.C"}("{pub / "gen_xrootd.list"}","{csv_path}")'],
                          mutating=False, timeout=6 * 3600)
        if rc != 0 and not self.dry_run:
            return False, f'gen sidecar: dump_delta_decays.C failed (rc {rc}): {out[-300:]}'
        rc, out = self.sh([sys.executable, str(HERE / 'pilot_check.py'), str(csv_path), '--br-gamma', br],
                          mutating=False)
        (pub / 'pilot_check.txt').write_text(out)
        if rc != 0 and not self.dry_run:
            return False, f'gen sidecar: pilot_check FAILED (see {pub / "pilot_check.txt"}):\n{out[-1500:]}'
        rc, out = self.sh([sys.executable, str(HERE / 'delta_br.py'), 'weights', '--in', str(csv_path),
                           '--out', str(w_path), '--br-gamma', br], mutating=False)
        if rc != 0 and not self.dry_run:
            return False, f'gen sidecar: weights failed: {out[-300:]}'
        dest = Path(manifest_dir(self.cfg)) / 'gen_sidecar'
        for src in (csv_path, w_path, pub / 'pilot_check.txt'):
            copy_to_pnfs(self.sh, src, str(dest / src.name))
        self.state['sidecar'] = dict(delta_decays=str(dest / csv_path.name), weights=str(dest / w_path.name),
                                     pilot_check=str(dest / 'pilot_check.txt'))
        self.log(f'gen sidecar: pilot_check passed; weights -> {dest / w_path.name}')
        return True, f'sidecar + weights published to {dest}'

    def write_manifest(self):
        """Campaign manifest on /pnfs: what was generated and where the outputs
        (and their lists) are, for the GPU SPINE step and for medulla."""
        dest = Path(manifest_dir(self.cfg))
        stages = {}
        for st in self.cfg['stages']:
            s = self.state['stages'][st['name']]
            stages[st['name']] = dict(fcl=st['fcl'], outdir=stage_dirs(self.cfg, st)['outdir'], good=s.get('good'),
                                      expected=s.get('expected'), attempts=s.get('attempts'),
                                      lists=s.get('published', {}))
        m = dict(name=self.cfg['name'], created=self.state['created'], finished=now(),
                 release=self.cfg['release'], qual=self.cfg['qual'], event_generator_list='NCRES',
                 overlay_inputdef=self.cfg['overlay_inputdef'], br_gamma=self.cfg['br_gamma'],
                 pdg_table=str(dest / 'ncdeltarad_pdg_table.txt'), config=str(dest / 'campaign_config.json'),
                 weight='per event from gen_sidecar/ncdeltarad_weights.csv, joined on (run, subrun, event); '
                        'computed by batch/ncdeltarad/delta_br.py from the GENIE Δ mass and decay channel',
                 sidecar=self.state.get('sidecar', {}),
                 stages=stages)
        local = self.dir / 'manifest.json'
        local.write_text(json.dumps(m, indent=1))
        for src, name in ((local, 'manifest.json'), (self.dir / 'ncdeltarad_pdg_table.txt', 'ncdeltarad_pdg_table.txt'),
                          (self.dir / 'config.json', 'campaign_config.json')):
            copy_to_pnfs(self.sh, src, str(dest / name))
        self.state['manifest'] = str(dest / 'manifest.json')
        self.notify(f"[{self.cfg['name']}] all production stages done", f"manifest: {self.state['manifest']}")

    def external(self, ext, s):
        fmt = dict(workdir=self.cfg['workdir'], pnfs_base=self.cfg['pnfs_base'],
                   persistent_base=self.cfg.get('persistent_base') or self.cfg['pnfs_base'], name=self.cfg['name'])
        if s['status'] == 'pending':
            rc, out = self.sh(ext['submit_cmd'].format(**fmt))
            if rc != 0:
                s['status'] = 'failed'
                self.notify(f"[{self.cfg['name']}] {ext['name']} submit failed", out[-500:])
                return 'failed'
            s.update(status='running', submitted=self.clock())
            return 'submitted'
        rc, out = self.sh(ext['check_cmd'].format(**fmt), mutating=False)
        if rc == 0:
            s['status'], s['finished'] = 'done', now()
            self.notify(f"[{self.cfg['name']}] {ext['name']} done")
            return 'advanced'
        if rc == 2:
            s['status'] = 'failed'
            self.notify(f"[{self.cfg['name']}] {ext['name']} failed", out[-500:])
            return 'failed'
        return 'running'


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def copy_to_pnfs(sh, src, dst):
    if sh.dry_run:
        sh.log(f'[dry-run] copy {src} -> {dst}')
        return
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    try:
        shutil.copyfile(src, dst)
    except OSError:
        rc, out = sh(['ifdh', 'cp', str(src), str(dst)])
        if rc != 0:
            raise SystemExit(f'cannot copy {src} -> {dst}: {out}')


def cmd_init(a):
    c = load_config(a.config)
    wd = Path(c['workdir'])
    if (wd / 'state.json').exists() and not a.force:
        raise SystemExit(f'{wd}/state.json exists (use --force to re-initialise)')
    wd.mkdir(parents=True, exist_ok=True)
    (wd / 'config.json').write_text(json.dumps(c, indent=1))
    sh = Shell(lambda m: print(m), a.dry_run)

    table_in = c['genie_pdg_table'] or os.path.join(os.environ.get('GENIE', ''),
                                                       'data/evgen/catalogues/pdg/genie_pdg_table.txt')
    text, n = delta_br.patch_table(open(table_in).read(), c['br_gamma'])
    if n != 2:
        raise SystemExit(f'{table_in}: patched {n} Nγ channels, expected 2')
    (wd / 'ncdeltarad_pdg_table.txt').write_text(text)
    fcldir = f"{c['pnfs_base']}/fcl"
    env = (HERE / 'ncdeltarad_genie_env.sh').read_text()
    env = env.replace('@NCDR_TABLE_PNFS@', f'{fcldir}/ncdeltarad_pdg_table.txt')
    env = env.replace('@NCDR_BR_GAMMA@', f"{c['br_gamma']:g}")
    (wd / 'ncdeltarad_genie_env.sh').write_text(env)
    gen_fcl = HERE / c['stages'][0]['fcl']
    if not gen_fcl.exists():
        raise SystemExit(f'{gen_fcl} not found (the gen fcl must live next to orchestrate.py)')
    for src in (wd / 'ncdeltarad_pdg_table.txt', wd / 'ncdeltarad_genie_env.sh', gen_fcl):
        copy_to_pnfs(sh, src, f'{fcldir}/{Path(src).name}')

    names = [s['name'] for s in c['stages']] + [e['name'] for e in c['external']]
    state = dict(created=now(), stages={n: dict(status='pending') for n in names})
    (wd / 'state.json').write_text(json.dumps(state, indent=1))
    Orchestrator(wd).write_xml()
    print(f'initialised {wd}: table {table_in} -> BR(Nγ)={c["br_gamma"]:g}; '
          f'XML {wd / (c["name"] + ".xml")}; staged fcl/env/table to {fcldir}')


def with_lock(wd, fn):
    with open(Path(wd) / '.lock', 'w') as lk:
        try:
            fcntl.flock(lk, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise SystemExit('another orchestrator instance holds the lock')
        return fn()


def cmd_tick(a):
    def go():
        o = Orchestrator(a.workdir, dry_run=a.dry_run)
        r = o.tick()
        o.save()
        return r
    print(with_lock(a.workdir, go))


def cmd_run(a):
    def go():
        while True:
            o = Orchestrator(a.workdir)
            r = o.tick()
            o.save()
            if r in ('done', 'failed'):
                return r
            # a fresh submission or check usually has follow-up work right away
            time.sleep(60 if r in ('submitted', 'advanced') else a.interval)
    print(with_lock(a.workdir, go))


def cmd_status(a):
    o = Orchestrator(a.workdir)
    print(f"campaign {o.cfg['name']}  (created {o.state['created']})")
    print(f"{'stage':10s} {'status':9s} {'att':>3s} {'good':>6s} {'expect':>6s}  queue / clusters")
    for n, s in o.state['stages'].items():
        q = ' '.join(f'{k}:{v}' for k, v in (s.get('queue') or {}).items()) or ' '.join(s.get('clusters', []))
        print(f"{n:10s} {s['status']:9s} {s.get('attempts', 0):3d} {s.get('good', ''):>6} "
              f"{s.get('expected', ''):>6}  {q}")


def cmd_retry(a):
    def go():
        o = Orchestrator(a.workdir)
        s = o.state['stages'][a.stage]
        s.update(status='pending', attempts=0, errors=0, clusters=[], queue={})
        o.save()
        o.log(f'{a.stage}: reset to pending by user')
    with_lock(a.workdir, go)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('init'); p.add_argument('config'); p.add_argument('--force', action='store_true')
    p.add_argument('--dry-run', action='store_true')
    p = sub.add_parser('status'); p.add_argument('workdir')
    p = sub.add_parser('tick'); p.add_argument('workdir'); p.add_argument('--dry-run', action='store_true')
    p = sub.add_parser('run'); p.add_argument('workdir'); p.add_argument('--interval', type=int, default=900)
    p = sub.add_parser('retry'); p.add_argument('workdir'); p.add_argument('stage')
    a = ap.parse_args()
    dict(init=cmd_init, status=cmd_status, tick=cmd_tick, run=cmd_run, retry=cmd_retry)[a.cmd](a)


if __name__ == '__main__':
    main()
