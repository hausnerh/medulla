#!/usr/bin/env python3
"""
SPINE outputs -> ML flat CAFs on /pnfs, for a finished ncdeltarad campaign.
Run on the gpvm (AL9 host, like the orchestrator) after SPINE inference:

    python batch/ncdeltarad/spine_caf.py $WD --h5 '<pnfs_base>/spine/hdf5_incoming/pass*/*.h5' --cafmaker-dir /exp/icarus/app/users/$USER/ml_cafmaker

Run it once with --dry-run: it checks everything below and says how many jobs and
events would be merged. The lists are rewritten on every run, so a merge must be one
pass: nothing is merged unless the pairing is clean (each event in one HDF5 only,
each job's events all in one HDF5 of its own, no HDF5 left over).

For every good caf job of the campaign (the job directories of the stage's files.list):
  1. pairs the job with the SPINE HDF5 file holding its events: (run, subrun,
     event) of the job's flat CAF vs the HDF5 run_info (SPINE writes one HDF5
     per larcv file, and nothing ties the larcv name to the caf job);
  2. in the SL7 container (config sl7_exec), as sbn_ml_cafmaker's
     tools/manager.py does it:
         merge_sources_simulation tmp.root <job's .caf.root> <hdf5>
         flatten_caf tmp.root <caf stem>.spine.flat.caf.root
     then `ifdh cp` of the flat ML CAF to <persistent_base>/spine/flatcaf/ and of
     the HDF5 to <pnfs_base>/spine/hdf5/<job>_<name> (scratch), unless it is
     already under <pnfs_base>/spine/ (larcv names repeat across jobs, so HDF5
     names do too);
  3. writes <persistent_base>/spine/{flatcaf,hdf5}_files.list and
     spine_manifest.json, and runs check_outputs.py (section 5 checks them).
Jobs with no events have nothing to merge and are listed, not failed.

merge_sources_simulation must run against the sbnana/sbnanaobj it was built
with: the setup is read from the build's CMakeCache.txt unless --setup is
given. flatten_caf comes from the same sbnana.
"""
import argparse
import glob
import json
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import orchestrate as O  # noqa: E402

UPS_SETUP = 'source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh'
# The HDF5 files can sit on /pnfs; nothing writes them here, so no lock is needed.
os.environ.setdefault('HDF5_USE_FILE_LOCKING', 'FALSE')


def setup_from_cmakecache(exe):
    """'setup sbnana <ver> -q <qual>' for the sbnana the merge build used."""
    cache = Path(exe).parent / 'CMakeCache.txt'
    if not cache.exists():
        return None
    text = cache.read_text(errors='replace')
    ver = re.search(r'/sbnana/(v\d+_\d+_\d+\w*)/', text, re.I)
    qual = re.search(r'/sbnana/v[\w]+/[^/\s]*?(e\d+)[^/\s]*?(prof|debug)', text, re.I)
    if not ver:
        obj = sorted(set(re.findall(r'/sbnanaobj/(v\d+_\d+_\d+\w*)/', text, re.I)))
        if obj:
            print(f'{cache}: built against sbnanaobj {obj} but names no sbnana: pass --setup with the '
                  f'sbnana (for flatten_caf) that uses that sbnanaobj')
        return None
    q = f'{qual.group(1)}:{qual.group(2)}' if qual else 'e26:prof'
    return f'{UPS_SETUP}\nsetup sbnana {ver.group(1)} -q {q}'


def h5_keys(path):
    """{(run, subrun, event)} of a SPINE HDF5 file (its run_info product)."""
    import h5py
    with h5py.File(path, 'r') as f:
        if 'run_info' not in f:
            raise SystemExit(f'{path}: no run_info dataset (keys: {sorted(f.keys())[:20]})')
        ri = f['run_info'][()]
        keys = set(zip(ri['run'].tolist(), ri['subrun'].tolist(), ri['event'].tolist()))
        if len(keys) != len(ri):
            raise SystemExit(f'{path}: {len(ri) - len(keys)} repeated events in run_info')
        return keys


def flat_caf_keys(path):
    import uproot
    f = uproot.open(path)
    if 'recTree' not in f:
        return set()
    a = f['recTree'].arrays(['rec.hdr.run', 'rec.hdr.subrun', 'rec.hdr.evt'], library='np')
    return set(zip(a['rec.hdr.run'].tolist(), a['rec.hdr.subrun'].tolist(), a['rec.hdr.evt'].tolist()))


def pair_jobs(jobs, h5s):
    """jobs: {job dir: caf keys}; h5s: {h5 path: keys}. Returns
    [(job, h5, n_caf, n_common)] -- each job with the HDF5 sharing most of
    its events --, the jobs with no HDF5 at all, and the jobs with no events."""
    pairs, unpaired, empty = [], [], []
    for job, keys in sorted(jobs.items()):
        if not keys:
            empty.append(job)
            continue
        best = max(h5s, key=lambda h: len(keys & h5s[h]), default=None)
        common = len(keys & h5s[best]) if best else 0
        if not common:
            unpaired.append(job)
        else:
            pairs.append((job, best, len(keys), common))
    return pairs, unpaired, empty


def pairing_problems(pairs, unpaired, h5s):
    """Why these pairs must not be merged: events in two HDF5 files (double
    processing), jobs with events missing, HDF5 files with no or several jobs."""
    problems = [f'{job}: NO SPINE output holds its events' for job in unpaired]
    seen = {}
    for h5, keys in sorted(h5s.items()):
        for k in sorted(keys):
            if k in seen:
                problems.append(f'event {k} is in both {seen[k]} and {h5}')
            seen.setdefault(k, h5)
    users = {}
    for job, h5, n, common in pairs:
        users.setdefault(h5, []).append(job)
        if common != n:
            problems.append(f'{job}: only {common}/{n} events in {h5}')
    problems += [f'{h5} holds events of several jobs: {jobs}' for h5, jobs in sorted(users.items()) if len(jobs) > 1]
    problems += [f'unused HDF5: {h5}' for h5 in sorted(set(h5s) - set(users))]
    return problems


def good_job_dirs(cfg, st):
    """Job directories of project.py's good outputs, as publish() finds them:
    files left in out/ by failed attempts are never read."""
    src = Path(O.stage_dirs(cfg, st)['bookdir']) / ('filesana.list' if st.get('ana') else 'files.list')
    good = [x.strip() for x in src.read_text().splitlines() if x.strip()] if src.exists() else []
    if not good:
        raise SystemExit(f'{src}: missing or empty, is the {st["name"]} stage done?')
    return sorted({Path(f).parent for f in good})


def h5_destination(pnfs_base, h5, job):
    """Where the HDF5 goes on scratch, or None if it is there already (the merge
    removes a destination before copying, so a file must never be its own)."""
    spine = Path(pnfs_base) / 'spine'
    if spine in Path(h5).parents:
        return None
    return str(spine / 'hdf5' / f'{job}_{os.path.basename(h5)}')


def merge_script(setup, merge_exec, caf, h5, tmp, flat, dest_flat, dest_h5):
    """Bash for the SL7 container: merge, flatten, copy the files to /pnfs (the
    HDF5 only if dest_h5 is given). dCache refuses to overwrite, so an existing
    destination is removed first."""
    q = shlex.quote
    body = [
        'export HDF5_USE_FILE_LOCKING=FALSE',
        f'{q(merge_exec)} {q(tmp)} {q(caf)} {q(h5)} || {{ echo "spine_caf: merge FAILED" >&2; exit 4; }}',
        f'flatten_caf {q(tmp)} {q(flat)} || {{ echo "spine_caf: flatten_caf FAILED" >&2; exit 5; }}',
        f'rm -f {q(tmp)}',
        'command -v ifdh >/dev/null || setup ifdhc || { echo "spine_caf: no ifdh" >&2; exit 6; }',
    ]
    for src, dst in ((flat, dest_flat), (h5, dest_h5)) if dest_h5 else ((flat, dest_flat),):
        body += [f'ifdh rm {q(dst)} >/dev/null 2>&1',
                 f'ifdh cp {q(src)} {q(dst)} || {{ echo "spine_caf: ifdh cp {src} FAILED" >&2; exit 6; }}']
    script = O.sl7_script(setup, ['true'])
    return script.replace('true\n', '\n'.join(body) + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('workdir', help='campaign workdir (config.json, state.json)')
    ap.add_argument('--h5', required=True, nargs='+', help='SPINE HDF5 files or globs (local paths)')
    ap.add_argument('--cafmaker-dir', help='prefix given to build_ml_cafmaker.sh: its merge build and setup')
    ap.add_argument('--merge-exec', help='sbn_ml_cafmaker build/merge_sources_simulation (another build)')
    ap.add_argument('--setup', help='UPS setup lines for merge + flatten_caf (default: from the build CMakeCache)')
    ap.add_argument('--dry-run', action='store_true', help='pair and write the scripts, run nothing')
    a = ap.parse_args()
    if a.cafmaker_dir:
        d = Path(a.cafmaker_dir)
        a.merge_exec = a.merge_exec or str(d / 'sbn_ml_cafmaker' / 'build' / 'merge_sources_simulation')
        if not a.setup and (d / 'setup_ml_cafmaker.sh').exists():
            a.setup = (d / 'setup_ml_cafmaker.sh').read_text().strip()
    if not a.merge_exec:
        ap.error('give --cafmaker-dir (build_ml_cafmaker.sh prefix) or --merge-exec')

    cfg = json.loads((Path(a.workdir) / 'config.json').read_text())
    base = cfg.get('persistent_base') or cfg['pnfs_base']
    caf_st = next(s for s in cfg['stages'] if s['name'] == 'caf')
    setup = a.setup or setup_from_cmakecache(a.merge_exec)
    if not setup:
        raise SystemExit(f'cannot read the sbnana setup from {Path(a.merge_exec).parent}/CMakeCache.txt: '
                         f'pass --setup "{UPS_SETUP}; setup sbnana <ver> -q e26:prof"')
    setup = setup.replace('; ', '\n')
    print(f'merge: {a.merge_exec}\nsetup:\n  ' + setup.replace('\n', '\n  '))

    h5_files = sorted({f for p in a.h5 for f in (glob.glob(p) or [p])})
    h5s = {f: h5_keys(f) for f in h5_files}
    jobs, cafs, problems = {}, {}, []
    for d in good_job_dirs(cfg, caf_st):
        nested = [f for f in d.glob('*.caf.root') if not f.name.endswith('.flat.caf.root')]
        flat = [f for f in d.glob('*.root') if not f.name.endswith('.caf.root') and not f.name.startswith('Supplemental-')]
        if len(nested) != 1 or len(flat) != 1:
            problems.append(f'{d.name}: {len(nested)} CAF / {len(flat)} flat CAF files')
            continue
        jobs[d.name], cafs[d.name] = flat_caf_keys(str(flat[0])), nested[0]
    pairs, unpaired, empty = pair_jobs(jobs, h5s)
    for job, h5, n, common in pairs:
        print(f'{job}: {common}/{n} events in {os.path.basename(h5)} ({len(h5s[h5])} events)')
    for job in empty:
        print(f'{job}: no events, nothing to merge')
    print(f'{len(pairs)} jobs with {sum(n for _, _, n, _ in pairs)} CAF events, {len(h5s)} HDF5 files with '
          f'{sum(len(k) for k in h5s.values())} events, {len(empty)} jobs without events')
    problems += pairing_problems(pairs, unpaired, h5s)
    if problems:
        print('\n'.join(problems[:50]))
        print(f'{len(problems)} problems: nothing merged')
        return 1

    work = Path(a.workdir) / 'spine_caf'
    work.mkdir(exist_ok=True)
    dest = Path(base) / 'spine'
    flat_list, h5_list, done, failed = [], [], [], []
    for job, h5, n, common in pairs:
        stem = cafs[job].name[:-len('.caf.root')]
        flat_local = work / f'{stem}.spine.flat.caf.root'
        # Only the flat ML CAF is kept long term; the HDF5 can be remade from larcv.
        dest_flat = dest / 'flatcaf' / flat_local.name
        dest_h5 = h5_destination(cfg['pnfs_base'], h5, job)
        script = work / f'merge_{job}.sh'
        script.write_text(merge_script(setup, a.merge_exec, str(cafs[job]), h5, str(work / f'{job}_tmp.root'),
                                       str(flat_local), str(dest_flat), dest_h5))
        if a.dry_run:
            print(f'[dry-run] {script}')
            continue
        cmd = cfg['sl7_exec'].format(script=shlex.quote(str(script))) if cfg.get('sl7_exec') else f'bash {script}'
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True, stdin=subprocess.DEVNULL)
        (work / f'merge_{job}.log').write_text(p.stdout + p.stderr)
        if p.returncode != 0:
            failed.append(job)
            tail = '\n    '.join((p.stdout + p.stderr).strip().splitlines()[-6:])
            print(f'{job}: FAILED (rc {p.returncode}), log {work / f"merge_{job}.log"}:\n    {tail}')
            continue
        flat_local.unlink(missing_ok=True)
        # merge_sources prints one "Matched Event" per CAF event found in the
        # HDF5 and "Found incomplete entry" when reading its products fails
        # (the event then keeps no SPINE products; rc is still 0).
        out = p.stdout + p.stderr
        matched, incomplete = out.count('Matched Event'), out.count('Found incomplete entry')
        flat_list.append(str(dest_flat))
        h5_list.append(dest_h5 or h5)
        done.append(dict(job=job, caf=str(cafs[job]), hdf5=dest_h5 or h5, flatcaf=str(dest_flat),
                         caf_events=n, paired_events=common, merge_matched=matched, merge_incomplete=incomplete))
        # Same StandardRecord class version as the input CAF but another
        # layout: ROOT keeps the input's layout and drops dlp/dlp_true from the
        # output (pilot, sbnanaobj v10_00_09 on v10_00_05_03 CAFs).
        clash = 'as the active class but a different checksum' in out
        bad = matched < common or incomplete or clash
        if bad:
            failed.append(job)
        print(f'{job}: merge matched {matched}/{n} CAF events, {incomplete} incomplete'
              + (', StandardRecord version clash with the input CAF (SPINE products dropped; rebuild '
                 'with build_ml_cafmaker.sh)' if clash else '')
              + (f'  <-- PROBLEM, see {work / f"merge_{job}.log"}' if bad else '') + f' -> {dest_flat}')
    if a.dry_run:
        return 0

    sh = O.Shell(print)
    for name, files in (('flatcaf', flat_list), ('hdf5', h5_list)):
        local = work / f'{name}_files.list'
        local.write_text(''.join(f + '\n' for f in files))
        O.copy_to_pnfs(sh, local, str(dest / local.name))
    man = work / 'spine_manifest.json'
    man.write_text(json.dumps(dict(merge_exec=a.merge_exec, setup=setup, jobs=done, failed=failed,
                                   empty=empty), indent=1))
    O.copy_to_pnfs(sh, man, str(dest / man.name))
    print(f'\n{len(done)} merged, {len(failed)} failed or incomplete, {len(empty)} without events -> {dest}')
    rc = subprocess.run([sys.executable, str(HERE / 'check_outputs.py'), base]).returncode
    return 1 if failed or rc else 0


if __name__ == '__main__':
    sys.exit(main())
