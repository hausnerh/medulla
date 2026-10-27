#!/usr/bin/env python3
"""
SPINE outputs -> ML flat CAFs on /pnfs, for a finished ncdeltarad campaign.
Run on the gpvm (AL9 host, like the orchestrator) after SPINE inference:

    python batch/ncdeltarad/spine_caf.py $WD --h5 '/exp/sbn/data/users/$USER/ncdr/spine/*.h5' --cafmaker-dir /exp/icarus/app/users/$USER/ml_cafmaker

For every caf job of the campaign (<persistent_base>/caf/out/<job>/):
  1. pairs the job with the SPINE HDF5 file holding its events: (run, subrun,
     event) of the job's flat CAF vs the HDF5 run_info (SPINE writes one HDF5
     per larcv file, and nothing ties the larcv name to the caf job);
  2. in the SL7 container (config sl7_exec), as sbn_ml_cafmaker's
     tools/manager.py does it:
         merge_sources_simulation tmp.root <job's .caf.root> <hdf5>
         flatten_caf tmp.root <caf stem>.spine.flat.caf.root
     then `ifdh cp` of the flat ML CAF and the HDF5 to
     <persistent_base>/spine/{flatcaf,hdf5}/;
  3. writes <persistent_base>/spine/{flatcaf,hdf5}_files.list and
     spine_manifest.json, and runs check_outputs.py (section 5 checks them).

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
        return set(zip(ri['run'].tolist(), ri['subrun'].tolist(), ri['event'].tolist()))


def flat_caf_keys(path):
    import uproot
    t = uproot.open(path)['recTree']
    a = t.arrays(['rec.hdr.run', 'rec.hdr.subrun', 'rec.hdr.evt'], library='np')
    return set(zip(a['rec.hdr.run'].tolist(), a['rec.hdr.subrun'].tolist(), a['rec.hdr.evt'].tolist()))


def pair_jobs(jobs, h5s):
    """jobs: {job dir: caf keys}; h5s: {h5 path: keys}. Returns
    [(job, h5, n_caf, n_common)] -- each job with the HDF5 sharing most of
    its events -- and the jobs with no HDF5 at all."""
    pairs, unpaired = [], []
    for job, keys in sorted(jobs.items()):
        best = max(h5s, key=lambda h: len(keys & h5s[h]), default=None)
        common = len(keys & h5s[best]) if best else 0
        if not common:
            unpaired.append(job)
        else:
            pairs.append((job, best, len(keys), common))
    return pairs, unpaired


def merge_script(setup, merge_exec, caf, h5, tmp, flat, dest_flat, dest_h5):
    """Bash for the SL7 container: merge, flatten, copy both files to /pnfs.
    dCache refuses to overwrite, so an existing destination is removed first."""
    q = shlex.quote
    body = [
        f'{q(merge_exec)} {q(tmp)} {q(caf)} {q(h5)} || {{ echo "spine_caf: merge FAILED" >&2; exit 4; }}',
        f'flatten_caf {q(tmp)} {q(flat)} || {{ echo "spine_caf: flatten_caf FAILED" >&2; exit 5; }}',
        f'rm -f {q(tmp)}',
        'command -v ifdh >/dev/null || setup ifdhc || { echo "spine_caf: no ifdh" >&2; exit 6; }',
    ]
    for src, dst in ((flat, dest_flat), (h5, dest_h5)):
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
    outdir = Path(O.stage_dirs(cfg, caf_st)['outdir'])
    setup = a.setup or setup_from_cmakecache(a.merge_exec)
    if not setup:
        raise SystemExit(f'cannot read the sbnana setup from {Path(a.merge_exec).parent}/CMakeCache.txt: '
                         f'pass --setup "{UPS_SETUP}; setup sbnana <ver> -q e26:prof"')
    setup = setup.replace('; ', '\n')
    print(f'merge: {a.merge_exec}\nsetup:\n  ' + setup.replace('\n', '\n  '))

    h5_files = sorted({f for p in a.h5 for f in (glob.glob(p) or [p])})
    h5s = {f: h5_keys(f) for f in h5_files}
    jobs, cafs = {}, {}
    for d in sorted(p for p in outdir.iterdir() if p.is_dir()):
        nested = [f for f in d.glob('*.caf.root') if not f.name.endswith('.flat.caf.root')]
        flat = [f for f in d.glob('*.root') if not f.name.endswith('.caf.root') and not f.name.startswith('Supplemental-')]
        if len(nested) != 1 or len(flat) != 1:
            print(f'skip {d.name}: {len(nested)} CAF / {len(flat)} flat CAF files')
            continue
        jobs[d.name], cafs[d.name] = flat_caf_keys(str(flat[0])), nested[0]
    pairs, unpaired = pair_jobs(jobs, h5s)
    for job, h5, n, common in pairs:
        print(f'{job}: {common}/{n} events in {os.path.basename(h5)} ({len(h5s[h5])} events)')
    for job in unpaired:
        print(f'{job}: NO SPINE output holds its events')
    used = {h for _, h, _, _ in pairs}
    for h in sorted(set(h5s) - used):
        print(f'unused HDF5: {h}')

    work = Path(a.workdir) / 'spine_caf'
    work.mkdir(exist_ok=True)
    dest = Path(base) / 'spine'
    flat_list, h5_list, done, failed = [], [], [], []
    for job, h5, n, common in pairs:
        stem = cafs[job].name[:-len('.caf.root')]
        flat_local = work / f'{stem}.spine.flat.caf.root'
        dest_flat, dest_h5 = dest / 'flatcaf' / flat_local.name, dest / 'hdf5' / os.path.basename(h5)
        script = work / f'merge_{job}.sh'
        script.write_text(merge_script(setup, a.merge_exec, str(cafs[job]), h5, str(work / f'{job}_tmp.root'),
                                       str(flat_local), str(dest_flat), str(dest_h5)))
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
        h5_list.append(str(dest_h5))
        done.append(dict(job=job, caf=str(cafs[job]), hdf5=str(dest_h5), flatcaf=str(dest_flat),
                         caf_events=n, paired_events=common, merge_matched=matched, merge_incomplete=incomplete))
        bad = matched < common or incomplete
        if bad:
            failed.append(job)
        print(f'{job}: merge matched {matched}/{n} CAF events, {incomplete} incomplete'
              f'{"  <-- PROBLEM, see " + str(work / f"merge_{job}.log") if bad else ""} -> {dest_flat}')
    if a.dry_run:
        return 0

    sh = O.Shell(print)
    for name, files in (('flatcaf', flat_list), ('hdf5', h5_list)):
        local = work / f'{name}_files.list'
        local.write_text(''.join(f + '\n' for f in files))
        O.copy_to_pnfs(sh, local, str(dest / local.name))
    man = work / 'spine_manifest.json'
    man.write_text(json.dumps(dict(merge_exec=a.merge_exec, setup=setup, jobs=done, failed=failed,
                                   unpaired=unpaired, unused_hdf5=sorted(set(h5s) - used)), indent=1))
    O.copy_to_pnfs(sh, man, str(dest / man.name))
    print(f'\n{len(done)} merged, {len(failed)} failed or incomplete, {len(unpaired)} without SPINE output -> {dest}')
    rc = subprocess.run([sys.executable, str(HERE / 'check_outputs.py'), base]).returncode
    return 1 if failed or unpaired or rc else 0


if __name__ == '__main__':
    sys.exit(main())
