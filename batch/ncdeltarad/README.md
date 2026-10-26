# Enhanced-radiative NC Δ sample: grid campaign

A dedicated BNB NC Δ(1232) overlay sample in which Δ→Nγ and Δ→Nπ⁰ are produced
1:1 at the Δ pole (nominal is 0.6% : 66%). Each event is weighted back to
nominal. The sample is for re-tuning the gOre 1γ1p cuts, which are limited by
signal MC statistics: there are 134 raw stage1 signal events in the run2 CV,
and the old cvext had the same rate per POT.

The chain is the ReCAF2026 overlay production at icaruscode `v10_06_00_13`. Only
the generator stage changes:

| Stage | fcl |
|---|---|
| decode | `decoder_run2_icarus.fcl` (decodes the overlay data's raw fragments: `daqTPCROI:*RAW`, `daqPMT`, `pmtbaselines`, `crthit`, `daqTrigger`) |
| gen | `simulation_genie_icarus_bnb_overlays_volDetEnclosure_fiducial_ncdeltarad.fcl` (this directory) |
| g4 | `larg4_icarus_cosmics_sce.fcl` |
| detsim | `detsim_2d_icarus_refactored_yzsim_overlay_Run4.fcl` |
| overlay | `overlay_waveforms.fcl` (merges MC into the data waveforms: `overlayTPCRaw*`, `overlayOpWaveforms`, `overlayCRTHit`) |
| stage0 | `stage0_run2_wcdnn_icarus_overlay.fcl` |
| stage1 | `stage1_run2_larcv_icarus_overlay.fcl` (its larcv output feeds SPINE) |
| caf | `cafmakerjob_icarus_detsim2d_systtools_and_fluxwgt_overlay.fcl` |

The overlay stage is icaruscode's `standard_overlay_overlay_icarus.fcl` step
(`fcl/standard/overlay/`). Without it, stage0 fails on the first event with
`ProductNotFound ... overlayOpWaveforms` (pilot).

The run2 overlay input (`..._offbeambnbminbias_compressed_data`) is compressed
*raw* data. It holds only artdaq fragments (`DAQEVB04` and
`ICARUSReprocessRaw`). The production overlay order is decode → gen → g4 →
detsim → overlay → stage0 → stage1, so the decoder runs first on the data
files and gen overlays onto its output. Without it, overlay fails with
`Event failed to find raw::RawDigit with label daqTPCROI:PHYSCRATEDATATPCWWRAW`
(pilot event dumps). The first stage, whatever its name, reads
`overlay_inputdef` with `gen_numjobs` jobs. Every stage up to and including
gen must take one file per job.

The SPINE step (reconstruction on the larcv files, then merging into the CAFs)
needs a GPU node, so it is not part of this campaign. The campaign's product
is the stage outputs on `/pnfs`, with file lists and a manifest (section 3a)
that the SPINE step reads.

## 1. Physics setup

**Process.** The gen fcl includes the production one and sets
`EventGeneratorList: "NCRES"`.
- `NCRES` is GENIE's global list (`EventGeneratorListAssembler.xml`). The
  AR23_20i tune still supplies the RES-NC model.
- POT stays per generated spill, so the sample is normalised by its own POT.
- The source is still the overlay `RootInput`, so every overlay event gets one
  NC RES interaction.

**Branching ratios.** GENIE R-3_04_02 has no fcl or XML knob for resonance
branching ratios: `BaryonResonanceDecayer::InhibitDecay` is a stub. It does
read `$GENIE_PDG_TABLE` before its default table (`PDGLibrary::LoadDBase`).
- `delta_br.py table` copies `genie_pdg_table.txt` and changes only two
  numbers: BR(Δ⁺→pγ) = BR(Δ⁰→nγ) = 0.663, which equals BR(Nπ⁰).
- `ncdeltarad_genie_env.sh` (the gen stage's `<initsource>`) exports
  `GENIE_PDG_TABLE`. It refuses to run unless the patched lines are present.
- **Only the gen stage may set it.** The CAF stage's GENIE reweighting reads
  the same table.

**The split depends on W.** GENIE re-evaluates the Δ⁺/Δ⁰ branching ratios at
every decay from the Δ invariant mass W (`EvolveDeltaBR`). The phase space
scales as p³ for Nπ and E_γ³·F_γ² for Nγ, and a channel closes below its
threshold. The resulting probabilities (`delta_br.py split`):

| W (GeV) | 1.10 | 1.15 | 1.233 | 1.30 | 1.40 | 1.60 |
|---|---|---|---|---|---|---|
| nominal Nγ | 0.034 | 0.011 | 0.006 | 0.0045 | 0.0032 | 0.0017 |
| enhanced Nγ | 0.79 | 0.55 | 0.40 | 0.33 | 0.26 | 0.16 |
| enhanced Nπ⁰ | 0.15 | 0.31 | 0.40 | 0.44 | 0.49 | 0.56 |
| enhanced Nπ± | 0.06 | 0.14 | 0.20 | 0.22 | 0.25 | 0.28 |

**Weight back to nominal.** The Nπ⁰ and Nπ± branching ratios keep their
nominal values, so both pion channels scale by the same factor. Every
non-radiative decay therefore gets one common weight:

```
radiative:      w = P_γ^nom(W) / P_γ^enh(W)            (≈ 0.015 at the pole)
non-radiative:  w = (1 − P_γ^nom(W)) / (1 − P_γ^enh(W))  (≈ 1.66 at the pole; ~25 near the Nπ threshold, where few events sit)
```

**Where the weight's inputs come from: the gen-level sidecar, not the CAF.**
W must be the Δ invariant mass GENIE decayed with. CAF truth does not carry
it. `SRTrueInteraction::w` is nugen's `MCNeutrino::W`, computed for a free
nucleon at rest (`GENIE2ART.cxx`). The pilot measured it against the Δ mass:
median +31 MeV, 1–99% range [−0.15, +0.49] GeV. GENIE's selected W
(`GTruth::fgW`) is not in the CAF either.

So after the gen stage, the orchestrator runs `dump_delta_decays.C` over the
good gen files. That writes one row per NC Δ decay, with (run, subrun,
event), the Δ PDG code, the Δ mass and the decay channel. It then runs
`pilot_check.py` as a gate and writes `ncdeltarad_weights.csv`
(`delta_br.py weights`). Medulla joins that file on (Run, Subrun, Evt). The
art event ID is unchanged from gen through CAF.

With the exact channel available, a strict π⁰:γ = 50/50 (scaling Nπ±
differently) would also be possible. It needs a per-channel weight; ask if
you want it.

**Using it in the analysis.**
- The sample is a complete NC Δ(1232) sample: all three channels, reweighted.
  Its `resnum == 0` events replace the CV's NC `resnum == 0` events. Drop the
  higher resonances it also contains; the CV covers those.
- Normalise it by its own POT and multiply by the sidecar weight.
- Statistics gain for radiative events is about 1/w: ×63 per POT, measured on
  the pilot.

**Geometry scan.** The gen fcl also sets `GeomScan: "flux: 10000 2.0 1"`.
- **Why:** GENIEHelper forces one global probability scale, built from the
  maximum cross-section times the maximum path length per material. The
  RES-NC cross-section flattens at high energy. Unlike with the Default list,
  that leaves no margin for the default box scan's underestimated path
  lengths. The pilot crashed with GENIE FATAL
  `Negative no interaction probability! (P = -40%)`.
- **Cost:** only GENIE throws. POT uses the same scale, so normalisation is
  unaffected.

## 2. Files

| File | Purpose |
|---|---|
| `simulation_genie_..._ncdeltarad.fcl` | gen fcl (NCRES only) |
| `delta_br.py` | GENIE Δ branching-ratio emulation, table patcher, `weight()`; `split` prints the table above; `weights` turns the sidecar into per-event weights |
| `ncdeltarad_genie_env.sh` | gen-stage `<initsource>`: sets and validates `GENIE_PDG_TABLE` |
| `campaign_run2.toml` | campaign config (paths, overlay SAM definition, stages, resources) |
| `orchestrate.py` | campaign orchestrator around larbatch `project.py` |
| `dump_delta_decays.C` | gallery macro: every NC Δ decay in gen files → sidecar CSV |
| `pilot_check.py` | sidecar validation: the gen-stage gate (also runnable by hand) |
| `../test_ncdeltarad.py` | tests: branching-ratio emulation, weight closure, patcher, pilot check, orchestrator state machine, /pnfs publishing + manifest |

## 3. Running

**Run the orchestrator on the AL9 host, not in the SL7 container.** Inside
the container, `jobsub_submit` fails with `HTCondorIOError: Failed
communication with collector`. larbatch is plain Python and runs on AL9 from
cvmfs. It also needs ICARUS's `experiment_utilities` module on `PYTHONPATH`.
The gen-stage sidecar dump still needs SL7 for ROOT and gallery, so it runs
in the container through `sl7_exec` (set in `campaign_run2.toml`).

`project.py` on AL9 needs three cvmfs directories on its paths: larbatch,
ICARUS's `experiment_utilities` (icarusutil `lib/`), and the SAM Python client.
`path_prepend` / `pythonpath_prepend` in `campaign_run2.toml` hold them. The
orchestrator applies them to every command it runs, including under tmux or
cron. Verified on icarusgpvm01 (2026-10-20): `project.py --status` works with
exactly those three directories. To run `project.py` by hand, export the same
directories (directories, not `.py` files).

```bash
# on the AL9 host
htgettoken -a htvaultprod.fnal.gov -i icarus
# init needs no release when genie_pdg_table points at the cvmfs table:
#   genie_pdg_table = "/cvmfs/larsoft.opensciencegrid.org/products/genie/v3_04_02a/Linux64bit+3.10-2.17-e26-prof/GENIE-Generator/data/evgen/catalogues/pdg/genie_pdg_table.txt"

cp batch/ncdeltarad/campaign_run2.toml my_campaign.toml   # fill in the FILL_ME values
python batch/ncdeltarad/orchestrate.py init my_campaign.toml
python batch/ncdeltarad/orchestrate.py tick <workdir> --dry-run   # print what it would submit
tmux new -s ncdr 'python batch/ncdeltarad/orchestrate.py run <workdir>'
python batch/ncdeltarad/orchestrate.py status <workdir>
```

**What `init` does.**
1. Patches the table from `$GENIE` (or `genie_pdg_table`).
2. Fills in and stages the env script, table and gen fcl to `<pnfs_base>/fcl`.
   That directory is the project's `<fcldir>`, and the table also ships with
   the job via `jobsub -f`.
3. Writes `state.json` and the project.py XML.

**What each tick does, per stage, in order:**
1. Re-renders the XML. `numjobs` equals the previous stage's good outputs, so
   no input file is silently dropped.
2. Refreshes the token non-interactively. If that fails, it waits and sends
   one notification.
3. Submits with `project.py --submit` and records the jobsub cluster IDs. A
   fresh submit (new campaign, `init --force`, `retry`) first runs
   `project.py --clean_one` for that stage only, because project.py refuses to
   submit into non-empty out/log/book directories. Makeups keep their outputs.
4. Refreshes the token again, then polls `jobsub_q`. A failed query is logged
   with its error (`<stage>: jobsub_q failed (N in a row)`), and the third
   failure in a row sends a notification. Commands run with stdin closed, so an
   auth prompt fails at once instead of hanging the loop until the 300 s
   timeout. Any job held for longer than `held_timeout_hours` is
   removed with `jobsub_rm`, even while other jobs of the stage run, and the
   makeup redoes it. A job held for exceeding its memory or disk request never
   runs: raise the request in the TOML and run `reconfig` before the makeup.
5. When the queue is empty, runs `--check` (`--checkana` for the CAF stage)
   and counts `files.list` / `filesana.list`.
6. Decides:
   - at least `min_good_fraction` good: done;
   - otherwise `--makeup`, up to `max_makeup` times;
   - then accept if at least `accept_fraction` is good, else the stage
     fails and the campaign stops.

7. Publishes the finished stage's file lists to `/pnfs` (section 3a).
8. **After gen only:** writes the sidecar and gates on it. It runs
   `dump_delta_decays.C` over xrootd, then `pilot_check.py`, then the
   weights. If the radiative fraction does not match the patched table, or
   the join key is not unique, the gen stage is marked failed and nothing is
   submitted to g4.

**Operations.**
- Optional external stages after `caf` are `submit_cmd`/`check_cmd` shell
  commands in the config. `check_cmd` exits 0 for done, 1 for running and 2
  for failed. None is configured: SPINE runs separately on a GPU node.
- State is written atomically, so `run` can be killed and restarted. `tick`
  also works from cron. `run` holds `.runlock` for its lifetime (one loop per
  workdir) but `.lock` only during each tick, so `reconfig`, `retry`, `adopt`
  and `recheck` work while the loop runs; they wait for a tick in progress.
- `retry <stage>` starts a stage over: its fresh submit cleans the stage's
  outputs first.
- `recheck <stage>` re-runs the check of a stage whose jobs already ran, with
  no resubmission. Use it when the check itself failed.
- `republish <stage>` recounts a done stage and rewrites its published lists
  and the manifest, with no jobs. Use it after `reconfig` changed the stage's
  `publish` globs (pilot: the grid renamed the flat CAFs to `<input>_<uuid>.root`).
- `check_outputs.py <persistent_base>` checks a finished campaign end to end:
  lists, CAF census, CAF truth against the gen sidecar, larcv entries.
- `adopt <stage> <jobid>` records a cluster that was submitted outside the
  loop as the stage's current attempt, so the next tick polls it instead
  of resubmitting.
- `reconfig <toml>` reloads the config into `config.json` and the XML without
  resetting the state. It also accepts inserted stages: the new stage and
  every stage after it go back to pending, while earlier stages keep their
  state. Removing or reordering stages needs `init --force`, which resets
  everything.
- `init`, `run` and `reconfig` first check that `project.py -h` starts with
  the configured `path_prepend` / `pythonpath_prepend`. A loop started from a
  shell without the manual exports otherwise fails every check
  (`ModuleNotFoundError: larbatch_posix`).
- `notify_cmd` (e.g. `mail`) fires on stage done/failed and when the token
  needs a refresh.

## 3a. What ends up on /pnfs

Intermediate stages stay on scratch. `stage1` and `caf` (`persistent = true`)
plus the manifest go to `persistent_base`:

```
<pnfs_base>/fcl/                     gen fcl, env script, patched table (job inputs)
<pnfs_base>/{decode,gen,g4,detsim,overlay,stage0}/
    out/<job>/...                    art outputs
    files.list                       good outputs (project.py --check)
<persistent_base>/stage1/
    out/<job>/...                    art + larcv outputs
    files.list                       good art outputs
    larcv_files.list                 larcv files of the good jobs      <- SPINE input
<persistent_base>/caf/
    out/<job>/...
    files.list                       good outputs (project.py --checkana)
    flatcaf_files.list               *.flat.caf.root                   <- merge target / medulla
    caf_files.list                   *.caf.root (non-flat)
<persistent_base>/gen_sidecar/
    delta_decays.csv                 every NC Δ decay: run, subrun, event, Δ pdg, Δ mass, channel
    ncdeltarad_weights.csv           per-event weight back to nominal   <- join on (run, subrun, event)
    pilot_check.txt                  the gate's output
<persistent_base>/manifest.json      release, NCRES, BR(Nγ), overlay definition, and per stage:
                                     outdir, good/expected, attempts, list paths and counts
<persistent_base>/ncdeltarad_pdg_table.txt, campaign_config.json   (provenance)
```

**How the lists are built.**
- The publish lists search only the output directories of the jobs that
  project.py counted as good. Leftovers of failed attempts are never listed.
- If a list has fewer files than there are good jobs (e.g. a larcv file is
  missing), the stage notification says so.
- The globs are config entries (`[[stage]] publish`). The larcv name comes
  from `superaMC` with `unique_filename`; confirm it on the pilot. Scratch is
  purged LRU, so set `persistent = true` on any other stage you want to keep.

## 4. Pilot first (`gen_numjobs = 10`), then scale

The orchestrator now runs the pilot check itself after gen (step 8 above).
By hand, on a local gen file:

```bash
echo $PWD/gen.root > gen.list
root -l -b -q 'batch/ncdeltarad/dump_delta_decays.C("gen.list","delta_decays.csv")'
python batch/ncdeltarad/pilot_check.py delta_decays.csv
```

**Required checks:**
- The radiative fraction per W bin matches the enhanced prediction. This
  also catches a `TDatabasePDG` filled before GENIE loaded its table.
- The join key is unique and every Δ has a two-body decay.

**Informational checks:** CAF `w` vs the Δ mass, the 30 MeV photon tag, and
struck nucleon vs Δ charge.

**Interactive smoke test (17 overlay files, 641 events), 2026-10-20.**
- **Enhanced run:** 537 NC Δ decays, χ²/ndf = 2.1/6 against the enhanced
  prediction (e.g. 0.458 ± 0.030 observed vs 0.425 predicted at the pole).
- **Nominal control:** χ²/ndf = 361/6, with about 0.8% radiative at the pole.
- **Photon tag:** exact (0% mis-tagged).
- **CAF `w`:** not the Δ mass (see section 1).

**Sizing.**
- **Overlay files:** the run2 off-beam min-bias files average about 38
  events. One file had only 9.
- **Signal yield (rough):** about 0.05 stage1 signal events per generated
  event, so about 1.9 per file.
- **Target:** ten times the CV's 134 stage1 signal events needs about 700
  files, i.e. about 700 gen jobs. Each gen job must take **one** overlay file
  (`maxfilesperjob = 1`): every file is its own art run, and larg4 aborts at the
  first event of a second run in the same job (pilot, 2026-10-22: G4Exception
  Tracking0011, "ProcessManager is NULL for particle = proton"). `init` refuses
  more unless `allow_multifile_gen = true`.
- **Measure, don't trust the estimate:** count stage1 category-0 events per
  job on the pilot.
- **Tape:** only about 1 in 3.5 files of the definition was on disk, so
  prestage the subset first (`samweb prestage-dataset`).

## 5. Not verified here (check on the pilot)

These were taken from the larbatch/jobsub conventions and could not be run
from a laptop:

- **`<initsource>`:** that it is sourced *after* the release setup, so `ifdh`
  is available. The table also arrives via `jobsub -f` into
  `$CONDOR_DIR_INPUT`, which the env script uses first.
- **`<ana>1</ana>` and `--checkana` for the CAF stage**, and the
  `filesana.list` name.
- **The job ID line `Use job id N.0@jobsubNN.fnal.gov`** in `project.py --submit`
  output. If it is missing, the orchestrator waits `fallback_wait_hours` and
  then checks.
- **The container image, memory/disk/lifetime and resource requests** in
  `campaign_run2.toml`. These are first guesses: take them from the
  production XMLs for the same fcls. Measured in the pilot: detsim and stage0
  peak at ~16 GB, so they request 20000 MB (detsim was held at 6000).
- **The overlay SAM definition** that matches the run2 production.

## 6. Still needed in medulla

- **The weight join:** join `gen_sidecar/ncdeltarad_weights.csv` on (Run,
  Subrun, Evt). The medulla trees already carry these.
- **A per-event weight in spineplot.** Today it only applies one POT scale per
  sample. `batch/gore_cut_opt.py` already takes explicit weights, so the cut
  re-tune can use the sample first.
- **The replacement logic:** drop CV NC `resnum == 0` events and add this
  sample's `resnum == 0` events with the sidecar weight.
