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
| gen | `simulation_genie_icarus_bnb_overlays_volDetEnclosure_fiducial_ncdeltarad.fcl` (this directory) |
| g4 | `larg4_icarus_cosmics_sce.fcl` |
| detsim | `detsim_2d_icarus_refactored_yzsim_overlay_Run4.fcl` |
| stage0 | `stage0_run2_wcdnn_icarus_overlay.fcl` |
| stage1 | `stage1_run2_larcv_icarus_overlay.fcl` (its larcv output feeds SPINE) |
| caf | `cafmakerjob_icarus_detsim2d_systtools_and_fluxwgt_overlay.fcl` |
| spine | your existing SPINE + merge-into-CAF step |

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

The weight needs only three inputs, all in CAF truth (`SRTrueInteraction`):

| Input | Source |
|---|---|
| W | `w` |
| Δ charge | `hitnuc`: in NC a struck proton gives Δ⁺, a struck neutron Δ⁰ |
| radiative or not | a primary photon above 30 MeV in an NC `resnum == 0` event |

The radiative tag works because GENIE decays the Δ before final-state
interactions, the photon is untouched by them, and π⁰ is still undecayed at
the GENIE level. Making π⁰ and γ exactly 50/50 would mean scaling the Nπ±
channel differently. The weight would then depend on π⁰ versus π± at
production, which final-state interactions hide from the CAF, so it would
need a per-event record carried from gen.

**Using it in the analysis.**
- The sample is a complete NC Δ(1232) sample: all three channels, reweighted.
  Its `resnum == 0` events replace the CV's NC `resnum == 0` events. Drop the
  higher resonances it also contains; the CV covers those.
- Normalise it by its own POT and multiply by w.
- Statistics gain for radiative events is about 1/w: ~65× per POT at the pole.

## 2. Files

| File | Purpose |
|---|---|
| `simulation_genie_..._ncdeltarad.fcl` | gen fcl (NCRES only) |
| `delta_br.py` | GENIE Δ branching-ratio emulation, table patcher, `weight()`; `split` prints the table above |
| `ncdeltarad_genie_env.sh` | gen-stage `<initsource>`: sets and validates `GENIE_PDG_TABLE` |
| `campaign_run2.toml` | campaign config (paths, overlay SAM definition, stages, resources) |
| `orchestrate.py` | campaign orchestrator around larbatch `project.py` |
| `dump_delta_decays.C` | gallery macro: every NC Δ decay in gen files → CSV |
| `pilot_check.py` | pilot validation from that CSV |
| `../test_ncdeltarad.py` | tests: branching-ratio emulation, weight closure, patcher, pilot check, orchestrator state machine |

## 3. Running

On an ICARUS gpvm, inside the SL7 container used for e26 releases:

```bash
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup icaruscode v10_06_00_13 -q e26:prof      # brings genie ($GENIE), gallery
setup larbatch                                  # project.py
htgettoken -a htvaultprod.fnal.gov -i icarus

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
3. Submits with `project.py --submit` and records the jobsub cluster IDs.
4. Polls `jobsub_q`. If every remaining job has been held for longer than
   `held_timeout_hours`, it runs `jobsub_rm` on them.
5. When the queue is empty, runs `--check` (`--checkana` for the CAF stage)
   and counts `files.list` / `filesana.list`.
6. Decides:
   - at least `min_good_fraction` good: done;
   - otherwise `--makeup`, up to `max_makeup` times;
   - then accept if at least `accept_fraction` is good, else the stage
     fails and the campaign stops.

**Operations.**
- External stages (SPINE, medulla) are `submit_cmd`/`check_cmd` shell
  commands in the config. `check_cmd` exits 0 for done, 1 for running and 2
  for failed.
- State is written atomically and a lock file prevents two instances, so
  `run` can be killed and restarted. `tick` also works from cron.
- `retry <stage>` resets a failed stage.
- `notify_cmd` (e.g. `mail`) fires on stage done/failed and when the token
  needs a refresh.

## 4. Pilot first (`gen_numjobs = 10`), then scale

Once the gen stage is done:

```bash
ls <pnfs_base>/gen/out/*/*.root > gen_files.list
root -l -b -q 'batch/ncdeltarad/dump_delta_decays.C("gen_files.list","delta_decays.csv")'
python batch/ncdeltarad/pilot_check.py delta_decays.csv
```

`pilot_check.py` tests the four assumptions the weighting rests on. It exits
non-zero if any fails:

1. **The patched table was used.** The radiative fraction per W bin must match
   the enhanced prediction, not the nominal one. This also catches a
   `TDatabasePDG` that was filled before GENIE loaded its table.
2. **CAF `w` equals the Δ mass GENIE decays with.** If it doesn't (e.g.
   binding), compute the weight from the Δ mass instead.
3. **The 30 MeV primary-photon tag reproduces the true decay channel.**
4. **`hitnuc` gives the Δ charge.**

**Sizing.** After the pilot has gone through CAF and medulla, measure the
number of stage1 category-0 events per job. Then
`gen_numjobs ≈ target / (that number)`. Each gen job uses one overlay input
file, so check the size of the overlay dataset.

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
  production XMLs for the same fcls.
- **The overlay SAM definition** that matches the run2 production.

## 6. Still needed in medulla

- **Truth branches** for the weight: true `w`, `hitnuc`, and the radiative tag
  (number of primary photons above 30 MeV).
- **A per-event weight in spineplot.** Today it only applies one POT scale per
  sample. `batch/gore_cut_opt.py` already takes explicit weights, so the cut
  re-tune can use the sample first.
- **The replacement logic:** drop CV NC `resnum == 0` events and add this
  sample with `delta_br.weight(...)`.
