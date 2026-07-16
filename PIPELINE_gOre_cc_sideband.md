# gOre CC Xγ+1p sideband — production pipeline (with `cc_sideband_category`)

End-to-end recipe for producing the **CC Xγ+1p sideband** data/MC plots
broken down by the new CC-oriented truth categorization
`vars::gOre::cc_sideband_category` (saved branch **`true_cc_sideband_category`**).

> **This pipeline is written for high-statistics grid production as the
> default.** Every heavy stage runs on the grid; the gpvm is used only to
> submit, `hadd`, and run the one lightweight systematics pass. A fully
> local, low-statistics run of the *same* chain is kept at the end
> ([Appendix A](#appendix-a--local-machinery-smoke-test)) purely as a
> machinery smoke-test to confirm the wiring before spending grid cycles.

The default `[[category]]` scheme lumps every CC interaction into one
`Other CC` bin, useless in this CC-dominated region. The new variable
mirrors the NC `[[category]]` breakdown for CC, splitting the Δ→Nγ signal
by **photon** multiplicity (1γ vs Nγ) instead of proton multiplicity, and
— because the region's reco cut (`cc_Xg_topology`) already requires a muon
— gating the topology bins on νμ so the `CC νₑ` bin is a **complete**
electron-neutrino sample.

| idx | category | definition (truth) |
|----|----|----|
| 0 | CC Δ→Nγ (1γ)          | νμ CC, Δ res (resnum 0), exactly 1 primary γ, no pions |
| 1 | CC Δ→Nγ (Nγ)          | νμ CC, Δ res, ≥2 primary γ, no pions |
| 2 | Other CC Xγ Post-FSI  | νμ CC, non-res, photon topology (no pions) |
| 3 | CC π⁰ (Δ Res)         | νμ CC, Δ res with a π⁰ |
| 4 | CC π⁰ (Other)         | νμ CC, non-res with a π⁰ |
| 5 | CC π±                 | νμ CC with a charged pion |
| 6 | CC νₑ                 | every νₑ CC (\|pdg\| == 12) |
| 7 | Other CC              | remaining CC (νμ with no tagged topology, ντ) |
| 8 | NC (all)              | any NC interaction (lumped) |
| 9 | Cosmic / Non-ν        | unmatched / not a neutrino |

See [Appendix B](#appendix-b--category-reference--tuning) to change the
definitions.

---

## Data flow (high-stats / grid)

```
 medulla.py --create-project           grid: N jobs, one sample-batch each
   -> project.db + systematics.toml       each self-builds this branch,
   medulla.py --launch-jobs               runs selection -> output_jobidNNNN.root
                                          copied to  $PROJECT/output/
        |
        | hadd $PROJECT/output/output_jobid*.root
        v
   output_gOre_1g1p.root   (full-stat selection, all samples)
        |
        | run_systematics  (once, on a gpvm)
        v
   output_gOre_1g1p_sys.root
        |
        | ifdh cp -> /pnfs ; launch_spineplot.sh  (grid, 16 GB)
        v
   figures on /pnfs  ->  ifdh cp back
```

Each grid job *also* runs a per-job systematics pass into
`output_systematics_jobid*.root`, but the gOre flow does **not** use those
(they come from an auto-generated template). Systematics is run once, from
the gOre-specific config, on the hadd-ed selection — matching the
`[input] path = 'output_gOre_1g1p.root'` contract in
`systematics/toml/gOre_1g1p_sidebands.toml`.

---

## Prerequisites (once per session)

```bash
# ICARUS environment
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup sbnana v10_01_02_01 -q e26:prof
setup cmake  v3_27_4

# grid token (needed to submit and for ifdh)
htgettoken -a htvaultprod.fnal.gov -i icarus
```

You do **not** need a local build for grid running — each job builds this
branch itself. A build is only needed for the smoke-test in Appendix A.

> **Fork + branch awareness — get this right or the run is wasted.** The
> grid clones `github.com/<gituser>/medulla`, checks out `<tag>`, and
> **builds `medulla` from that ref**. So `<tag>` must be the branch that
> actually contains `cc_sideband_category`; otherwise the built binary has
> no `true_cc_sideband_category` var and the selection errors on the branch
> your `project.db` requests. That code currently lives **only** on
> `worktree-gOre-cc-sideband-category` — `feature/hausnerh_gOre_1g1p` does
> **not** have it yet. Every `medulla.py` call below therefore uses:
>
> ```bash
> TAG=worktree-gOre-cc-sideband-category   # branch that has cc_sideband_category
> GITUSER=hausnerh                         # your fork owner
> ```
>
> Once you merge this into `feature/hausnerh_gOre_1g1p`, set
> `TAG=feature/hausnerh_gOre_1g1p` instead. `--gituser` is also what the
> pre-flight validates against; omit it and the default `justinjmueller` is
> checked and you get *"Tag '…' does not exist"*.

---

## Stage 1 — grid selection (high statistics)

`--create-project` and `--launch-jobs` are **two separate invocations**
(the launch guard is evaluated against the project state *before*
creation). `--batch-size` is files-per-job — tune it so each job fits in
`--memory`.

> **Resource flags belong on the launch call.** `--memory / --disk /
> --lifetime` are consumed by the jobsub submission, so they only take
> effect on `--launch-jobs`. Passing them to `--create-project` silently
> does nothing and your jobs run with the defaults (1800 MB / 1h / 25 GB).

```bash
PROJ=/pnfs/icarus/scratch/users/$USER/gOre_cc_sideband
# TAG / GITUSER set in Prerequisites above.

# 1a. build the project (splits samples into jobs, writes project.db).
#     --tag/--gituser here only satisfy the pre-flight ref check; project
#     creation reads the LOCAL toml, so run it from a checkout that has
#     the cc_sideband_category branch line.
python3 batch/medulla.py --experiment icarus \
    --project-dir $PROJ \
    --create-project --toml selection/toml/gOre_1g1p_sidebands.toml \
    --batch-size 20 \
    --tag $TAG --gituser $GITUSER

# 1b. launch all pending jobs — resources go HERE:
python3 batch/medulla.py --experiment icarus \
    --project-dir $PROJ \
    --tag $TAG --gituser $GITUSER \
    --launch-jobs --memory 4000 --disk 20 --lifetime 3h
```

Confirm the branch actually made it into a stored job config *before* a big
launch — this catches a `project.db` built from a checkout that lacked the
branch line (the grid rewrites `[general] output = "output"`, so job files
are `output_jobidNNNN.root`).

`project.db` lives on `/pnfs`, and SQLite cannot open a DB in place on
dCache (`Error: stepping, disk I/O error (10)`), so copy it to a local
POSIX path first, then query the copy. Use `ifdh cp` — a plain `cp` off
`/pnfs` can fail with `Operation not permitted` when the NFS credential
has lapsed:

```bash
ifdh cp $PROJ/project.db /tmp/project.db
sqlite3 /tmp/project.db \
  "SELECT cfg FROM configuration LIMIT 1;" | grep cc_sideband_category
# for the stage2/3 rerun, grep selected_cc_Xg1p_stage3 instead
```

> If either the `ifdh cp` above or the `--launch-jobs` step fails on
> `/pnfs` access (`Operation not permitted`), refresh credentials —
> `kinit` (for the NFS `/pnfs` mount) and
> `htgettoken -a htvaultprod.fnal.gov -i icarus` (bearer token) — then
> retry. `--launch-jobs` and the status refresh copy `project.db` off
> `/pnfs` internally, so they need `/pnfs` access healthy too.

---

## Stage 2 — monitor and gather

Re-running `medulla.py` on an existing project refreshes job status from
the output directory (an output ≥1 KB marks a job `completed`):

```bash
python3 batch/medulla.py --experiment icarus --project-dir $PROJ \
    --tag $TAG --gituser $GITUSER
```

When the jobs are in, hadd the per-job **selection** outputs into the
single file the systematics stage expects:

```bash
hadd -f build/output_gOre_1g1p.root $PROJ/output/output_jobid*.root
# verify the CC tree + new branch survived the hadd:
root -l -b -q -e 'TFile f("build/output_gOre_1g1p.root"); \
  f.Get("events/cv/selected_cc_Xg1p_stage1")->Print();' 2>&1 | grep cc_sideband_category
```

---

## Stage 3 — systematics (once, on a gpvm)

`run_systematics` is the C++ binary and is not the memory hog (that is the
Python spineplot pass in Stage 4), so run it once on the full hadded
sample. It reads `output_gOre_1g1p.root` from the working directory.

```bash
cd build   # where output_gOre_1g1p.root lives; [input] path is relative
./systematics/run_systematics ../systematics/toml/gOre_1g1p_sidebands.toml
# -> output_gOre_1g1p_sys.root
cd ..
```

The `true_cc_sideband_category` column is copied through automatically —
the systematics config references trees by name/origin and does not
re-declare the branch list.

*(If Stage 3 itself grows too heavy at full stats, it can instead be run
per-batch on the grid by passing `--systematic
systematics/toml/gOre_1g1p_sidebands.toml` in Stage 1a and hadding the
`output_systematics_jobid*.root` files — but the single-pass route above
is the gOre default.)*

---

## Stage 4 — sideband plots on the grid (heavy)

The CC config `gOre_cc_Xg1p_stage1_datamc.toml` already points at
`category_branch = 'true_cc_sideband_category'` with the 10 CC labels +
Data overlay. At full stats the spineplot systematics pass peaks ~11 GB
RSS and is OOM-killed on a shared gpvm, so it runs on the grid via the
branch's fork-aware payload (requests 16 GB, defaults already target this
config / branch / fork).

Grid nodes cannot read `/exp`, so stage the systematics ROOT on dCache
first:

```bash
OUT=/pnfs/icarus/scratch/users/$USER/CCSidebandPlots
ifdh cp build/output_gOre_1g1p_sys.root $OUT/output_gOre_1g1p_sys.root

for cfg in gOre_cc_Xg1p_stage1_datamc \
           gOre_cc_Xg1p_stage2_datamc \
           gOre_cc_Xg1p_stage3_datamc ; do
  ./batch/launch_spineplot.sh \
      --input=$OUT/output_gOre_1g1p_sys.root \
      --output=$OUT/$cfg \
      --config=$cfg \
      --tag=$TAG --gituser=$GITUSER
done
# Three stages of the CC sideband, mirroring the actual selection cut chain
# but keeping the muon (cc_Xg_topology) for sideband perspective:
#   stage1 = preselection
#   stage2 = + pi0_rejection
#   stage3 = + pi0_rejection + egamma_separation
# --config MUST be overridden per stage; --tag MUST be $TAG (the built-in
# default feature/hausnerh_gOre_1g1p lacks these configs and the new
# true_cc_sideband_category branch).
```

---

## Stage 5 — retrieve figures

```bash
ifdh cp -r /pnfs/icarus/scratch/users/$USER/CCSidebandPlots \
           /exp/icarus/app/users/$USER/plots/gOre_cc_Xg1p
```

Variables plotted (each stacked by `cc_sideband_category`):
leading-shower KE / start dE/dx / photon-softmax / primary-softmax /
directional & axial spread, reco Δ-mass (and |Δm−1232|), leading-proton
primary-softmax and KE.

---

## Appendix A — local machinery smoke-test

Run the whole chain on one node against a **small** input to prove the
wiring (branch present, systematics copies it through, plots render)
before committing grid resources. This is a correctness check, not a
physics result — the stats will be poor.

```bash
# 0. build once
mkdir -p build && cd build && cmake .. && make -j4 && cd ..

# 1. shrink the input: point ONE sample at a single flat.root and disable
#    the rest, e.g. copy the toml and edit [[sample]] paths/disable:
cp selection/toml/gOre_1g1p_sidebands.toml /tmp/gOre_smoke.toml
#   (set every [[sample]] disable=true except cv; point cv 'path' at one file)

# 2. selection -> writes gOre_1g1p_sidebands.root
./build/selection/medulla /tmp/gOre_smoke.toml
mv gOre_1g1p_sidebands.root build/output_gOre_1g1p.root

# 3. systematics -> output_gOre_1g1p_sys.root
cd build && ./systematics/run_systematics ../systematics/toml/gOre_1g1p_sidebands.toml && cd ..

# 4. plots locally (the all-in-one driver; INCLUDE_HEAVY pulls in the CC
#    config, fine at smoke-test stats):
INCLUDE_HEAVY=1 OUTBASE=/tmp/gOre_plots \
  spineplot/run_gOre_1g1p_all.sh build/output_gOre_1g1p_sys.root
#   -> /tmp/gOre_plots/gOre_cc_Xg1p/stage1_presel_datamc/
```

To run just the one CC config directly:

```bash
python3 spineplot/spineplot.py \
  spineplot/configurations/analyses/icarus/gOre_cc_Xg1p_stage1_datamc.toml
```

---

## Appendix B — category reference / tuning

Definitions live in one function,
`selection/include/gOre/vars_gOre.h::cc_sideband_category`. Photon counting
is over post-FSI primaries above the gOre threshold, so π⁰ decay photons
(the π⁰ is the primary) are not counted. Precedence is photon topology →
π⁰ → π± → νₑ. To change binning or ordering, reorder the `return` branches
and update the matching `category_labels` / `category_colors` /
`category_assignment` in
`spineplot/configurations/analyses/icarus/gOre_cc_Xg1p_stage1_datamc.toml`.

**After any edit to the variable you must rebuild** (`make -j4`) and
re-run Stages 1–4; on the grid this happens automatically because each job
rebuilds the branch you pass via `--tag`.

---

## Notes / gotchas

- **`--tag` must contain the var, and match the code you built the project
  from** — the grid builds `medulla` from `--tag`. Point it at
  `worktree-gOre-cc-sideband-category` (not `feature/hausnerh_gOre_1g1p`,
  which lacks `cc_sideband_category`) until the two are merged.
- **`--gituser hausnerh` on every `medulla.py` call** — the pre-flight tag
  check and the grid clone both use it; the default `justinjmueller` does
  not have this branch.
- **Resources on `--launch-jobs`, not `--create-project`** — `--memory /
  --disk / --lifetime` are ignored by project creation; on the launch call
  they map to the jobsub request.
- **Two calls for create + launch** — `--create-project` then
  `--launch-jobs`; combining them raises `FileNotFoundError`.
- **`/pnfs`, not `/exp`, for grid I/O** — job inputs/outputs and the
  spineplot payload's input must be on dCache.
- **Data (onbeam)** has no truth; the config forces every onbeam event to
  category 10 via `[samples.onbeam.precompute]`. MC events with no matched
  truth neutrino go to bin 9 (`Cosmic / Non-ν`) via `fillna = 9`.
- **NC stages unaffected** — the `selected_1g1p_stage{1,2,3}` plots still
  use the default `true_category`; only the CC sideband uses the new branch.
