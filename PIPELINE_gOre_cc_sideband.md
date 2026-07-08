# gOre CC Xγ+1p sideband — production pipeline (with `cc_sideband_category`)

End-to-end recipe for producing the **CC Xγ+1p sideband** data/MC plots
broken down by the new CC-oriented truth categorization
`vars::gOre::cc_sideband_category` (saved branch **`true_cc_sideband_category`**).

The default `[[category]]` scheme lumps every CC interaction into one
`Other CC` bin, which is useless in this CC-dominated region. The new
variable mirrors the NC `[[category]]` breakdown for CC, splitting the
Δ→Nγ signal by **photon** multiplicity (1γ vs Nγ) instead of proton
multiplicity, and — because the region's reco cut (`cc_Xg_topology`)
already requires a muon — gating the topology bins on νμ so the `CC νₑ`
bin is a **complete** electron-neutrino sample.

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

Photon counting is over post-FSI primaries above the gOre threshold, so
π⁰ decay photons (the π⁰ is the primary) are not counted. Precedence is
photon topology → π⁰ → π± → νₑ; to change it, reorder the `return`
branches in `selection/include/gOre/vars_gOre.h::cc_sideband_category`
and update the labels in the spineplot config.

---

## Step 0 — environment + build (REQUIRED after this change)

The new category is a C++ variable, so the selection binary **must be
rebuilt** before the branch appears in the output.

```bash
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup sbnana v10_01_02_01 -q e26:prof
setup cmake  v3_27_4

cd <medulla>            # repo root
mkdir -p build && cd build
cmake .. && make -j4
```

Build products (README §Build):
- `build/selection/medulla`         — TOML-configured selection
- `build/systematics/run_systematics` — systematics / reweighting

Sanity-check the new var is registered before a long run:

```bash
./build/selection/medulla selection/toml/gOre_1g1p_sidebands.toml 2>&1 | head
# then confirm the branch exists on the CC tree:
root -l -b -q -e 'TFile f("gOre_1g1p_sidebands.root"); \
  f.Get("events/cv/selected_cc_Xg1p_stage1")->Print();' 2>&1 | grep cc_sideband_category
```

---

## Step 1 — selection

Runs the 3 NC stages + the CC sideband tree over the CAF samples. The
CC tree `selected_cc_Xg1p_stage1` now carries `true_cc_sideband_category`.

```bash
# from the repo root, XRootD token required for /pnfs input
./build/selection/medulla selection/toml/gOre_1g1p_sidebands.toml
# -> writes ./gOre_1g1p_sidebands.root  ([general] output = "gOre_1g1p_sidebands")
```

Grid (recommended for full statistics — see Step 1b) produces one file
per job; `hadd` them into the single file the systematics stage expects:

```bash
hadd -f build/output_gOre_1g1p.root <grid-outputs>/*.root
```

> The systematics config reads `output_gOre_1g1p.root`
> (`systematics/toml/gOre_1g1p_sidebands.toml [input] path`). A local
> single-shot run writes `gOre_1g1p_sidebands.root`; rename or `hadd` it
> to `output_gOre_1g1p.root` before Step 2.

### Step 1b — selection on the grid (optional, high stats)

`batch/medulla.py` can run this branch on the grid. `--memory / --disk /
--lifetime` size the jobs; `--tag` picks the git ref; and `--gituser`
picks the GitHub owner to clone from **and** to validate `--tag` against.

> **`--gituser` is required for a fork branch.** The grid job clones from
> `github.com/<gituser>/medulla` and the pre-flight check validates the
> tag against the same repo. It defaults to `justinjmueller`, so a branch
> that only lives on your fork fails with *"Tag '…' does not exist"* unless
> you pass `--gituser hausnerh`. (`--project-dir` is also required on every
> invocation.)

```bash
# create the project (note --project-dir AND --gituser):
python3 batch/medulla.py --experiment icarus \
    --project-dir /pnfs/icarus/scratch/users/$USER/gOre_1g1p_sidebands \
    --create-project --toml selection/toml/gOre_1g1p_sidebands.toml \
    --batch-size <N> --tag feature/hausnerh_gOre_1g1p --gituser hausnerh \
    --memory 4000 --disk 20 --lifetime 2h

# then launch:
python3 batch/medulla.py --experiment icarus \
    --project-dir /pnfs/icarus/scratch/users/$USER/gOre_1g1p_sidebands \
    --tag feature/hausnerh_gOre_1g1p --gituser hausnerh --launch-jobs
```

---

## Step 2 — systematics (multisim + detsys + variation weights)

```bash
./build/systematics/run_systematics systematics/toml/gOre_1g1p_sidebands.toml
# reads  output_gOre_1g1p.root
# writes output_gOre_1g1p_sys.root
```

The `true_cc_sideband_category` column is copied through automatically —
the systematics config references trees by name/origin and does not
re-declare the branch list, so no edit is needed there.

---

## Step 3 — sideband plots (spineplot)

The CC sideband config `gOre_cc_Xg1p_stage1_datamc.toml` already points at
`category_branch = 'true_cc_sideband_category'` with the 10 CC labels +
Data overlay. It is the **heavy** config (~18.9k MC events; full
systematics peak ~11 GB RSS) and is OOM-killed on a shared gpvm, so the
all-in-one driver skips it unless `INCLUDE_HEAVY=1` and it is run on a
≥16 GB node (or the grid).

**Option A — high-memory interactive node** (≥16 GB):

```bash
INCLUDE_HEAVY=1 \
OUTBASE=/exp/icarus/app/users/$USER/plots \
  spineplot/run_gOre_1g1p_all.sh build/output_gOre_1g1p_sys.root

# -> PDFs/PNGs under $OUTBASE/gOre_cc_Xg1p/stage1_presel_datamc/
```

**Option B — grid (recommended for the heavy config).** The branch ships a
dedicated, fork-aware spineplot payload that requests 16 GB. Stage the
systematics ROOT on dCache first (grid nodes cannot read `/exp`):

```bash
htgettoken -a htvaultprod.fnal.gov -i icarus
ifdh cp build/output_gOre_1g1p_sys.root \
        /pnfs/icarus/scratch/users/$USER/CCSidebandPlots/output_gOre_1g1p_sys.root

# defaults already target this config / branch / fork (hausnerh):
./batch/launch_spineplot.sh \
    --input=/pnfs/icarus/scratch/users/$USER/CCSidebandPlots/output_gOre_1g1p_sys.root \
    --output=/pnfs/icarus/scratch/users/$USER/CCSidebandPlots
```

To run just that one config directly on a fat node (bypassing the driver):

```bash
python3 spineplot/spineplot.py \
  spineplot/configurations/analyses/icarus/gOre_cc_Xg1p_stage1_datamc.toml
```

Output variables plotted (each stacked by `cc_sideband_category`):
leading-shower KE / start dE/dx / photon-softmax / primary-softmax /
directional & axial spread, reco Δ-mass (and |Δm−1232|), leading-proton
primary-softmax and KE.

---

## Notes / gotchas

- **Rebuild is mandatory** — skipping Step 0 leaves the old output with
  no `true_cc_sideband_category` branch; spineplot then errors on the
  missing `category_branch` (or, for onbeam, only the precompute-set
  value survives).
- **Data (onbeam)** has no truth; the config forces every onbeam event to
  category 10 via `[samples.onbeam.precompute]`. MC events with no matched
  truth neutrino go to bin 9 (`Cosmic / Non-ν`) via `fillna = 9`.
- **Consistency with the NC config** — the NC stage plots still use the
  default `true_category`; only the CC sideband uses the new branch.
- **Tuning categories** — definitions live in one function
  (`cc_sideband_category`); after editing it you must rebuild (Step 0)
  and re-run Steps 1–3.
