#!/bin/bash

#######################################################################
# Grid payload: run a spineplot config on a systematics ROOT file and
# copy the resulting figures back to dCache. Companion to submit.sh
# (which runs the selection); this one runs ONLY spineplot.
#
# Why a separate payload: spineplot is pure Python (numpy, pandas,
# matplotlib, uproot, toml) — no cmake/make. But those packages are not
# in the sbnana CVMFS stack, so the payload builds a throwaway venv and
# pip-installs them at job start. It also needs lots of RAM: the
# cc_Xg1p_stage1 sideband (~18,912 MC events × 3 sys trees) peaks at
# ~11 GB in the develop-line Systematic.process path. Request >=16 GB
# at submission (see launch_spineplot.sh).
#
# Usage (invoked by jobsub via launch_spineplot.sh):
#   submit_spineplot.sh --input=PNFS_ROOT --output=PNFS_DIR \
#                       [--config=BASENAME] [--tag=REF] [--gituser=USER]
#
# Arguments:
#   --input=PNFS_ROOT  : dCache (/pnfs) path to the systematics ROOT
#                        file (e.g. output_gOre_1g1p_sys.root). MUST be
#                        on /pnfs — grid nodes cannot read /exp BlueArc.
#   --output=PNFS_DIR  : dCache directory to copy the figures into.
#   --config=BASENAME  : spineplot config basename under
#                        spineplot/configurations/analyses/icarus/
#                        (default: gOre_cc_Xg1p_stage1_datamc).
#   --tag=REF          : git ref to checkout (default:
#                        feature/hausnerh_gOre_1g1p).
#   --gituser=USER     : GitHub user to clone from (default: hausnerh).
#######################################################################

set -uo pipefail

# Initialize variables
INPUT=""
OUTPUT="/pnfs/icarus/scratch/users/hhausner/CCSidebandPlots"
CONFIG="gOre_cc_Xg1p_stage1_datamc"
TAG="feature/hausnerh_gOre_1g1p"
GITUSER="hausnerh"

usage() {
    grep '^#' "$0" | sed 's/^#//'
    exit 1
}

# Parse arguments (support both --flag=value and --flag value)
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input=*)   INPUT="${1#*=}";   shift ;;
    --input)     INPUT="$2";        shift 2 ;;
    --output=*)  OUTPUT="${1#*=}";  shift ;;
    --output)    OUTPUT="$2";       shift 2 ;;
    --config=*)  CONFIG="${1#*=}";  shift ;;
    --config)    CONFIG="$2";       shift 2 ;;
    --tag=*)     TAG="${1#*=}";     shift ;;
    --tag)       TAG="$2";          shift 2 ;;
    --gituser=*) GITUSER="${1#*=}"; shift ;;
    --gituser)   GITUSER="$2";      shift 2 ;;
    -h|--help)   usage ;;
    --)          shift; break ;;
    *) echo "Unknown option: $1" >&2; usage ;;
  esac
done

#######################################################################
# Check for required arguments
#######################################################################
missing_args=()
[[ -z "$INPUT" ]]  && missing_args+=("--input")
[[ -z "$OUTPUT" ]] && missing_args+=("--output")
if [[ ${#missing_args[@]} -gt 0 ]]; then
    echo "Error: Missing required argument(s): ${missing_args[*]}" >&2
    usage
fi

#######################################################################
# Initial Setup
#######################################################################

# IFDH options
export IFDH_CP_MAXRETRIES=2
export IFDH_WEB_TIMEOUT=100

# Setup CVMFS area + the sbnana stack. We do NOT need sbnana itself for
# spineplot, but `setup sbnana ... e26:prof` puts a modern python3
# (>=3.9) on PATH, which the throwaway venv is built from.
# The ups/cvmfs setup scripts reference unbound variables (e.g.
# `Options[@]` in setup_icarus.sh), which trip `set -u`; relax nounset
# for the environment setup only, then restore it for our own logic.
set +u
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup sbnana v10_01_02_01 -q e26:prof
ups active
set -u

echo "Using python: $(which python3)  ($(python3 --version 2>&1))"

#######################################################################
# Get the code
#######################################################################
git clone https://github.com/${GITUSER}/medulla.git
cd medulla
git checkout ${TAG}

#######################################################################
# Build a throwaway venv with the spineplot dependencies
#######################################################################
python3 -m venv .plotvenv
# shellcheck disable=SC1091
set +u                          # activate scripts may reference unbound vars
source .plotvenv/bin/activate
set -u
python3 -m pip install --upgrade pip
# Full spineplot dependency set. analysis.py eagerly imports every artist
# at module load, so scipy (spectra/efficiency/ternary) and scikit-learn
# (confusion/roc) are required even for a plain SpineSpectra1D config.
# uproot pulls in awkward automatically.
python3 -m pip install numpy scipy pandas matplotlib uproot toml scikit-learn
echo "venv packages:"
python3 -m pip list 2>/dev/null | grep -iE 'numpy|scipy|pandas|matplotlib|uproot|toml|scikit-learn'

#######################################################################
# Prestage the input systematics ROOT file from dCache
#######################################################################
mkdir -p data
echo "Copying input ROOT file: $INPUT"
ifdh cp "$INPUT" data/input_sys.root
ls -lrth data/

#######################################################################
# Run spineplot
#######################################################################
CFG_SRC="spineplot/configurations/analyses/icarus/${CONFIG}.toml"
if [[ ! -f "$CFG_SRC" ]]; then
    echo "Error: config not found at $CFG_SRC" >&2
    exit 1
fi

# Rewrite the config's `[output] path` to a local directory (the
# committed value points at /exp/... which does not exist on the node).
# Each gOre wrapper has exactly one top-level `path = ` line (the
# [output] block), so this sed is unambiguous.
mkdir -p plots
CFG_LOCAL="job_${CONFIG}.toml"
sed -E "s|^path = .*|path = 'plots'|" "$CFG_SRC" > "spineplot/$CFG_LOCAL"

# spineplot resolves its `[[this_includes]]` paths relative to CWD, so
# run from inside spineplot/ (matches run_gOre_1g1p_all.sh).
INPUT_ABS="$(pwd)/data/input_sys.root"
( cd spineplot && python3 spineplot.py --config "$CFG_LOCAL" --input "$INPUT_ABS" )
RC=$?

echo "spineplot exit code: $RC"
ls -lrth spineplot/plots/ 2>/dev/null

#######################################################################
# Copy the figures back to dCache
#######################################################################
shopt -s nullglob
made=(spineplot/plots/*.png spineplot/plots/*.pdf)
if [[ ${#made[@]} -eq 0 ]]; then
    echo "Error: spineplot produced no figures; not copying anything." >&2
    exit 1
fi

echo "Copying ${#made[@]} figure(s) to $OUTPUT"
ifdh mkdir_p "$OUTPUT" 2>/dev/null || true
for f in "${made[@]}"; do
    echo "  -> $(basename "$f")"
    ifdh cp "$f" "$OUTPUT/$(basename "$f")"
done

echo "Done."
exit "$RC"
