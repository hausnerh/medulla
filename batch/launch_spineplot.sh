#!/bin/bash

#######################################################################
# Launch a single grid job that runs spineplot on a systematics ROOT
# file and copies the figures to dCache. Standalone analog of
# utilities.launch_jobsub, but for plotting: there is exactly one job,
# so no project.db / jobs table is needed.
#
# The cc_Xg1p_stage1 sideband needs ~11 GB RAM (see
# submit_spineplot.sh), hence --memory=16000MB and DEDICATED/
# OPPORTUNISTIC only (OFFSITE dropped: 16 GB is hard to match offsite
# and the runtime pip install wants reliable egress).
#
# Prerequisites:
#   1. A valid token:   htgettoken -a htvaultprod.fnal.gov -i icarus
#   2. The input systematics ROOT file staged on dCache (/pnfs), NOT
#      /exp. Grid nodes cannot read BlueArc /exp. Stage it once with:
#        ifdh cp /exp/icarus/data/users/hhausner/SPINE/medulla/build/output_gOre_1g1p_sys.root \
#                /pnfs/icarus/scratch/users/hhausner/CCSidebandPlots/output_gOre_1g1p_sys.root
#
# Usage:
#   ./launch_spineplot.sh [--input=PNFS_ROOT] [--output=PNFS_DIR]
#                         [--config=BASENAME] [--tag=REF] [--gituser=USER]
#
# Defaults:
#   --input   = /pnfs/icarus/scratch/users/hhausner/CCSidebandPlots/output_gOre_1g1p_sys.root
#   --output  = /pnfs/icarus/scratch/users/hhausner/CCSidebandPlots
#   --config  = gOre_cc_Xg1p_stage1_datamc
#   --tag     = feature/hausnerh_gOre_1g1p
#   --gituser = hausnerh
#######################################################################

set -uo pipefail

EXP="icarus"
INPUT="/pnfs/icarus/scratch/users/hhausner/CCSidebandPlots/output_gOre_1g1p_sys.root"
OUTPUT="/pnfs/icarus/scratch/users/hhausner/CCSidebandPlots"
CONFIG="gOre_cc_Xg1p_stage1_datamc"
TAG="feature/hausnerh_gOre_1g1p"
GITUSER="hausnerh"

usage() {
    grep '^#' "$0" | sed 's/^#//'
    exit 1
}

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
    *) echo "Unknown option: $1" >&2; usage ;;
  esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAYLOAD="$SCRIPT_DIR/submit_spineplot.sh"
if [[ ! -f "$PAYLOAD" ]]; then
    echo "Error: payload not found at $PAYLOAD" >&2
    exit 1
fi

cmd=(
    jobsub_submit
    -G "$EXP"
    -N 1
    --memory=16000MB
    --disk=25GB
    --expected-lifetime=2h
    --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC
    "--append_condor_requirements='(TARGET.HAS_Singularity==true)'"
    --singularity-image=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest
    "file://$PAYLOAD"
    --
    "--input=$INPUT"
    "--output=$OUTPUT"
    "--config=$CONFIG"
    "--tag=$TAG"
    "--gituser=$GITUSER"
)

echo "[INFO] -- Launching spineplot grid job:"
printf '  %s\n' "${cmd[*]}"
read -r -p "Confirm job launch? [Y/N] " resp
if [[ "${resp,,}" != "y" ]]; then
    echo "[INFO] -- User aborted job launch."
    exit 0
fi

"${cmd[@]}"
