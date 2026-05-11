#!/usr/bin/env bash
#
# Run every gOre 1γ1p spineplot wrapper against the v2 systematics file.
# Covers 9 configs: stage{1,2,3} × {MC-only, datamc, datamc_sideband}.
#
# Usage:
#   ./run_gOre_1g1p_all.sh [<input.root>]
#
# Defaults:
#   INPUT  = $INPUT env var, or the v2 hadded systematics file at
#            /exp/icarus/data/users/hhausner/SPINE/medulla/build/output_gOre_systematics_v2.root
#   OUTBASE= /exp/icarus/app/users/hhausner/plots/gOre_1g1p
#   LOGDIR = /tmp
#
# Behaviour:
#   - mkdir -p every output subdir before running
#   - one config per call, sequential (matplotlib is not thread-safe)
#   - per-config log under $LOGDIR/spineplot_<cfg>.log
#   - failure on one config does NOT abort the rest
#   - summary table at end: status + onbeam survival grep
#
# Exit code: number of failed configs (0 = all green).

set -uo pipefail

# -------- Inputs --------
INPUT="${1:-${INPUT:-/exp/icarus/data/users/hhausner/SPINE/medulla/build/output_gOre_systematics_v2.root}}"
OUTBASE="${OUTBASE:-/exp/icarus/app/users/hhausner/plots}"
LOGDIR="${LOGDIR:-/tmp}"

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input ROOT file does not exist: $INPUT" >&2
    exit 2
fi

# Resolve script dir → spineplot/, then locate configurations relative
# to it so the script works from anywhere.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_DIR="$SCRIPT_DIR/configurations/analyses/icarus"
MAIN_PY="$SCRIPT_DIR/spineplot.py"

if [[ ! -f "$MAIN_PY" ]]; then
    echo "ERROR: spineplot.py not found at $MAIN_PY" >&2
    exit 2
fi

# -------- Config list --------
# Tuple: <config_basename>:<output_subdir_relative_to_OUTBASE>
# Output subdir matches `[output] path` inside each TOML so mkdir -p
# pre-creates the right tree.
CONFIGS=(
    "gOre_1g1p_stage1:gOre_1g1p/stage1_presel"
    "gOre_1g1p_stage2:gOre_1g1p/stage2_pi0_rej"
    "gOre_1g1p_stage3:gOre_1g1p/stage3_egam_sep"
    "gOre_1g1p_stage1_datamc:gOre_1g1p/stage1_presel_datamc"
    "gOre_1g1p_stage2_datamc:gOre_1g1p/stage2_pi0_rej_datamc"
    "gOre_1g1p_stage3_datamc:gOre_1g1p/stage3_egam_sep_datamc"
    "gOre_1g1p_stage1_datamc_sideband:gOre_1g1p/stage1_presel_datamc_sideband"
    "gOre_1g1p_stage2_datamc_sideband:gOre_1g1p/stage2_pi0_rej_datamc_sideband"
    "gOre_1g1p_stage3_datamc_sideband:gOre_1g1p/stage3_egam_sep_datamc_sideband"
    # CC + Xγ + 1p sideband — added 2026-05-08, depends on the
    # `selected_cc_Xg1p` tree (requires selection + systematics rerun).
    "gOre_cc_Xg1p_stage1_datamc:gOre_cc_Xg1p/stage1_presel_datamc"
    "gOre_cc_Xg1p_stage2_datamc:gOre_cc_Xg1p/stage2_kinematic_datamc"
    "gOre_cc_Xg1p_stage3_datamc:gOre_cc_Xg1p/stage3_shower_id_datamc"
)

# -------- Pre-flight --------
echo "Input  : $INPUT"
echo "Output : $OUTBASE/<stage>"
echo "Logs   : $LOGDIR/spineplot_<cfg>.log"
echo "Configs: ${#CONFIGS[@]}"
echo

for entry in "${CONFIGS[@]}"; do
    out_subdir="${entry##*:}"
    mkdir -p "$OUTBASE/$out_subdir" || {
        echo "ERROR: cannot create $OUTBASE/$out_subdir" >&2
        exit 2
    }
done

# -------- Run --------
declare -a STATUS
declare -a SURVIVAL
FAIL_COUNT=0

for entry in "${CONFIGS[@]}"; do
    cfg="${entry%%:*}"
    cfg_path="$CONFIG_DIR/$cfg.toml"
    log="$LOGDIR/spineplot_$cfg.log"

    if [[ ! -f "$cfg_path" ]]; then
        echo ">>> SKIP $cfg (config not found at $cfg_path)"
        STATUS+=("MISSING")
        SURVIVAL+=("-")
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi

    echo ">>> Running $cfg"
    if python "$MAIN_PY" --config "$cfg_path" --input "$INPUT" > "$log" 2>&1; then
        STATUS+=("OK")
        # Pull final onbeam survival count if the wrapper has data overlay.
        line="$(grep "Sample 'onbeam' presel:" "$log" | tail -1 || true)"
        if [[ -n "$line" ]]; then
            SURVIVAL+=("${line#*presel: }")
        else
            SURVIVAL+=("MC-only")
        fi
    else
        STATUS+=("FAIL")
        SURVIVAL+=("see log")
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
done

# -------- Summary --------
echo
echo "================ summary ================"
printf "%-45s %-8s %s\n" "config" "status" "onbeam"
printf "%-45s %-8s %s\n" "---------------------------------------------" "--------" "------"
for i in "${!CONFIGS[@]}"; do
    cfg="${CONFIGS[$i]%%:*}"
    printf "%-45s %-8s %s\n" "$cfg" "${STATUS[$i]}" "${SURVIVAL[$i]}"
done
echo "========================================="
echo "Failed: $FAIL_COUNT / ${#CONFIGS[@]}"

exit "$FAIL_COUNT"
