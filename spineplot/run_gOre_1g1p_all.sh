#!/usr/bin/env bash
#
# Run every gOre 1γ1p spineplot wrapper against the systematics file.
# Covers 7 configs: stage{1,2,3} × {MC-only, datamc} + cc_Xg1p stage1 datamc.
#
# Usage:
#   ./run_gOre_1g1p_all.sh [<input.root>]
#
# Defaults:
#   INPUT  = $INPUT env var, or build/output_gOre_1g1p_sys.root relative
#            to the repo root (the name systematics/toml/gOre_1g1p_sidebands.toml
#            writes to).
#   OUTBASE= /exp/icarus/app/users/hhausner/plots
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
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
INPUT="${1:-${INPUT:-$REPO_ROOT/build/output_gOre_1g1p_sys.root}}"
OUTBASE="${OUTBASE:-/exp/icarus/app/users/hhausner/plots}"
LOGDIR="${LOGDIR:-/tmp}"

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input ROOT file does not exist: $INPUT" >&2
    exit 2
fi

# Locate configurations relative to SCRIPT_DIR (set above) so the
# script works from anywhere.
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
    # NC 1γ1p — stage cuts baked into per-stage selection trees
    # `selected_1g1p_stage{1,2,3}` (driver: gOre_1g1p_sidebands.toml).
    "gOre_1g1p_stage1:gOre_1g1p/stage1_presel"
    "gOre_1g1p_stage2:gOre_1g1p/stage2_pi0_rej"
    "gOre_1g1p_stage3:gOre_1g1p/stage3_egam_sep"
    "gOre_1g1p_stage1_datamc:gOre_1g1p/stage1_presel_datamc"
    "gOre_1g1p_stage2_datamc:gOre_1g1p/stage2_pi0_rej_datamc"
    "gOre_1g1p_stage3_datamc:gOre_1g1p/stage3_egam_sep_datamc"
    # CC + Xγ + 1p sideband — only `_stage1` tree exists in the v3
    # systematics output. Stage 2 and 3 wrappers retired until matching
    # selection trees are added to `selection/toml/gOre_1g1p_sidebands.toml`.
    "gOre_cc_Xg1p_stage1_datamc:gOre_cc_Xg1p/stage1_presel_datamc"
    # Δ-mass sideband wrappers retired — selection now bakes the
    # mass cut into `selected_1g1p_stage{2,3}` via `pi0_rejection`, so
    # the cut can't be inverted at the spineplot layer. Re-add once a
    # `selected_1g1p_sideband_stage<N>` tree is provided by selection.
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
    # cd into spineplot/ so the configs' relative `[[this_includes]]`
    # paths (`configurations/common/styles.toml`) resolve.
    if (cd "$SCRIPT_DIR" && python "$MAIN_PY" --config "$cfg_path" --input "$INPUT") > "$log" 2>&1; then
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
