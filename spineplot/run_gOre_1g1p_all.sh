#!/usr/bin/env bash
#
# Run the gOre 1γ1p spineplot wrappers against the systematics file.
#
# Default set (6 configs, light): stage{1,2,3} × {MC-only, datamc}.
# These run comfortably on a standard gpvm.
#
# Heavy set (1 config, opt-in): gOre_cc_Xg1p_stage1_datamc. The cc
# sideband has ~18,912 MC events; with all three sys trees the
# develop-line Systematic.process path peaks at ~11 GB RSS and gets
# OOM-killed on a shared gpvm. It is SKIPPED unless INCLUDE_HEAVY=1,
# and should be run on a high-memory node or grid job (>=16 GB).
#
# Usage:
#   ./run_gOre_1g1p_all.sh [<input.root>]              # 6 light configs
#   INCLUDE_HEAVY=1 ./run_gOre_1g1p_all.sh [<input>]   # + cc sideband
#
# Defaults:
#   INPUT  = $INPUT env var, or build/output_gOre_1g1p_sys.root relative
#            to the repo root (the name systematics/toml/gOre_1g1p_sidebands.toml
#            writes to).
#   OUTBASE= /exp/icarus/app/users/hhausner/plots
#   LOGDIR = /tmp
#
# Behaviour:
#   - each output subdir is wiped (rm *.pdf *.png) the first time the
#     run targets it, so the dir reflects the current run only
#   - each wrapper toml is templated to $LOGDIR/<cfg>.toml with its
#     `[output] path` rewritten under $OUTBASE/<subdir>, then run
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
# Light set — runs on any node. Tuple: <config_basename>:<output_subdir>.
# Output subdir matches `[output] path` inside each TOML.
CONFIGS=(
    # NC 1γ1p — stage cuts baked into per-stage selection trees
    # `selected_1g1p_stage{1,2,3}` (driver: gOre_1g1p_sidebands.toml).
    "gOre_1g1p_stage1:gOre_1g1p/stage1_presel"
    "gOre_1g1p_stage2:gOre_1g1p/stage2_pi0_rej"
    "gOre_1g1p_stage3:gOre_1g1p/stage3_egam_sep"
    "gOre_1g1p_stage1_datamc:gOre_1g1p/stage1_presel_datamc"
    "gOre_1g1p_stage2_datamc:gOre_1g1p/stage2_pi0_rej_datamc"
    "gOre_1g1p_stage3_datamc:gOre_1g1p/stage3_egam_sep_datamc"
    # Δ-mass sideband wrappers retired — selection now bakes the
    # mass cut into `selected_1g1p_stage{2,3}` via `pi0_rejection`, so
    # the cut can't be inverted at the spineplot layer. Re-add once a
    # `selected_1g1p_sideband_stage<N>` tree is provided by selection.
)

# Heavy set — only run when INCLUDE_HEAVY=1 (needs a high-mem node /
# grid job, ~11 GB RSS). cc_Xg1p stage1 has ~18,912 MC events; full
# systematics OOM-kill a standard gpvm. See the header of
# gOre_cc_Xg1p_stage1_datamc.toml.
HEAVY_CONFIGS=(
    "gOre_cc_Xg1p_stage1_datamc:gOre_cc_Xg1p/stage1_presel_datamc"
    "gOre_cc_Xg1p_stage2_datamc:gOre_cc_Xg1p/stage2_pi0_rej_datamc"
    "gOre_cc_Xg1p_stage3_datamc:gOre_cc_Xg1p/stage3_egam_sep_datamc"
)

if [[ "${INCLUDE_HEAVY:-0}" == "1" ]]; then
    CONFIGS+=("${HEAVY_CONFIGS[@]}")
fi

# -------- Pre-flight --------
echo "Input  : $INPUT"
echo "Output : $OUTBASE/<stage>"
echo "Logs   : $LOGDIR/spineplot_<cfg>.log"
echo "Configs: ${#CONFIGS[@]}"
if [[ "${INCLUDE_HEAVY:-0}" != "1" ]]; then
    echo "Heavy  : cc_Xg1p_stage1_datamc SKIPPED (set INCLUDE_HEAVY=1 on a >=16 GB node to include it)"
fi
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
declare -A WIPED_DIRS
FAIL_COUNT=0

for entry in "${CONFIGS[@]}"; do
    cfg="${entry%%:*}"
    out_subdir="${entry##*:}"
    cfg_path="$CONFIG_DIR/$cfg.toml"
    log="$LOGDIR/spineplot_$cfg.log"
    out_dir="$OUTBASE/$out_subdir"
    templated_cfg="$LOGDIR/$cfg.toml"

    if [[ ! -f "$cfg_path" ]]; then
        echo ">>> SKIP $cfg (config not found at $cfg_path)"
        STATUS+=("MISSING")
        SURVIVAL+=("-")
        FAIL_COUNT=$((FAIL_COUNT + 1))
        continue
    fi

    # Wipe stale outputs the FIRST time this run targets each subdir,
    # then let subsequent configs targeting the same subdir append
    # (used by the memory-split cc_Xg1p_stage1_datamc_part{1,2} pair).
    if [[ -z "${WIPED_DIRS[$out_dir]:-}" ]]; then
        find "$out_dir" -maxdepth 1 -type f \( -name '*.pdf' -o -name '*.png' \) -delete 2>/dev/null || true
        WIPED_DIRS[$out_dir]=1
    fi

    # Template the wrapper: rewrite `[output] path = ...` to land under
    # $OUTBASE/<out_subdir>. This is what makes $OUTBASE actually
    # override the hard-coded toml output path. sed-style escape for
    # the replacement so a path with slashes is safe.
    out_dir_escaped="$(printf '%s' "$out_dir" | sed 's:[\\/&]:\\&:g')"
    sed -E "s|^path = .*|path = '${out_dir_escaped}'|" "$cfg_path" > "$templated_cfg"

    echo ">>> Running $cfg  ->  $out_dir"
    # cd into spineplot/ so the configs' relative `[[this_includes]]`
    # paths (`configurations/common/styles.toml`) resolve.
    if (cd "$SCRIPT_DIR" && python "$MAIN_PY" --config "$templated_cfg" --input "$INPUT") > "$log" 2>&1; then
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
