#!/bin/bash
#######################################################################
# One-shot, unattended CC Xγ+1p sideband production.
#
#   selection (grid) -> monitor -> hadd -> systematics -> stage
#   -> plots stage1/2/3 (grid) -> monitor -> fetch figures
#
# Monitors the grid itself: polls the queue, RELEASES held jobs with
# bumped memory/disk, and on any failure FETCHES the job logs into a
# local debug dir. Also monitors CREDENTIALS: on lapsed kerberos or
# bearer token it logs the exact fix commands and busy-waits until the
# user renews them (works foreground or under nohup).
#
# Designed to be launched and left alone:
#   nohup ./batch/run_cc_sideband.sh > run_cc_sideband.nohup 2>&1 &
#
# Resume mode: skip Stage 1 (selection submit) and continue from an
# existing project's outputs, e.g. if a previous run had a handful of
# stragglers that never completed:
#   ./batch/run_cc_sideband.sh --resume /pnfs/icarus/scratch/users/$USER/gOre_cc_sideband_20260716_134324
#   ./batch/run_cc_sideband.sh --resume /pnfs/.../gOre_cc_sideband_20260716_134324/output/
#
# Run from the medulla repo root, on a gpvm, with a VALID token +
# kerberos ticket (kinit ; htgettoken -a htvaultprod.fnal.gov -i icarus).
#
# NOTE: the jobsub_q / condor_qedit / jobsub_fetchlog calls are
# best-effort and may need tweaking for your jobsub_lite version — each
# is wrapped so a mismatch degrades gracefully instead of aborting. Read
# the CONFIG block and adjust before first use.
#######################################################################
set -uo pipefail

########################## CONFIG #####################################
EXPERIMENT=icarus
TAG=worktree-gOre-cc-sideband-category      # branch that has cc_sideband_category
GITUSER=hausnerh                            # fork owner (clone + tag validation)
USERNAME=${USER:-hhausner}

SEL_TOML=selection/toml/gOre_1g1p_sidebands.toml
SYS_TOML=systematics/toml/gOre_1g1p_sidebands.toml

BATCH_SIZE=20                               # CAF files per selection job
MEMORY_MB=4000                              # selection job memory
DISK_GB=2000                                # selection job disk (GB).
                                            #  WARNING: 2000 GB = 2 TB/job is very large.
                                            #  If jobs sit Idle (unmatchable), lower it.
LIFETIME=8h                                 # selection job lifetime
HELD_MEMORY_MB=8000                         # bumped resources when RELEASING a held job
HELD_DISK_GB=2000
POLL=300                                    # seconds between queue polls
CRED_WAIT_MAX=14400                         # seconds to wait for creds to come back (4h)

STAMP=$(date +%Y%m%d_%H%M%S)
OUT=/pnfs/icarus/scratch/users/$USERNAME/CCSidebandPlots
SEL_HADD=$PWD/build/output_gOre_1g1p.root
SYS_ROOT=$PWD/build/output_gOre_1g1p_sys.root
CONFIGS=(gOre_cc_Xg1p_stage1_datamc gOre_cc_Xg1p_stage2_datamc gOre_cc_Xg1p_stage3_datamc)
#######################################################################

# ---- argument parsing ------------------------------------------------
RESUME_PROJ=""
PLOTS_ONLY=0
SYS_ROOT_OVERRIDE=""
usage(){
  cat <<EOF
Usage: $0 [--resume PROJECT_DIR] [--plots-only [SYS_ROOT]] [--help]

  --resume PROJECT_DIR   Resume from an existing project on /pnfs. Accepts
                         either the project root or its output/ subdir.
                         Skips Stage 1 (selection submit + wait); proceeds
                         with whatever selection outputs exist under
                         PROJECT_DIR/output/. Useful when a run has
                         stragglers that never completed and you want to
                         push on with partial stats.

  --plots-only [SYS_ROOT]
                         Skip selection/hadd/systematics entirely and go
                         straight to Stage 4-5 (stage the sys ROOT + submit
                         the three plot jobs + retrieve). Reuses an existing
                         systematics ROOT — the local build one by default,
                         or the path you pass. Use this after a completed
                         systematics run so you do NOT recompute it.

Without a flag the script creates a fresh timestamped project and runs the
full pipeline end-to-end.
EOF
}
while [[ $# -gt 0 ]]; do
    case "$1" in
        --resume)      RESUME_PROJ="${2:-}"; shift 2 ;;
        --resume=*)    RESUME_PROJ="${1#*=}"; shift ;;
        --plots-only)  PLOTS_ONLY=1; shift
                       # optional non-flag arg = sys ROOT path
                       if [[ $# -gt 0 && "$1" != -* ]]; then SYS_ROOT_OVERRIDE="$1"; shift; fi ;;
        --plots-only=*) PLOTS_ONLY=1; SYS_ROOT_OVERRIDE="${1#*=}"; shift ;;
        -h|--help)     usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

# Resolve PROJ + adopt a matching stamp when resuming so the log/debug
# names line up with the original project directory.
if [[ -n "$RESUME_PROJ" ]]; then
    RESUME_PROJ=${RESUME_PROJ%/}
    [[ "$(basename "$RESUME_PROJ")" == "output" ]] && RESUME_PROJ=$(dirname "$RESUME_PROJ")
    PROJ="$RESUME_PROJ"
    if [[ "$(basename "$PROJ")" =~ _([0-9]{8}_[0-9]{6})$ ]]; then
        STAMP="${BASH_REMATCH[1]}_resume$(date +%H%M%S)"
    else
        STAMP="$(date +%Y%m%d_%H%M%S)_resume"
    fi
else
    PROJ=/pnfs/icarus/scratch/users/$USERNAME/gOre_cc_sideband_$STAMP
fi

# Reuse an explicit sys ROOT path if one was passed to --plots-only.
[[ -n "$SYS_ROOT_OVERRIDE" ]] && SYS_ROOT="$SYS_ROOT_OVERRIDE"

DEBUG_DIR=$PWD/cc_sideband_debug_$STAMP
LOGFILE=$PWD/run_cc_sideband_$STAMP.log
mkdir -p "$DEBUG_DIR" build

log(){ echo "[$(date -u +%H:%M:%S)] $*" | tee -a "$LOGFILE"; }
die(){ log "FATAL: $*"; exit 1; }

# ---- environment / preflight ----------------------------------------
setup_env(){
    set +u   # ups/cvmfs setup scripts reference unbound vars
    source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
    setup sbnana v10_01_02_01 -q e26:prof
    setup cmake  v3_27_4
    set -u
}

# Non-fatal, silent-on-success credential check. Returns 0 iff BOTH the
# kerberos ticket and a vault/bearer token look valid. On failure logs
# exactly what the user needs to run.
check_creds(){
    local kok=1 tok=1
    klist -s 2>/dev/null || kok=0
    if command -v httokendecode >/dev/null 2>&1; then
        httokendecode -H >/dev/null 2>&1 || tok=0
    else
        local uid; uid=$(id -u)
        [[ -s "${BEARER_TOKEN_FILE:-/run/user/$uid/bt_u$uid}" ]] || tok=0
    fi
    if [[ "$kok" -eq 0 || "$tok" -eq 0 ]]; then
        log "CREDS: kerberos=$([[ $kok -eq 1 ]] && echo OK || echo MISSING) bearer=$([[ $tok -eq 1 ]] && echo OK || echo MISSING)"
        log "CREDS: renew in another shell (or here) and this script will pick it up:"
        log "         kinit"
        log "         htgettoken -a htvaultprod.fnal.gov -i icarus"
        return 1
    fi
    return 0
}

# Loop until check_creds passes or timeout. Works under nohup — no
# interactive prompt; the user just renews creds in another shell.
wait_for_creds(){
    local max_wait=${1:-$CRED_WAIT_MAX} waited=0
    while ! check_creds; do
        if [[ "$waited" -ge "$max_wait" ]]; then
            log "CREDS: timed out waiting for valid credentials after ${max_wait}s"
            return 1
        fi
        log "CREDS: rechecking in 60s (waited ${waited}s of max ${max_wait}s)"
        sleep 60; waited=$((waited+60))
    done
    log "CREDS: valid credentials detected — resuming."
    return 0
}

ensure_creds(){ check_creds || wait_for_creds; }

preflight(){
    [[ -f "$SEL_TOML" ]] || die "run me from the medulla repo root (missing $SEL_TOML)"
    [[ -x build/systematics/run_systematics ]] \
        || log "WARN: build/systematics/run_systematics not found — build medulla before the systematics step."
    ensure_creds || die "no valid credentials at startup — nothing to do until they are renewed"
}

parse_jobid(){ grep -oE '[0-9]+\.[0-9]+@[A-Za-z0-9._-]+' | head -1; }

# hadd $SEL_HADD from a list of input files, capturing hadd's real exit
# status through the tail|tee pipe.
do_hadd(){ hadd -f "$SEL_HADD" "$@" 2>&1 | tail -5 | tee -a "$LOGFILE"; return "${PIPESTATUS[0]}"; }

# Report which CC event trees exist under events/full/ in the sys ROOT.
# TFile::ls() lists only top-level keys, so Get() each one explicitly.
# Echoes e.g. "s1=1 s2=1 s3=1".
sys_stage_report(){ # sysroot
    root -l -b -q -e "TFile f(\"$1\"); printf(\"s1=%d s2=%d s3=%d\n\", \
        f.Get(\"events/full/selected_cc_Xg1p_stage1\")!=0, \
        f.Get(\"events/full/selected_cc_Xg1p_stage2\")!=0, \
        f.Get(\"events/full/selected_cc_Xg1p_stage3\")!=0);" 2>/dev/null \
      | grep -oE 's[123]=[01]' | tr '\n' ' '
}

fetch_logs(){ # cluster schedd label
    log "$3: fetching job logs -> $DEBUG_DIR"
    jobsub_fetchlog -G "$EXPERIMENT" --jobid "$1@$2" --destdir  "$DEBUG_DIR" >/dev/null 2>&1 \
     || jobsub_fetchlog -G "$EXPERIMENT" --jobid "$1@$2" --dest-dir "$DEBUG_DIR" >/dev/null 2>&1 \
     || log "$3: jobsub_fetchlog failed — run: jobsub_fetchlog -G $EXPERIMENT --jobid $1@$2 --destdir $DEBUG_DIR"
}

# Poll the queue until this cluster drains. Releases held jobs with
# bumped resources; pauses (busy-wait, logging the fix) when credentials
# lapse. Drain is only declared after two consecutive empty reads so a
# transient jobsub_q blip does not end the wait early.
monitor_cluster(){ # cluster schedd label
    local cluster="$1" schedd="$2" label="$3"
    local q active held jid stable=0
    log "$label: monitoring cluster $cluster@$schedd (poll every ${POLL}s)"
    while true; do
        if ! check_creds; then
            wait_for_creds || { log "$label: giving up on this cluster (no creds)"; return 1; }
        fi
        if ! q=$(jobsub_q -G "$EXPERIMENT" 2>/dev/null); then
            log "$label: jobsub_q failed this poll; retrying in ${POLL}s"
            sleep "$POLL"; continue
        fi
        active=$(printf '%s\n' "$q" | grep -Ec "(^|[[:space:]])${cluster}\.[0-9]+@")
        if [[ "$active" -eq 0 ]]; then
            stable=$((stable + 1))
            [[ "$stable" -ge 2 ]] && { log "$label: queue drained"; return 0; }
            sleep "$POLL"; continue
        fi
        stable=0
        held=$(printf '%s\n' "$q" | awk -v c="$cluster" '$1 ~ (c "\\.[0-9]+@") && $5=="H"{print $1}')
        if [[ -n "$held" ]]; then
            for jid in $held; do
                log "$label: releasing held $jid (mem=${HELD_MEMORY_MB}MB disk=${HELD_DISK_GB}GB)"
                condor_qedit -name "$schedd" "${jid%@*}" RequestMemory "$HELD_MEMORY_MB"             >/dev/null 2>&1 || true
                condor_qedit -name "$schedd" "${jid%@*}" RequestDisk   "$((HELD_DISK_GB*1024*1024))" >/dev/null 2>&1 || true
                jobsub_release -G "$EXPERIMENT" --jobid "$jid" >/dev/null 2>&1 || true
            done
        fi
        log "$label: $active active$([[ -n "$held" ]] && echo ', released held')"
        sleep "$POLL"
    done
}

# Best-effort: read the expected job count from project.db. dCache can
# lag right after create, so retry a few times and just set NJOBS empty
# if we cannot read it (a subsequent step just skips its count check).
read_project_njobs(){
    NJOBS=""
    local db="$DEBUG_DIR/project.db" attempt
    for attempt in 1 2 3 4 5; do
        if ifdh cp "$PROJ/project.db" "$db" >/dev/null 2>&1; then
            NJOBS=$(sqlite3 "$db" "SELECT COUNT(*) FROM jobs;" 2>/dev/null)
            [[ "${NJOBS:-0}" -gt 0 ]] 2>/dev/null && return 0
        fi
        NJOBS=""
        sleep 4
    done
    return 1
}

# ---- pipeline phases -------------------------------------------------
submit_selection(){
    log "=== Stage 1: grid selection ==="
    # Robust, LOCAL check that the checkout has the stage2/3 trees (the stored
    # project config is derived from this toml). Avoids reading project.db back
    # off /pnfs, which dCache has not made consistent yet right after create.
    grep -q 'selected_cc_Xg1p_stage3' "$SEL_TOML" \
        || die "local $SEL_TOML lacks the stage3 tree — git pull the $TAG branch first"

    log "Creating project $PROJ (batch-size $BATCH_SIZE)"
    python3 batch/medulla.py --experiment "$EXPERIMENT" --project-dir "$PROJ" \
        --create-project --toml "$SEL_TOML" --batch-size "$BATCH_SIZE" \
        --tag "$TAG" --gituser "$GITUSER" 2>&1 | tee -a "$LOGFILE" \
        || die "create-project failed"

    if read_project_njobs; then
        log "Project has $NJOBS jobs."
    else
        log "WARN: could not read project.db job count (dCache lag?); continuing without an expected-count check."
    fi

    log "Launching ${NJOBS:-?} jobs (mem=${MEMORY_MB}MB disk=${DISK_GB}GB lifetime=$LIFETIME)"
    local out jid
    out=$(yes | python3 batch/medulla.py --experiment "$EXPERIMENT" --project-dir "$PROJ" \
              --tag "$TAG" --gituser "$GITUSER" --launch-jobs \
              --memory "$MEMORY_MB" --disk "$DISK_GB" --lifetime "$LIFETIME" 2>&1)
    echo "$out" | tee -a "$LOGFILE"
    jid=$(echo "$out" | parse_jobid)
    [[ -n "$jid" ]] || die "could not parse jobsub id from launch output"
    SEL_SCHEDD=${jid#*@}
    SEL_CLUSTER=${jid%@*}; SEL_CLUSTER=${SEL_CLUSTER%.*}
    log "Selection cluster = $SEL_CLUSTER @ $SEL_SCHEDD"
}

gather_and_merge(){
    if [[ -n "${SEL_CLUSTER:-}" ]]; then
        monitor_cluster "$SEL_CLUSTER" "$SEL_SCHEDD" selection
    fi

    log "=== Stage 2: gather + hadd ==="
    [[ -n "${SEL_CLUSTER:-}" ]] || log "(resume mode: proceeding with outputs already on /pnfs)"
    ensure_creds || die "no credentials to read $PROJ/output"

    # If we did not come from submit_selection (resume), try to learn the
    # expected job count so we can log a completeness ratio.
    [[ -z "${NJOBS:-}" ]] && read_project_njobs || true

    local outs
    mapfile -t outs < <(ifdh ls "$PROJ/output" 2>/dev/null | grep -oE 'output_jobid[0-9]+\.root' | sort -u)
    local done=${#outs[@]}
    log "Selection produced $done / ${NJOBS:-?} job outputs."
    if [[ -n "${NJOBS:-}" && "$done" -lt "$NJOBS" ]]; then
        log "WARN: $((NJOBS - done)) selection job(s) missing output — continuing with partial stats."
        [[ -n "${SEL_CLUSTER:-}" ]] && fetch_logs "$SEL_CLUSTER" "$SEL_SCHEDD" selection
    fi
    [[ "$done" -gt 0 ]] || die "no selection outputs produced — see logs in $DEBUG_DIR"

    # Prefer the NFS /pnfs path (fast, no local copy; works with valid
    # creds). Fall back to ifdh-cp'ing each output local and hadd'ing those
    # if the NFS read is denied (dCache perms / xrootd path quirks).
    local nfs=() f
    for f in "${outs[@]}"; do nfs+=("$PROJ/output/$f"); done
    log "hadd -> $SEL_HADD ($done files) — via NFS /pnfs"
    if do_hadd "${nfs[@]}" && [[ -s "$SEL_HADD" ]]; then
        log "hadd via NFS /pnfs succeeded"
    else
        log "NFS hadd failed — falling back to ifdh-cp local then hadd"
        local tmp="build/_haddtmp_$STAMP" loc=()
        mkdir -p "$tmp"
        for f in "${outs[@]}"; do
            ifdh cp "$PROJ/output/$f" "$tmp/$f" >/dev/null 2>&1 || die "ifdh cp $f failed (creds?)"
            loc+=("$tmp/$f")
        done
        do_hadd "${loc[@]}" && [[ -s "$SEL_HADD" ]] || { rm -rf "$tmp"; die "local hadd failed"; }
        rm -rf "$tmp"
    fi
    [[ -s "$SEL_HADD" ]] || die "hadd produced no $SEL_HADD"
}

do_systematics(){
    log "=== Stage 3: systematics ==="
    ( cd build && ./systematics/run_systematics "../$SYS_TOML" ) 2>&1 | tail -15 | tee -a "$LOGFILE" \
        || die "run_systematics failed"
    [[ -s "$SYS_ROOT" ]] || die "run_systematics produced no $SYS_ROOT"
    local chk; chk=$(sys_stage_report "$SYS_ROOT")
    log "sys ROOT CC event trees (events/full): ${chk:-unknown}"
    echo "$chk" | grep -q 's3=1' \
        || die "stage3 tree not in events/full of $SYS_ROOT — systematics did not produce it"
}

submit_plots(){
    log "=== Stage 4: stage sys ROOT + plot stage1/2/3 on grid ==="
    [[ -s "$SYS_ROOT" ]] || die "sys ROOT not found: $SYS_ROOT (run systematics, or pass --plots-only <path>)"
    log "sys ROOT CC event trees (events/full): $(sys_stage_report "$SYS_ROOT" || echo unknown)"
    ensure_creds || die "no credentials to stage sys ROOT to $OUT"
    ifdh cp "$SYS_ROOT" "$OUT/output_gOre_1g1p_sys.root" >/dev/null 2>&1 \
        || die "failed to stage $SYS_ROOT to $OUT"
    log "Staged sys ROOT to $OUT/output_gOre_1g1p_sys.root"

    local plot_jobs=() cfg out jid pschedd pcluster
    for cfg in "${CONFIGS[@]}"; do
        log "Launching plot job: $cfg"
        out=$(yes | ./batch/launch_spineplot.sh \
                --input="$OUT/output_gOre_1g1p_sys.root" --output="$OUT/$cfg" \
                --config="$cfg" --tag="$TAG" --gituser="$GITUSER" 2>&1)
        echo "$out" | tee -a "$LOGFILE"
        jid=$(echo "$out" | parse_jobid)
        [[ -z "$jid" ]] && { log "WARN: could not parse jobid for $cfg"; continue; }
        pschedd=${jid#*@}
        pcluster=${jid%@*}; pcluster=${pcluster%.*}
        plot_jobs+=("$pcluster $pschedd $cfg")
    done

    local rec
    for rec in "${plot_jobs[@]}"; do
        set -- $rec   # cluster schedd cfg
        monitor_cluster "$1" "$2" "plot:$3"
        if ifdh ls "$OUT/$3" 2>/dev/null | grep -qE '\.(png|pdf)$'; then
            log "plot:$3 produced figures in $OUT/$3"
        else
            log "WARN: plot:$3 produced no figures — fetching logs"
            fetch_logs "$1" "$2" "plot:$3"
        fi
    done
}

retrieve(){
    log "=== Stage 5: retrieve figures ==="
    ensure_creds || die "no credentials to pull figures back"
    local dst="$PWD/CCSidebandPlots_$STAMP" cfg
    for cfg in "${CONFIGS[@]}"; do
        mkdir -p "$dst/$cfg"
        ifdh cp -D "$OUT/$cfg"/*.png "$dst/$cfg/" >/dev/null 2>&1 || true
        ifdh cp -D "$OUT/$cfg"/*.pdf "$dst/$cfg/" >/dev/null 2>&1 || true
    done
    log "Figures -> $dst ; debug logs (if any) -> $DEBUG_DIR ; full log -> $LOGFILE"
}

main(){
    log "CC sideband run starting (stamp $STAMP)"
    setup_env
    preflight
    if [[ "$PLOTS_ONLY" -eq 1 ]]; then
        log "Plots-only mode: reusing sys ROOT $SYS_ROOT (skipping selection/hadd/systematics)"
        submit_plots
        retrieve
        log "DONE."
        return
    fi
    [[ -n "$RESUME_PROJ" ]] && log "Resume mode: PROJ=$PROJ (skipping Stage 1)"
    if [[ -z "$RESUME_PROJ" ]]; then
        submit_selection
    fi
    gather_and_merge
    do_systematics
    submit_plots
    retrieve
    log "DONE."
}
main "$@"
