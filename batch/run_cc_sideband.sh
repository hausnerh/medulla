#!/bin/bash
#######################################################################
# One-shot, unattended CC Xγ+1p sideband production — SPLIT-ENVIRONMENT.
#
# On EL9 gpvms jobsub_lite runs NATIVELY but the medulla/SPINE runtime
# (ifdh, root, hadd, run_systematics) only works inside the SL7 apptainer
# container, and the two cannot coexist in one shell. So this script runs
# NATIVELY (jobsub_submit/_q/_release/_fetchlog, launch_spineplot.sh,
# medulla.py --launch-jobs all bare) and shells into SL7 via SL7RUN for
# every ifdh / hadd / run_systematics / root / sqlite3 / --create-project
# step. The long monitor loop is native; it only dips into SL7 at phase
# boundaries.
#
#   selection (grid) -> monitor -> hadd -> systematics -> stage
#   -> plots stage1/2/3 (grid) -> monitor -> fetch figures
#
# Launch and leave it (run natively, NOT inside sl7_container):
#   nohup ./batch/run_cc_sideband.sh > run_cc_sideband.nohup 2>&1 &
#
# FIRST TIME: validate the container plumbing before a real run:
#   ./batch/run_cc_sideband.sh --check-env
#
# Requires a valid token + kerberos ticket (kinit ; htgettoken -a
# htvaultprod.fnal.gov -i icarus) — both are visible inside SL7 via the
# /run/user bind.
#
# NOTE: jobsub_q / condor_qedit / jobsub_fetchlog parsing is best-effort
# and may need tweaking for your jobsub_lite version; each is wrapped so a
# mismatch degrades gracefully. Read the CONFIG block before first use.
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
DISK_GB=2000                                # selection job disk (GB); 2 TB is very
                                            #  large — lower if jobs stay Idle.
LIFETIME=8h                                 # selection job lifetime
HELD_MEMORY_MB=8000                         # bumped resources when RELEASING a held job
HELD_DISK_GB=2000
POLL=300                                    # seconds between queue polls
CRED_WAIT_MAX=14400                         # seconds to wait for creds to return (4h)

# ---- SL7 container (apptainer) for ifdh / hadd / run_systematics / root ----
APPTAINER=/cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer
SL7_IMAGE=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:jsl
# Binds copied from start_SL7dev_jsl.sh (adjust if yours differ).
SL7_BINDS=(--pid --ipc
    -B "/etc/hosts,/run/user/$(id -u),/tmp,/opt,/cvmfs,/exp/,/pnfs/,/nashome"
    -B "/etc/profile.d/jobsub_lite.sh,/etc/condor/,/etc/grid-security/")
# Command(s) run inside SL7 to set up the medulla/SPINE runtime (ifdh, root,
# hadd, run_systematics). If `setup_spine` isn't defined in a non-interactive
# login shell, change this to e.g.
#   'source /exp/'"$(id -ng)"'/data/users/vito/podman/profile_jsl && setup_spine'
SL7_SETUP='setup_spine'
# Setup to run in the NATIVE shell before jobsub (empty = jobsub already on PATH).
NATIVE_SETUP=''

STAMP=$(date +%Y%m%d_%H%M%S)
OUT=/pnfs/icarus/scratch/users/$USERNAME/CCSidebandPlots
SEL_HADD=$PWD/build/output_gOre_1g1p.root
SYS_ROOT=$PWD/build/output_gOre_1g1p_sys.root
CONFIGS=(gOre_cc_Xg1p_stage1_datamc gOre_cc_Xg1p_stage2_datamc gOre_cc_Xg1p_stage3_datamc)
REPO="$PWD"
#######################################################################

# Run a single command string inside the SL7 container, from the repo dir,
# with the SPINE env sourced. Everything ifdh/root/hadd/systematics goes here.
SL7RUN(){
    "$APPTAINER" exec "${SL7_BINDS[@]}" --home "$HOME":"$HOME" --pwd "$REPO" \
        "$SL7_IMAGE" /bin/bash -lc "{ $SL7_SETUP ; } >/dev/null 2>&1 ; cd '$REPO' && { $1 ; }"
}
# Run a native command with any NATIVE_SETUP applied (for jobsub).
NAT(){ if [[ -n "$NATIVE_SETUP" ]]; then bash -lc "{ $NATIVE_SETUP ; } >/dev/null 2>&1 ; $1"; else bash -lc "$1"; fi; }

# ---- argument parsing ------------------------------------------------
RESUME_PROJ=""; PLOTS_ONLY=0; SYS_ROOT_OVERRIDE=""; CHECK_ENV=0
usage(){
  cat <<EOF
Usage: $0 [--resume PROJECT_DIR] [--plots-only [SYS_ROOT]] [--check-env] [--help]

  --resume PROJECT_DIR   Skip Stage 1; continue from an existing project's
                         outputs on /pnfs (root or its output/ subdir).
  --plots-only [SYS_ROOT]
                         Skip selection/hadd/systematics; stage an existing
                         sys ROOT (local build one by default) and plot.
  --check-env            Run only the SL7/native plumbing smoke tests and exit.
                         Run this once before your first real run.

Run this NATIVELY (not inside sl7_container). jobsub is native; ifdh/hadd/
run_systematics/root run in SL7 via the SL7RUN wrapper (see CONFIG block).
EOF
}
while [[ $# -gt 0 ]]; do
    case "$1" in
        --resume)      RESUME_PROJ="${2:-}"; shift 2 ;;
        --resume=*)    RESUME_PROJ="${1#*=}"; shift ;;
        --plots-only)  PLOTS_ONLY=1; shift
                       if [[ $# -gt 0 && "$1" != -* ]]; then SYS_ROOT_OVERRIDE="$1"; shift; fi ;;
        --plots-only=*) PLOTS_ONLY=1; SYS_ROOT_OVERRIDE="${1#*=}"; shift ;;
        --check-env)   CHECK_ENV=1; shift ;;
        -h|--help)     usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ -n "$RESUME_PROJ" ]]; then
    RESUME_PROJ=${RESUME_PROJ%/}
    [[ "$(basename "$RESUME_PROJ")" == "output" ]] && RESUME_PROJ=$(dirname "$RESUME_PROJ")
    PROJ="$RESUME_PROJ"
    if [[ "$(basename "$PROJ")" =~ _([0-9]{8}_[0-9]{6})$ ]]; then STAMP="${BASH_REMATCH[1]}_resume$(date +%H%M%S)"
    else STAMP="$(date +%Y%m%d_%H%M%S)_resume"; fi
else
    PROJ=/pnfs/icarus/scratch/users/$USERNAME/gOre_cc_sideband_$STAMP
fi
[[ -n "$SYS_ROOT_OVERRIDE" ]] && SYS_ROOT="$SYS_ROOT_OVERRIDE"

DEBUG_DIR=$PWD/cc_sideband_debug_$STAMP
LOGFILE=$PWD/run_cc_sideband_$STAMP.log
mkdir -p "$DEBUG_DIR" build

log(){ echo "[$(date -u +%H:%M:%S)] $*" | tee -a "$LOGFILE"; }
die(){ log "FATAL: $*"; exit 1; }

# ---- credentials (native) --------------------------------------------
check_creds(){
    local kok=1 tok=1 uid; uid=$(id -u)
    klist -s 2>/dev/null || kok=0
    [[ -s "${BEARER_TOKEN_FILE:-/run/user/$uid/bt_u$uid}" ]] || tok=0
    if [[ "$kok" -eq 0 || "$tok" -eq 0 ]]; then
        log "CREDS: kerberos=$([[ $kok -eq 1 ]] && echo OK || echo MISSING) bearer=$([[ $tok -eq 1 ]] && echo OK || echo MISSING)"
        log "CREDS: renew in another shell and this script will pick it up:"
        log "         kinit ; htgettoken -a htvaultprod.fnal.gov -i icarus"
        return 1
    fi
    return 0
}
wait_for_creds(){
    local waited=0
    while ! check_creds; do
        [[ "$waited" -ge "$CRED_WAIT_MAX" ]] && { log "CREDS: timed out after ${CRED_WAIT_MAX}s"; return 1; }
        log "CREDS: rechecking in 60s (waited ${waited}s)"; sleep 60; waited=$((waited+60))
    done
    log "CREDS: valid — resuming."; return 0
}
ensure_creds(){ check_creds || wait_for_creds; }

# ---- helpers ---------------------------------------------------------
parse_jobid(){ grep -oE '[0-9]+\.[0-9]+@[A-Za-z0-9._-]+' | head -1; }

# hadd (SL7)
do_hadd(){ SL7RUN "hadd -f '$SEL_HADD' $*" 2>&1 | tail -5 | tee -a "$LOGFILE"; return "${PIPESTATUS[0]}"; }

# report which CC event trees exist under events/full/ (SL7 root via a temp macro
# to avoid nested-quote hell). Echoes "s1=1 s2=1 s3=1".
sys_stage_report(){ # sysroot
    local id=cc$$_$RANDOM mac=/tmp/$id.C
    printf 'void %s(){TFile f("%s");printf("s1=%%d s2=%%d s3=%%d\\n",f.Get("events/full/selected_cc_Xg1p_stage1")!=0,f.Get("events/full/selected_cc_Xg1p_stage2")!=0,f.Get("events/full/selected_cc_Xg1p_stage3")!=0);}\n' "$id" "$1" > "$mac"
    SL7RUN "root -l -b -q '$mac'" 2>/dev/null | grep -oE 's[123]=[01]' | tr '\n' ' ' || true
    rm -f "$mac"
}
# entry count of a tree (-1 if absent). SL7 root via temp macro.
tree_entries(){ # sysroot treepath
    local id=te$$_$RANDOM mac=/tmp/$id.C out
    printf 'void %s(){TFile f("%s");TTree*t=(TTree*)f.Get("%s");printf("N=%%lld\\n",t?t->GetEntries():-1);}\n' "$id" "$1" "$2" > "$mac"
    out=$(SL7RUN "root -l -b -q '$mac'" 2>/dev/null | grep -oE 'N=-?[0-9]+' | head -1 | cut -d= -f2)
    rm -f "$mac"; echo "${out:-}"
}
# best-effort job count from project.db (SL7 ifdh+sqlite3), retried past dCache lag.
read_project_njobs(){
    NJOBS=""; local db="$DEBUG_DIR/project.db" attempt
    for attempt in 1 2 3 4 5; do
        if SL7RUN "ifdh cp '$PROJ/project.db' '$db'" >/dev/null 2>&1; then
            NJOBS=$(SL7RUN "sqlite3 '$db' 'SELECT COUNT(*) FROM jobs;'" 2>/dev/null | tr -dc '0-9')
            [[ "${NJOBS:-0}" -gt 0 ]] 2>/dev/null && return 0
        fi
        NJOBS=""; sleep 4
    done
    return 1
}
# fetch job logs (native jobsub_fetchlog)
fetch_logs(){ # cluster schedd label
    log "$3: fetching job logs -> $DEBUG_DIR"
    NAT "jobsub_fetchlog -G $EXPERIMENT --jobid '$1@$2' --destdir '$DEBUG_DIR'" >/dev/null 2>&1 \
     || NAT "jobsub_fetchlog -G $EXPERIMENT --jobid '$1@$2' --dest-dir '$DEBUG_DIR'" >/dev/null 2>&1 \
     || log "$3: jobsub_fetchlog failed — run it by hand for $1@$2"
}

# ---- monitor (native jobsub) -----------------------------------------
monitor_cluster(){ # cluster schedd label
    local cluster="$1" schedd="$2" label="$3" q active held jid stable=0
    log "$label: monitoring cluster $cluster@$schedd (poll ${POLL}s)"
    while true; do
        if ! check_creds; then wait_for_creds || { log "$label: giving up (no creds)"; return 1; }; fi
        if ! q=$(NAT "jobsub_q -G $EXPERIMENT" 2>/dev/null); then
            log "$label: jobsub_q failed; retry in ${POLL}s"; sleep "$POLL"; continue
        fi
        active=$(printf '%s\n' "$q" | grep -Ec "(^|[[:space:]])${cluster}\.[0-9]+@")
        if [[ "$active" -eq 0 ]]; then
            stable=$((stable+1)); [[ "$stable" -ge 2 ]] && { log "$label: queue drained"; return 0; }
            sleep "$POLL"; continue
        fi
        stable=0
        held=$(printf '%s\n' "$q" | awk -v c="$cluster" '$1 ~ (c "\\.[0-9]+@") && $5=="H"{print $1}')
        if [[ -n "$held" ]]; then
            for jid in $held; do
                log "$label: releasing held $jid (mem=${HELD_MEMORY_MB}MB disk=${HELD_DISK_GB}GB)"
                NAT "condor_qedit -name '$schedd' '${jid%@*}' RequestMemory $HELD_MEMORY_MB" >/dev/null 2>&1 || true
                NAT "condor_qedit -name '$schedd' '${jid%@*}' RequestDisk $((HELD_DISK_GB*1024*1024))" >/dev/null 2>&1 || true
                NAT "jobsub_release -G $EXPERIMENT --jobid '$jid'" >/dev/null 2>&1 || true
            done
        fi
        log "$label: $active active$([[ -n "$held" ]] && echo ', released held')"
        sleep "$POLL"
    done
}

# ---- pipeline phases -------------------------------------------------
preflight(){
    [[ -f "$SEL_TOML" ]] || die "run me from the medulla repo root (missing $SEL_TOML)"
    [[ -x "$APPTAINER" ]] || die "apptainer not found at $APPTAINER — fix APPTAINER in CONFIG"
    ensure_creds || die "no valid credentials at startup"
}

submit_selection(){
    log "=== Stage 1: grid selection ==="
    grep -q 'selected_cc_Xg1p_stage3' "$SEL_TOML" \
        || die "local $SEL_TOML lacks the stage3 tree — git pull the $TAG branch first"

    log "Creating project $PROJ (batch-size $BATCH_SIZE) [SL7]"
    SL7RUN "python3 batch/medulla.py --experiment $EXPERIMENT --project-dir '$PROJ' \
        --create-project --toml '$SEL_TOML' --batch-size $BATCH_SIZE \
        --tag '$TAG' --gituser '$GITUSER'" 2>&1 | tee -a "$LOGFILE" \
        || die "create-project failed"

    if read_project_njobs; then log "Project has $NJOBS jobs."
    else log "WARN: could not read project.db job count (dCache lag?); continuing."; fi

    log "Launching ${NJOBS:-?} jobs [native jobsub] (mem=${MEMORY_MB}MB disk=${DISK_GB}GB lifetime=$LIFETIME)"
    local out jid
    out=$(NAT "yes | python3 batch/medulla.py --experiment $EXPERIMENT --project-dir '$PROJ' \
              --tag '$TAG' --gituser '$GITUSER' --launch-jobs \
              --memory $MEMORY_MB --disk $DISK_GB --lifetime $LIFETIME" 2>&1)
    echo "$out" | tee -a "$LOGFILE"
    jid=$(echo "$out" | parse_jobid)
    [[ -n "$jid" ]] || die "could not parse jobsub id from launch output"
    SEL_SCHEDD=${jid#*@}; SEL_CLUSTER=${jid%@*}; SEL_CLUSTER=${SEL_CLUSTER%.*}
    log "Selection cluster = $SEL_CLUSTER @ $SEL_SCHEDD"
}

gather_and_merge(){
    [[ -n "${SEL_CLUSTER:-}" ]] && monitor_cluster "$SEL_CLUSTER" "$SEL_SCHEDD" selection
    log "=== Stage 2: gather + hadd ==="
    [[ -n "${SEL_CLUSTER:-}" ]] || log "(resume mode: using outputs already on /pnfs)"
    ensure_creds || die "no credentials to read $PROJ/output"
    [[ -z "${NJOBS:-}" ]] && read_project_njobs || true

    local outs
    mapfile -t outs < <(SL7RUN "ifdh ls '$PROJ/output'" 2>/dev/null | grep -oE 'output_jobid[0-9]+\.root' | sort -u)
    local done=${#outs[@]}
    log "Selection produced $done / ${NJOBS:-?} job outputs."
    if [[ -n "${NJOBS:-}" && "$done" -lt "$NJOBS" ]]; then
        log "WARN: $((NJOBS-done)) job(s) missing output — continuing with partial stats."
        [[ -n "${SEL_CLUSTER:-}" ]] && fetch_logs "$SEL_CLUSTER" "$SEL_SCHEDD" selection
    fi
    [[ "$done" -gt 0 ]] || die "no selection outputs — see $DEBUG_DIR"

    # ifdh-cp each output local, then hadd (all in SL7). Avoids NFS/xrootd perms.
    local tmp="build/_haddtmp_$STAMP" loc=() f
    SL7RUN "mkdir -p '$tmp'"
    for f in "${outs[@]}"; do
        SL7RUN "ifdh cp '$PROJ/output/$f' '$tmp/$f'" >/dev/null 2>&1 || die "ifdh cp $f failed (creds?)"
        loc+=("$tmp/$f")
    done
    log "hadd -> $SEL_HADD ($done files) [SL7]"
    do_hadd "${loc[@]/#/$REPO/}" || { SL7RUN "rm -rf '$tmp'"; die "hadd failed"; }
    SL7RUN "rm -rf '$tmp'"
    [[ -s "$SEL_HADD" ]] || die "hadd produced no $SEL_HADD"
}

do_systematics(){
    log "=== Stage 3: systematics [SL7] ==="
    SL7RUN "cd build && ./systematics/run_systematics '../$SYS_TOML'" 2>&1 | tail -15 | tee -a "$LOGFILE" \
        || die "run_systematics failed"
    [[ -s "$SYS_ROOT" ]] || die "run_systematics produced no $SYS_ROOT"
    local chk; chk=$(sys_stage_report "$SYS_ROOT")
    log "sys ROOT CC event trees (events/full): ${chk:-unknown}"
    echo "$chk" | grep -q 's3=1' || die "stage3 tree not in events/full of $SYS_ROOT"
}

submit_plots(){
    log "=== Stage 4: stage sys ROOT + plot stage1/2/3 on grid ==="
    [[ -s "$SYS_ROOT" ]] || die "sys ROOT not found: $SYS_ROOT (run systematics, or --plots-only <path>)"
    log "sys ROOT CC event trees (events/full): $(sys_stage_report "$SYS_ROOT")"
    ensure_creds || die "no credentials to stage sys ROOT"
    SL7RUN "ifdh mkdir_p '$OUT'" >/dev/null 2>&1 || true
    SL7RUN "ifdh rm '$OUT/output_gOre_1g1p_sys.root'" >/dev/null 2>&1 || true
    if ! SL7RUN "ifdh cp '$SYS_ROOT' '$OUT/output_gOre_1g1p_sys.root'" 2> "$DEBUG_DIR/stage_ifdh.err"; then
        log "ifdh cp error: $(tail -3 "$DEBUG_DIR/stage_ifdh.err" 2>/dev/null | tr '\n' ' ')"
        die "failed to stage $SYS_ROOT to $OUT (see $DEBUG_DIR/stage_ifdh.err)"
    fi
    log "Staged sys ROOT to $OUT/output_gOre_1g1p_sys.root"

    local plot_jobs=() cfg out jid pschedd pcluster stg n
    for cfg in "${CONFIGS[@]}"; do
        stg=$(echo "$cfg" | grep -oE 'stage[0-9]+')
        n=$(tree_entries "$SYS_ROOT" "events/full/selected_cc_Xg1p_$stg")
        if [[ "${n:-0}" -le 0 ]] 2>/dev/null; then
            log "SKIP $cfg: events/full/selected_cc_Xg1p_$stg has ${n:-0} entries — nothing to plot."; continue
        fi
        log "Launching plot job: $cfg ($n MC events) [native jobsub]"
        out=$(NAT "yes | ./batch/launch_spineplot.sh --input='$OUT/output_gOre_1g1p_sys.root' \
                --output='$OUT/$cfg' --config='$cfg' --tag='$TAG' --gituser='$GITUSER'" 2>&1)
        echo "$out" | tee -a "$LOGFILE"
        jid=$(echo "$out" | parse_jobid)
        [[ -z "$jid" ]] && { log "WARN: could not parse jobid for $cfg"; continue; }
        pschedd=${jid#*@}; pcluster=${jid%@*}; pcluster=${pcluster%.*}
        plot_jobs+=("$pcluster $pschedd $cfg")
    done

    local rec
    for rec in "${plot_jobs[@]}"; do
        set -- $rec
        monitor_cluster "$1" "$2" "plot:$3"
        if SL7RUN "ifdh ls '$OUT/$3'" 2>/dev/null | grep -qE '\.(png|pdf)$'; then
            log "plot:$3 produced figures in $OUT/$3"
        else
            log "WARN: plot:$3 produced no figures — fetching logs"; fetch_logs "$1" "$2" "plot:$3"
        fi
    done
}

retrieve(){
    log "=== Stage 5: retrieve figures ==="
    ensure_creds || die "no credentials to pull figures"
    local dst="$REPO/CCSidebandPlots_$STAMP" cfg fp bn
    mkdir -p "$dst"
    for cfg in "${CONFIGS[@]}"; do
        mkdir -p "$dst/$cfg"
        # enumerate via ifdh ls (SL7), copy each by basename — no shell glob on /pnfs.
        while IFS= read -r fp; do
            [[ -z "$fp" ]] && continue
            bn=$(basename "$fp")
            SL7RUN "ifdh cp '$OUT/$cfg/$bn' '$dst/$cfg/$bn'" >/dev/null 2>&1 || true
        done < <(SL7RUN "ifdh ls '$OUT/$cfg'" 2>/dev/null | grep -E '\.(png|pdf)$')
    done
    log "Figures -> $dst ; debug logs -> $DEBUG_DIR ; full log -> $LOGFILE"
}

# ---- environment smoke test ------------------------------------------
check_env(){
    log "=== --check-env: validating native + SL7 plumbing ==="
    local ok=1
    log "[native] jobsub_q ..."
    NAT "jobsub_q -G $EXPERIMENT >/dev/null 2>&1" && log "  native jobsub_q OK" || { log "  native jobsub_q FAILED (NATIVE_SETUP?)"; ok=0; }
    log "[native] python3 has toml (needed by medulla.py --launch-jobs) ..."
    NAT "python3 -c 'import toml' >/dev/null 2>&1" && log "  native python3 toml OK" || { log "  native python3 lacks 'toml' — run from your .venv, or pip install toml"; ok=0; }
    log "[SL7] setup_spine + which ifdh/root/hadd ..."
    SL7RUN "command -v ifdh && command -v root && command -v hadd" >/dev/null 2>&1 \
        && log "  SL7 ifdh/root/hadd OK" || { log "  SL7 tools NOT found (SL7_SETUP / image?)"; ok=0; }
    log "[SL7] ifdh ls scratch ..."
    SL7RUN "ifdh ls /pnfs/icarus/scratch/users/$USERNAME/ >/dev/null 2>&1" \
        && log "  SL7 ifdh ls OK" || log "  SL7 ifdh ls FAILED (creds? bind?)"
    [[ "$ok" -eq 1 ]] && log "check-env: PASS" || log "check-env: FAIL — fix the items above before a real run"
    return $((1-ok))
}

main(){
    log "CC sideband run starting (stamp $STAMP)"
    preflight
    if [[ "$CHECK_ENV" -eq 1 ]]; then check_env; exit $?; fi
    if [[ "$PLOTS_ONLY" -eq 1 ]]; then
        log "Plots-only mode: reusing sys ROOT $SYS_ROOT"
        submit_plots; retrieve; log "DONE."; return
    fi
    [[ -n "$RESUME_PROJ" ]] && log "Resume mode: PROJ=$PROJ (skipping Stage 1)"
    [[ -z "$RESUME_PROJ" ]] && submit_selection
    gather_and_merge
    do_systematics
    submit_plots
    retrieve
    log "DONE."
}
main "$@"
