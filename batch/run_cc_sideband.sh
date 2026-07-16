#!/bin/bash
#######################################################################
# One-shot, unattended CC Xγ+1p sideband production.
#
#   selection (grid) -> monitor -> hadd -> systematics -> stage
#   -> plots stage1/2/3 (grid) -> monitor -> fetch figures
#
# Monitors the grid itself: polls the queue, RELEASES held jobs with
# bumped memory/disk, and on any failure FETCHES the job logs into a
# local debug dir. Designed to be launched and left alone:
#
#   nohup ./batch/run_cc_sideband.sh > run_cc_sideband.nohup 2>&1 &
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

STAMP=$(date +%Y%m%d_%H%M%S)
PROJ=/pnfs/icarus/scratch/users/$USERNAME/gOre_cc_sideband_$STAMP
OUT=/pnfs/icarus/scratch/users/$USERNAME/CCSidebandPlots
SEL_HADD=$PWD/build/output_gOre_1g1p.root
SYS_ROOT=$PWD/build/output_gOre_1g1p_sys.root
DEBUG_DIR=$PWD/cc_sideband_debug_$STAMP
LOGFILE=$PWD/run_cc_sideband_$STAMP.log
CONFIGS=(gOre_cc_Xg1p_stage1_datamc gOre_cc_Xg1p_stage2_datamc gOre_cc_Xg1p_stage3_datamc)
# xrootd door used to hadd per-job outputs off dCache without the NFS
# /pnfs data path (which can EPERM). Adjust if your door differs.
XROOTD_DOOR="root://fndca1.fnal.gov:1094"
#######################################################################

mkdir -p "$DEBUG_DIR" build
log(){ echo "[$(date -u +%H:%M:%S)] $*" | tee -a "$LOGFILE"; }
die(){ log "FATAL: $*"; exit 1; }

setup_env(){
    set +u   # ups/cvmfs setup scripts reference unbound vars
    source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
    setup sbnana v10_01_02_01 -q e26:prof
    setup cmake  v3_27_4
    set -u
}

preflight(){
    [[ -f "$SEL_TOML" ]] || die "run me from the medulla repo root (missing $SEL_TOML)"
    [[ -x build/systematics/run_systematics ]] \
        || log "WARN: build/systematics/run_systematics not found — build medulla before the systematics step."
    if ! klist -s 2>/dev/null; then
        log "WARN: no valid kerberos ticket (klist -s failed); /pnfs may EPERM."
        log "      Fix: kinit ; htgettoken -a htvaultprod.fnal.gov -i icarus"
    fi
}

parse_jobid(){ grep -oE '[0-9]+\.[0-9]+@[A-Za-z0-9._-]+' | head -1; }

fetch_logs(){ # cluster schedd label
    log "$3: fetching job logs -> $DEBUG_DIR"
    jobsub_fetchlog -G "$EXPERIMENT" --jobid "$1@$2" --destdir  "$DEBUG_DIR" >/dev/null 2>&1 \
     || jobsub_fetchlog -G "$EXPERIMENT" --jobid "$1@$2" --dest-dir "$DEBUG_DIR" >/dev/null 2>&1 \
     || log "$3: jobsub_fetchlog failed — run: jobsub_fetchlog -G $EXPERIMENT --jobid $1@$2 --destdir $DEBUG_DIR"
}

# Poll the queue until this cluster drains. Releases held jobs with
# bumped resources. A transient jobsub_q failure is retried, and drain is
# only declared after two consecutive empty reads (avoids a blip ending
# the wait early).
monitor_cluster(){ # cluster schedd label
    local cluster="$1" schedd="$2" label="$3"
    local q active held jid stable=0
    log "$label: monitoring cluster $cluster@$schedd (poll every ${POLL}s)"
    while true; do
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

# ---- pipeline phases -------------------------------------------------
submit_selection(){
    log "=== Stage 1: grid selection ==="
    log "Creating project $PROJ (batch-size $BATCH_SIZE)"
    python3 batch/medulla.py --experiment "$EXPERIMENT" --project-dir "$PROJ" \
        --create-project --toml "$SEL_TOML" --batch-size "$BATCH_SIZE" \
        --tag "$TAG" --gituser "$GITUSER" 2>&1 | tee -a "$LOGFILE" \
        || die "create-project failed"

    local db="$DEBUG_DIR/project.db"
    ifdh cp "$PROJ/project.db" "$db" >/dev/null 2>&1 \
        || die "cannot read $PROJ/project.db (refresh creds: kinit ; htgettoken)"
    NJOBS=$(sqlite3 "$db" "SELECT COUNT(*) FROM jobs;" 2>/dev/null)
    sqlite3 "$db" "SELECT cfg FROM configuration LIMIT 1;" 2>/dev/null \
        | grep -q selected_cc_Xg1p_stage3 \
        || die "stage3 tree missing from project config — is the local checkout on $TAG updated?"
    log "Project has $NJOBS jobs; stage2/3 config present."

    log "Launching $NJOBS jobs (mem=${MEMORY_MB}MB disk=${DISK_GB}GB lifetime=$LIFETIME)"
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
    monitor_cluster "$SEL_CLUSTER" "$SEL_SCHEDD" selection

    log "=== Stage 2: gather + hadd ==="
    local outs
    mapfile -t outs < <(ifdh ls "$PROJ/output" 2>/dev/null | grep -oE 'output_jobid[0-9]+\.root' | sort -u)
    local done=${#outs[@]}
    log "Selection produced $done / ${NJOBS:-?} job outputs."
    if [[ -n "${NJOBS:-}" && "$done" -lt "$NJOBS" ]]; then
        log "WARN: $((NJOBS - done)) selection job(s) missing output — fetching logs, continuing with partial stats."
        fetch_logs "$SEL_CLUSTER" "$SEL_SCHEDD" selection
    fi
    [[ "$done" -gt 0 ]] || die "no selection outputs produced — see logs in $DEBUG_DIR"

    local xrd=() f
    for f in "${outs[@]}"; do xrd+=("${XROOTD_DOOR}${PROJ}/output/${f}"); done
    log "hadd -> $SEL_HADD ($done files, via xrootd)"
    hadd -f "$SEL_HADD" "${xrd[@]}" 2>&1 | tail -5 | tee -a "$LOGFILE" \
        || die "hadd failed (adjust XROOTD_DOOR, or ifdh-cp the outputs local first)"
    [[ -s "$SEL_HADD" ]] || die "hadd produced no $SEL_HADD"
}

do_systematics(){
    log "=== Stage 3: systematics ==="
    ( cd build && ./systematics/run_systematics "../$SYS_TOML" ) 2>&1 | tail -15 | tee -a "$LOGFILE" \
        || die "run_systematics failed"
    [[ -s "$SYS_ROOT" ]] || die "run_systematics produced no $SYS_ROOT"
    local trees
    trees=$(root -l -b -q -e "TFile f(\"$SYS_ROOT\"); f.ls();" 2>/dev/null \
            | grep -oE 'selected_cc_Xg1p_stage[123]' | sort -u | tr '\n' ' ')
    log "sys ROOT CC trees: ${trees:-none}"
    echo "$trees" | grep -q selected_cc_Xg1p_stage3 \
        || die "stage3 tree not in $SYS_ROOT — selection/systematics did not produce it"
}

submit_plots(){
    log "=== Stage 4: stage sys ROOT + plot stage1/2/3 on grid ==="
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
    local dst="$PWD/CCSidebandPlots_$STAMP" cfg
    for cfg in "${CONFIGS[@]}"; do
        mkdir -p "$dst/$cfg"
        ifdh cp -D "$OUT/$cfg"/*.png "$dst/$cfg/" >/dev/null 2>&1 || true
        ifdh cp -D "$OUT/$cfg"/*.pdf "$dst/$cfg/" >/dev/null 2>&1 || true
    done
    log "Figures -> $dst ; debug logs (if any) -> $DEBUG_DIR ; full log -> $LOGFILE"
}

main(){
    log "CC sideband one-shot run starting (stamp $STAMP)"
    setup_env
    preflight
    submit_selection
    gather_and_merge
    do_systematics
    submit_plots
    retrieve
    log "DONE."
}
main "$@"
