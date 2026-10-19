#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=24
#SBATCH --account=def-tobi
#SBATCH -t 23:55:01
# Ask SLURM to send SIGUSR1 to this batch shell 15 min before the wall-time
# kill, so the next segment can be queued while this one is still running.
#SBATCH --signal=B:USR1@900

# 8 ranks x 24 OMP threads/rank (= 192 cores, one rank per NUMA domain on
# Trillium's 8x24 topology per runModular_trillium_calibrate.sh's numactl
# output), NOT a single big-thread-count rank.
#
# Why: measured locally (2-socket/32-physical-core dev box, same ITER case)
# that modular's own per-rank OMP load balance is fine up to its physical
# core count but gets volatile and erratic past it (mover_push max/avg rose
# from ~1.0 at 8-32 threads/rank to spikes of 1.67 at 64 - three thread
# counts past this box's 32 physical cores), and separately that modular's
# per-thread scaling efficiency falls behind legacy's well before that even
# at "safe" thread counts (legacy: ~97% scaling efficiency 8->64 threads;
# modular: ~50%, plateauing hard past 32). Legacy showed neither effect
# (flat ~1.04-1.10 imbalance at every thread count tested, up to 64).
# Splitting the node into more, smaller-thread-count ranks avoids both:
# every rank's own OMP team stays comfortably inside the region where
# modular scales well, and MPI (not OpenMP) does the work of filling the
# rest of the node.
#
# This 8x24 split is carried over from the calibration sweep's NUMA-aligned
# candidate, NOT independently re-verified at ITER scale with multiple MPI
# ranks (today's local testing only exercised single-rank/many-OMP-thread
# configs) - rerun runModular_trillium_calibrate.sh on Trillium once
# available to confirm/tune this against the other candidate splits before
# treating it as final.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=true
# OMP_PLACES=threads, not cores: an overnight affinity sweep (OMP=48,
# oversubscribed on the dev box) measured threads giving the best imbalance
# ratio (1.091) vs cores (1.113) and sockets (1.84-2.06, much worse - avoid
# entirely). Difference is modest at cores vs threads but free to take.
export OMP_PLACES=threads
export I_MPI_PIN_DOMAIN=omp

echo "OMP_NUM_THREADS=$OMP_NUM_THREADS, MPI ranks=$SLURM_NTASKS"

# ---------------------------------------------------------------------------
# Auto-resubmission chain
#
# Every segment (wall-time kill, crash, or normal end of nsteps) queues the
# next one with --dependency=afterany, which restarts from the newest
# *complete* backup. Submit the first segment normally (sbatch runModular.sh).
# Its restart flag in conditions.inp is used as-is. Later segments force
# "yes".
#
# To stop the chain: `touch STOP_CHAIN` in the run directory (checked both
# when queuing and when a segment starts), or scancel the pending job.
# The chain also stops by itself after MAX_SEGMENTS segments, when no valid
# backup exists, or after MAX_NO_PROGRESS consecutive segments that restart
# from the same backup (= crashing before writing a new one).
#
# Caveats from the code (mod_restart.f90 / mod_simulation.f90):
#   - flag_restart=1 reads only Backups/DATA.BAK/, but backups alternate
#     BAK/BAK2, so the newest complete one is copied into BOTH dirs before
#     restarting (a crash during the restarted run's first write to BAK
#     then still leaves BAK2 intact).
#   - A restarted run executes nsteps MORE steps from the backup's step; it
#     does not stop at tmax. MAX_SEGMENTS is the run-length control here.
#   - np_dup is re-applied on every restart, so it is forced to 1. here.
#   - Steps between the last backup and the end of a segment are redone,
#     so Output/*.dat time series get those steps appended twice.
# ---------------------------------------------------------------------------
cd "${SLURM_SUBMIT_DIR:-.}" || exit 1

SEGMENT=${SEGMENT:-1}
MAX_SEGMENTS=${MAX_SEGMENTS:-30}
MAX_NO_PROGRESS=${MAX_NO_PROGRESS:-2}
CHAIN_SCRIPT=${CHAIN_SCRIPT:-runModular.sh}
STOP_FILE=STOP_CHAIN
STATE_FILE=.chain_state
COND=input_dir/conditions.inp
BAK=Backups/DATA.BAK
BAK2=Backups/DATA.BAK2
NRANKS=$SLURM_NTASKS

echo "Chain segment $SEGMENT/$MAX_SEGMENTS, job $SLURM_JOB_ID"

# Rewrite the "Restart (yes/no/...)" line of conditions.inp:
# set_restart_line <flag> <np_dup>  (OMP/MPI fields always follow this job)
set_restart_line() {
  sed -i -E "s|^[[:space:]]*[^[:space:]]+[[:space:]]+[^[:space:]]+[[:space:]]+[^[:space:]]+[[:space:]]+[^[:space:]]+([[:space:]]*! Restart \()|$1 $OMP_NUM_THREADS $NRANKS $2\1|" "$COND"
  grep "! Restart (" "$COND"
}

# True if dir holds a nonempty particles<r>.bak for every rank
backup_ok() {
  local r
  for ((r = 0; r < NRANKS; r++)); do
    [[ -s $1/particles$r.bak ]] || return 1
  done
}

# Newest mtime (epoch s) over the dir's rank files
backup_mtime() {
  stat -c %Y "$1"/particles*.bak | sort -n | tail -1
}

if [[ -e $STOP_FILE ]]; then
  echo "$STOP_FILE present - not running."
  exit 0
fi

if (( SEGMENT == 1 )); then
  # Keep the user's restart flag and np_dup and set OMP/MPI to match the job
  read -r flag _ _ dup _ < <(grep "! Restart (" "$COND")
  set_restart_line "$flag" "$dup"
  echo "0 0" > "$STATE_FILE"
else
  # Candidate backups, newest first
  mapfile -t cands < <(
    for d in "$BAK" "$BAK2"; do
      backup_ok "$d" && echo "$(backup_mtime "$d") $d"
    done | sort -rn | awk '{print $2}')

  # If the previous segment died mid-backup (last "Backing up" with no
  # "Done!" after it), the newest dir is the partial one: skip it.
  if [[ -n $PREV_DUMP && -f $PREV_DUMP ]] && \
     grep -E "Backing up simulation data|^[[:space:]]*Done!" "$PREV_DUMP" | tail -1 | grep -q "Backing up"; then
    echo "Previous segment was killed during a backup write; skipping ${cands[0]}"
    cands=("${cands[@]:1}")
  fi

  if (( ${#cands[@]} == 0 )); then
    echo "No complete backup found in $BAK or $BAK2 - stopping chain."
    exit 1
  fi
  src=${cands[0]}
  sig=$(backup_mtime "$src")

  read -r last_sig no_prog < "$STATE_FILE" 2>/dev/null || { last_sig=0; no_prog=0; }
  if [[ $sig == "$last_sig" ]]; then
    no_prog=$((no_prog + 1))
    echo "No new backup since last restart ($no_prog/$MAX_NO_PROGRESS)."
    if (( no_prog >= MAX_NO_PROGRESS )); then
      echo "Run keeps failing before writing a backup - stopping chain."
      exit 1
    fi
  else
    no_prog=0
  fi
  echo "$sig $no_prog" > "$STATE_FILE"

  # Make both backup dirs hold the chosen backup (copy via tmp + mv so a
  # kill mid-copy never leaves a truncated file under the real name)
  echo "Restarting from $src (mtime $(date -d @"$sig"))"
  for dst in "$BAK" "$BAK2"; do
    [[ $dst == "$src" ]] && continue
    for ((r = 0; r < NRANKS; r++)); do
      cp -p "$src/particles$r.bak" "$dst/particles$r.bak.tmp" && \
        mv "$dst/particles$r.bak.tmp" "$dst/particles$r.bak" || exit 1
    done
  done

  set_restart_line yes 1.
fi

RESUBMITTED=0
resubmit() {
  (( RESUBMITTED )) && return
  RESUBMITTED=1
  if [[ -e $STOP_FILE ]]; then
    echo "$STOP_FILE present - not resubmitting."
    return
  fi
  if (( SEGMENT >= MAX_SEGMENTS )); then
    echo "Reached MAX_SEGMENTS=$MAX_SEGMENTS - not resubmitting."
    return
  fi
  sbatch --dependency=afterany:"$SLURM_JOB_ID" \
    --export=ALL,SEGMENT=$((SEGMENT + 1)),MAX_SEGMENTS=$MAX_SEGMENTS,MAX_NO_PROGRESS=$MAX_NO_PROGRESS,CHAIN_SCRIPT=$CHAIN_SCRIPT,PREV_DUMP=$DUMP \
    "$CHAIN_SCRIPT"
}
trap 'echo "Wall time approaching - queuing next segment"; resubmit' USR1

# One log per segment; dump.modular always points at the current one
DUMP=dump.modular.$SLURM_JOB_ID
ln -sfn "$DUMP" dump.modular

# Background + wait so the USR1 trap can fire while mpirun runs
mpirun -np "$SLURM_NTASKS" ./build/run_min > "$DUMP" 2>&1 &
MPI_PID=$!
while true; do
  wait "$MPI_PID"
  rc=$?
  kill -0 "$MPI_PID" 2>/dev/null || break
done
echo "run_min exited with status $rc"

resubmit
exit "$rc"
