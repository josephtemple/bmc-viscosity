#!/bin/bash
#
# Politely launch a LAMMPS run on a shared, scheduler-less cluster node,
# inside a detachable tmux session. Output streams live to the attached
# pane and is also mirrored to stdout.log for later reference.
#
# Usage: bash run.sh <input.in> <run_name>

set -euo pipefail

INPUT_FILE="${1:?Usage: bash run.sh <input.in> <run_name>}"
RUN_NAME="${2:?Usage: bash run.sh <input.in> <run_name>}"
OUTDIR="/data/jtemple/${RUN_NAME}"
LOGFILE="${OUTDIR}/lammps_run.log"
SESSION_NAME="lammps_${RUN_NAME}"

mkdir -p "${OUTDIR}"

# ---------------------------------------------------------------
# 1. See how many cores exist, and how many are already in use.
# ---------------------------------------------------------------
TOTAL_CORES=$(nproc)
LOAD_1MIN=$(awk '{print $1}' /proc/loadavg)
BUSY_CORES=$(printf '%.0f' "${LOAD_1MIN}")
FREE_CORES=$(( TOTAL_CORES - BUSY_CORES ))

echo "Total cores on this node: ${TOTAL_CORES}"
echo "1-min load average:       ${LOAD_1MIN}"
echo "Estimated free cores:     ${FREE_CORES}"

# ---------------------------------------------------------------
# 2. Decide how many ranks to actually use.
# ---------------------------------------------------------------
MAX_FRACTION=0.5
CUSHION=2

SAFE_CORES=$(( FREE_CORES - CUSHION ))
CAP_CORES=$(printf '%.0f' "$(echo "${TOTAL_CORES} * ${MAX_FRACTION}" | bc)")

if [ "${SAFE_CORES}" -lt 1 ]; then
    echo "Not enough free cores right now (estimated ${FREE_CORES} free)."
    echo "Check 'top' or 'htop' manually before running -- refusing to launch."
    exit 1
fi

NPROCS=$(( SAFE_CORES < CAP_CORES ? SAFE_CORES : CAP_CORES ))

echo "Launching with ${NPROCS} MPI ranks (out of ${TOTAL_CORES} total cores)."

# ---------------------------------------------------------------
# 3. Check for a name collision before launching.
# ---------------------------------------------------------------
if tmux has-session -t "${SESSION_NAME}" 2>/dev/null; then
    echo "A tmux session named '${SESSION_NAME}' already exists."
    echo "Attach with: tmux attach -t ${SESSION_NAME}"
    echo "Or pick a different run_name."
    exit 1
fi

# ---------------------------------------------------------------
# 4. Launch inside tmux, at low CPU + I/O priority. Output streams
#    live to the pane and is mirrored to stdout.log via tee.
# ---------------------------------------------------------------
tmux new-session -d -s "${SESSION_NAME}" \
    "nice -n 15 ionice -c2 -n7 mpirun --mca mtl ^ofi -np ${NPROCS} lmp -in ${INPUT_FILE} -log ${LOGFILE} \
     2>&1 | tee ${OUTDIR}/stdout.log"

echo "Launched inside tmux session: ${SESSION_NAME}"
echo "Log file:    ${LOGFILE}"
echo "Stdout:      ${OUTDIR}/stdout.log (also visible live via tmux attach)"
echo ""
echo "To watch it live:     tmux attach -t ${SESSION_NAME}"
echo "  (detach again with Ctrl+b then d -- this does NOT stop the run)"
echo "To check it's alive:  tmux has-session -t ${SESSION_NAME} && echo running"
echo "To stop it cleanly:   tmux send-keys -t ${SESSION_NAME} C-c"
echo "To kill the session:  tmux kill-session -t ${SESSION_NAME}"