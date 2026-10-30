#!/usr/bin/env bash
# Run INSIDE one LSF allocation so all configurations use the same node/resources.
set -euo pipefail
cd "$(dirname "$0")"
: "${LSB_JOBID:?Submit this script with bsub; do not run on the login node}"
suite=${1:-check}
runs=${2:-10}
[[ "$runs" =~ ^[1-9][0-9]*$ ]]
export POLUS_BIND=legacy
export OMP_DYNAMIC=FALSE
tag="v2_${suite}_${LSB_JOBID}"
completed=0
run_config() {
    bash ./polus_job.sh "$1" "$2" "$3" "$4" "$runs" "$tag"
    completed=$((completed+runs))
    printf 'PROGRESS %s/%s completed\n' "$completed" "$total"
}
case "$suite" in
    check)
        [[ "$runs" == 1 ]] || { echo 'Use check 1 for the preliminary run' >&2; exit 1; }
        total=6
        export SAVE_SOLUTIONS=0
        run_config seq 400 600 1
        run_config omp 400 600 1
        run_config omp 400 600 4
        run_config omp 400 600 16
        run_config seq 800 1200 1
        run_config omp 800 1200 32
        ;;
    full)
        total=$((16*runs))
        export SAVE_SOLUTIONS=1
        for n in 20 40 80; do run_config seq "$n" "$n" 1; done
        for p in 1 4 16; do run_config omp 80 80 "$p"; done
        run_config seq 400 600 1
        for p in 2 4 8 16; do run_config omp 400 600 "$p"; done
        run_config seq 800 1200 1
        for p in 4 8 16 32; do run_config omp 800 1200 "$p"; done
        ;;
    *) echo 'Usage: polus_suite_v2.sh check 1 | full 10' >&2; exit 1;;
esac
printf 'COMPLETE %s/%s; results/%s\n' "$completed" "$total" "$tag"
