#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mode=${1:?mode}; m=${2:?M}; n=${3:?N}; threads=${4:?threads}
runs=${5:-10}; tag=${6:-manual}
[[ "$mode" == seq || "$mode" == omp ]]
for value in "$m" "$n" "$threads" "$runs"; do
    [[ "$value" =~ ^[1-9][0-9]*$ ]]
done
[[ "$tag" =~ ^[A-Za-z0-9_.-]+$ ]]
export OMP_NUM_THREADS="$threads"
export OMP_DYNAMIC=FALSE
if [[ ${POLUS_BIND:-modern} == legacy ]]; then
    export OMP_PROC_BIND=TRUE
    unset OMP_PLACES
else
    export OMP_PROC_BIND=spread
    export OMP_PLACES=threads
fi
export LC_ALL=C
root="results/${tag}/${mode}_${m}_${n}_${threads}_${LSB_JOBID:-local}"
mkdir -p "$root"
exec > >(tee "$root/job.out") 2> >(tee "$root/job.err" >&2)
printf 'host=%s job=%s mode=%s grid=%sx%s threads=%s runs=%s\n' \
    "$(hostname)" "${LSB_JOBID:-local}" "$mode" "$m" "$n" "$threads" "$runs"
printf 'OMP_NUM_THREADS=%s OMP_PROC_BIND=%s OMP_PLACES=%s\n' \
    "$OMP_NUM_THREADS" "$OMP_PROC_BIND" "${OMP_PLACES:-unset}"
if command -v lscpu >/dev/null; then lscpu; fi
if command -v taskset >/dev/null; then taskset -pc $$ || true; fi
if command -v sha256sum >/dev/null; then sha256sum solver_seq.c solver_omp.c solver_seq solver_omp; fi
printf 'M,N,threads,iterations,time,step_norm,residual_norm\n' > "$root/runs.csv"
for ((run=1;run<=runs;run++)); do
    # One solution per configuration; all later runs skip CSV output.
    output=-
    if ((run==1)) && [[ ${SAVE_SOLUTIONS:-1} == 1 ]]; then output="$root/solution.csv"; fi
    if [[ "$mode" == seq ]]; then
        ./solver_seq "$m" "$n" 1e-8 "$output" > "$root/run_${run}.out" 2> "$root/run_${run}.err"
    else
        ./solver_omp "$m" "$n" "$threads" 1e-8 "$output" > "$root/run_${run}.out" 2> "$root/run_${run}.err"
        awk -v expected="$threads" '$1=="Threads" && $3==expected {ok=1} END {exit !ok}' "$root/run_${run}.out"
    fi
    grep -q '^Converged     : yes$' "$root/run_${run}.out"
    grep -q '^Timing scope  : iterative_solver_only_v2$' "$root/run_${run}.out"
    awk -F, '$1=="RESULT" {sub(/^RESULT,/, "");print;found++} END {exit found!=1}' \
        "$root/run_${run}.out" >> "$root/runs.csv"
    printf 'Completed %s %sx%s p=%s: %s/%s; ' "$mode" "$m" "$n" "$threads" "$run" "$runs"
    grep '^Algorithm time:' "$root/run_${run}.out"
done
awk -F, 'NR==1 {print;next} !seen || $5<best {best=$5;line=$0;seen=1} END {if(seen) print line;else exit 1}' \
    "$root/runs.csv" > "$root/best.csv"
cat "$root/best.csv"
