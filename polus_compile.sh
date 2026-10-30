#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
CC=${CC:-gcc}
command -v "$CC" >/dev/null
"$CC" --version
"$CC" -O3 -std=c11 -Wall -Wextra -Wpedantic solver_seq.c -o solver_seq -lm
"$CC" -O3 -std=c11 -Wall -Wextra -Wpedantic -fopenmp solver_omp.c -o solver_omp -lm
