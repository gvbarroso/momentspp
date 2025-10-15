#!/usr/bin/env bash
[[ -f ~/.bashrc ]] && source ~/.bashrc
set -euo pipefail

# Generate job list: one line per (order, model, P)
generate_jobs() {
  local n=$(tail -n +2 params.csv | wc -l)
  for (( o=50; o<=1000; o+=50 )); do # parallelize over Order first to optimize time
    for (( i=1; i<=n; i++ )); do
      echo "$o $i 16" # double (16-digit) precision
      echo "$o $i 32" # mpreal (32-digit) precision
    done
  done
}

# Worker function: runs a single momentspp job
run_job() {
  local o="$1"
  local m="$2"
  local p="$3"
  (
    cd "model_${m}"
    momentspp param=../opt.bpp F="model_${m}.yaml" O="$o" V=1 P="$p" NT=1 &> "log_O_${o}_P${p}.txt"
    echo "done: Order $o, model $m, P=$p"
  )
}

export -f run_job

# Launch jobs with Order-first priority and max 16 parallel jobs
generate_jobs | sort -k1,1n -k2,2n | xargs -n 3 -P 16 bash -c 'run_job "$@"' _

