#!/usr/bin/env bash
set -euo pipefail

dt=8
dT=10
ncp=5
t_f=900
RUN_DIR="${t_f}_${dt}c"

export SV_DIR="$RUN_DIR" dt ncp dT t_f
BASE_DIR="/Users/kevinnag/Documents/School/2026research/Williams-Otto-RTO/rto-rec"
#SRC_CSV="${BASE_DIR}/wo.csv"
#DST_DIR="${BASE_DIR}/sims/${RUN_DIR}"

#mkdir -p "$DST_DIR"
#mkdir -p "$DST_DIR/logs"
#cp "$SRC_CSV" "$DST_DIR/wo.csv"

#for i in $(seq 1 $(((t_f+dt-1)/dt))); do
#    python3 "${BASE_DIR}/wo-rec.py"
#done

python3 "${BASE_DIR}/fa.py"
python3 "${BASE_DIR}/test.py"

