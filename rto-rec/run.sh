#!/usr/bin/env bash
dt=8
dT=10
ncp=5
t_f=16
RUN_DIR="${t_f}_${dt}"

BASE_DIR="/Users/kevinnag/Documents/School/2026research/Williams-Otto-RTO/rto-rec"
export SV_DIR="$RUN_DIR" dt ncp dT t_f
SRC_CSV="${BASE_DIR}/data/wo.csv"
DST_DIR="${BASE_DIR}/sims/${RUN_DIR}"
python3 "${BASE_DIR}/fa.py"
mkdir -p "$DST_DIR"
mkdir -p "$DST_DIR/logs"
cp "$SRC_CSV" "$DST_DIR/wo.csv"

for i in $(seq 1 $(((t_f+dt-1)/dt))); do
   python3 "${BASE_DIR}/wo-rec.py"
   python3 "${BASE_DIR}/rto.py"
done



