#!/usr/bin/env bash
set -euo pipefail

RUN_DIR="rec_5_fwd"
N_RUNS=180
export SV_DIR="$RUN_DIR"

SRC_CSV="/Users/kevinnag/Documents/School/2026research/Williams-Otto-RTO/williams-otto-rec/wo.csv"
DST_DIR="/Users/kevinnag/Documents/School/2026research/Williams-Otto-RTO/williams-otto-rec/sims/${RUN_DIR}"

mkdir -p "$DST_DIR"
mkdir -p "$DST_DIR/logs"
cp "$SRC_CSV" "$DST_DIR/wo.csv"


for i in $(seq 1 "$N_RUNS"); do
    python3 /Users/kevinnag/Documents/School/2026research/Williams-Otto-RTO/williams-otto-rec/wo-rec.py
done
