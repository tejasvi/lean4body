#!/bin/bash
# The runs of Table 2 at the threshold 19/32 = 0.59375 (the bound 3/4 follows a fortiori), then
# the end-to-end sampling test of their accepted leaves.  verify.py appends to its log, so the
# script removes the old logs first.
cd "$(dirname "$0")"
PY=${PY:-python}
rm -f direct_run_19_32.log chart_run_19_32.log
$PY verify.py direct --thresh 19/32 --budget 300000 --workers 16 --log direct_run_19_32.log --json direct_run_19_32.json --leaves direct_leaves_19_32.npy
$PY verify.py chart  --thresh 19/32 --budget 300000 --workers 16 --log chart_run_19_32.log  --json chart_run_19_32.json --leaves chart_leaves_19_32.npy
$PY check_leaves.py direct direct_leaves_19_32.npy 4000 3000 --workers 16 --seed 777 > check_leaves_direct_19_32.log 2>&1
$PY check_leaves.py chart  chart_leaves_19_32.npy  4000 3000 --workers 16 --seed 777 > check_leaves_chart_19_32.log 2>&1
