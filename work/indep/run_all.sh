#!/bin/bash
# The runs of Table 2 at the threshold 3/4: the verifier on the direct region and on the chart,
# 16 workers each, then the end-to-end sampling test of the accepted leaves.  verify.py appends
# to its log, so the script removes the old logs first.
cd "$(dirname "$0")"
PY=${PY:-python}
rm -f direct_run.log chart_run.log
$PY verify.py direct --workers 16 --log direct_run.log --json direct_run.json --leaves direct_leaves.npy
$PY verify.py chart  --workers 16 --log chart_run.log  --json chart_run.json --leaves chart_leaves.npy
$PY check_leaves.py direct direct_leaves.npy 4000 2000 --workers 16 > check_leaves_direct.log 2>&1
$PY check_leaves.py chart  chart_leaves.npy  4000 2000 --workers 16 > check_leaves_chart.log 2>&1
