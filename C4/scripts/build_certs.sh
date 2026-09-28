#!/bin/bash
# Build the project so that the kernel's memory stays well below 128 GB.  Lake checks up to one
# module per hardware thread at once, and a module needs up to ~6 GB while the kernel evaluates
# its largest certificate.  So the 52 modules whose largest certificate takes more than ~8 s are
# built first, in the batches of `heavy_batches.txt` (written by `mkcert.py`; at most ~40 GB each,
# by the estimate 0.4 GB + 0.125 GB per second of kernel time), then the other certificate
# modules, and then everything else.  After the certificates and at the end, the script prints
# the elapsed time and, with `times`, the processor time (user and system) of the script itself
# and of the processes it started.
set -e
cd "$(dirname "$0")/.."
while read -r mods; do
  targets=""
  for m in $mods; do targets="$targets C4Cert.$m"; done
  lake build $targets
done < scripts/heavy_batches.txt
lake build C4Cert
echo "certificates built after $SECONDS s; processor time of the script, then of its children:"
times
lake build
echo "everything built after $SECONDS s; processor time of the script, then of its children:"
times
