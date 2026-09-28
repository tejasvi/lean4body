"""Negative control (a test, NOT part of the proof): run the same branch-and-bound with a threshold
BELOW the numerically known maximum of min_y trS on Z (~0.59) on the initial cells that contain
the numerical maximiser.  The claim is then false, so the verifier must report unproved leaves.
For some failing leaves we print min_y trS at the root above the leaf centre (floats), to show
that they sit where trS really exceeds the threshold.

usage: negative_control.py   (writes to stdout)
"""
import sys
import time

import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import verify
from formulas import prims_direct, prims_chart, certificate, trS, optimal_Y, Ftil, region_scaled


def root_minY(mode, B):
    """float: root of F above the centre of B (bisection on the linear x-range), min_y trS there"""
    p = [0.5 * (lo + hi) for lo, hi in B]
    prims = prims_direct if mode == 'direct' else prims_chart
    f = lambda t: Ftil(prims(*p, t))
    if mode == 'direct':
        a, b, c = p
        lo = max((a * a - 2 * b - b * b) / (2 * a), (1 - 2 * a * c - c * c) / (2 * a), (c - a) / (2 * b),
                 -(1 - b) / (2 * c) if c > 0 else -1.0, -1.0)
        hi = min((a - c) / 2, (1 - b) / (2 * a), 1.0)
    else:
        b, al, ga = p
        u = 1 + b * al
        lo = max(al / 2 - (2 + b) / (2 * u), -ga - (b * ga * ga + al) / (2 * u))
        hi = min((al - ga) / 2, -(1 + al) / (2 * u))
    if not lo < hi or f(lo) * f(hi) > 0:
        return None
    for _ in range(100):
        m = 0.5 * (lo + hi)
        if f(lo) * f(m) <= 0:
            hi = m
        else:
            lo = m
    t = 0.5 * (lo + hi)
    P = prims(*p, t)
    if min(region_scaled(P)) < 0:
        return None
    cert = certificate(P)
    return t, trS(cert, optimal_Y(cert))


def run(mode, cell, thresh, budget):
    boxes = dict(verify.initial_boxes(mode))
    t0 = time.time()
    st = verify.process_initial((mode, cell, boxes[cell], thresh, budget, False))
    print(f"{mode} cell {cell} = {boxes[cell]}, threshold {thresh}, node budget {budget}: "
          f"{st['nodes']} nodes, accepted {st['accepted']}, discarded {st['disc_region'] + st['disc_F']}, "
          f"UNPROVED {len(st['fail'])}  ({time.time() - t0:.1f}s)")
    vals = []
    for f in st['fail']:
        r = root_minY(mode, f['box'])
        if r is not None:
            vals.append((r[1], f['box'], r[0]))
    vals.sort(key=lambda v: -v[0])
    if vals:
        above = sum(1 for v in vals if v[0] > float(eval(thresh)))
        print(f"  min_y trS at the root above the centre of the unproved leaves: max {vals[0][0]:.5f}, "
              f"min {vals[-1][0]:.5f}; {above}/{len(vals)} exceed the threshold")
        print(f"  e.g. leaf {vals[0][1]}, dependent coord {vals[0][2]:.6f}, min_y trS {vals[0][0]:.5f}")


if __name__ == '__main__':
    run('direct', (2, 3, 8), '9/16', 4000)
    run('chart', (7, 11, 8), '9/16', 4000)
    run('direct', (2, 3, 8), '37/64', 4000)
