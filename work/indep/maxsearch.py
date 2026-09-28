"""Float search for the maximum of min_y trS over Z (NOT part of the proof; it only tells how
tight the certified bounds are).  x (xi) is the root of F in the region above (a,b,c)
((b,al,ga)), found by bisection; min over y of the quadratic trS is explicit.  Local search:
compass/pattern search with shrinking steps, started from a coarse grid's best points and from
given points, clamped to the parameter domain.

usage: maxsearch.py direct|chart
"""
import itertools
import sys

import numpy as np

import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from formulas import prims_direct, prims_chart, certificate, trS, optimal_Y, Ftil, region_scaled

DOM = {'direct': [(0.5, 1.75), (0.125, 1.0), (0.0, 1.75)],
       'chart': [(1e-12, 0.125), (-2.0, 1.0), (0.0, 3.0)]}


def xrange_(mode, p):
    if mode == 'direct':
        a, b, c = p
        lows = [(a * a - 2 * b - b * b) / (2 * a), (1 - 2 * a * c - c * c) / (2 * a), (c - a) / (2 * b), -1.0]
        if c > 0:
            lows.append(-(1 - b) / (2 * c))
        highs = [(a - c) / 2, (1 - b) / (2 * a), 1.0]
    else:
        b, al, ga = p
        u = 1 + b * al
        lows = [al / 2 - (2 + b) / (2 * u), -ga - (b * ga * ga + al) / (2 * u)]
        highs = [(al - ga) / 2, -(1 + al) / (2 * u)]
    return max(lows), min(highs)


def value(mode, p):
    """min_y trS at the Z point above p, or None"""
    prims = prims_direct if mode == 'direct' else prims_chart
    lo, hi = xrange_(mode, p)
    if not lo < hi:
        return None
    f = lambda t: Ftil(prims(*p, t))
    try:
        flo, fhi = f(lo), f(hi)
    except (ValueError, ZeroDivisionError):
        return None
    if flo * fhi > 0:
        return None
    for _ in range(80):
        m = 0.5 * (lo + hi)
        fm = f(m)
        if flo * fm <= 0:
            hi = m
        else:
            lo, flo = m, fm
    t = 0.5 * (lo + hi)
    P = prims(*p, t)
    if min(region_scaled(P)) < -1e-12:
        return None
    cert = certificate(P)
    return trS(cert, optimal_Y(cert)), t


def clamp(mode, p):
    return tuple(min(max(v, lo), hi) for v, (lo, hi) in zip(p, DOM[mode]))


def pattern_search(mode, p, step=0.02, tol=1e-9):
    best = value(mode, p)
    if best is None:
        return None
    fbest = best[0]
    dirs = [d for d in itertools.product((-1, 0, 1), repeat=3) if any(d)]
    while step > tol:
        improved = False
        for d in dirs:
            q = clamp(mode, tuple(v + step * di for v, di in zip(p, d)))
            r = value(mode, q)
            if r is not None and r[0] > fbest:
                p, fbest, best = q, r[0], r
                improved = True
                break
        if not improved:
            step /= 2
    return fbest, p, best[1]


def main(mode):
    dom = DOM[mode]
    n = 28
    grids = [np.linspace(lo, hi, n) for lo, hi in dom]
    vals = []
    for p in itertools.product(*grids):
        r = value(mode, p)
        if r is not None:
            vals.append((r[0], p))
    vals.sort(key=lambda v: -v[0])
    print(f"{mode}: coarse grid {n}^3, {len(vals)} Z points; best {vals[0][0]:.6f} at {vals[0][1]}")
    starts = [v[1] for v in vals[:25]]
    extra = {'direct': [(1.695482, 0.957925, 0.558907), (0.64831, 0.31780, 0.50495)],
             'chart': [(0.122146, -0.587611, 1.078179), (0.125, -0.576, 1.068)]}[mode]
    starts += extra
    res = []
    for s in starts:
        r = pattern_search(mode, s)
        if r is not None:
            res.append(r)
    res.sort(key=lambda r: -r[0])
    seen = []
    for f, p, t in res:
        if any(max(abs(u - v) for u, v in zip(p, q)) < 1e-3 for q in seen):
            continue
        seen.append(p)
        print(f"  local max {f:.7f} at params {tuple(round(v, 6) for v in p)}, dependent {t:.6f}")
        if len(seen) >= 6:
            break


if __name__ == '__main__':
    main(sys.argv[1])
