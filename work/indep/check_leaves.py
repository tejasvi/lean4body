"""End-to-end test of a finished verifier run (a test, NOT part of the proof).

Samples points of Z using only the literal SPEC formulas (formulas.literal: positions -> R, r, s,
areas, lam, w, D, tau, beta, F, region) in mpmath, and checks for each sampled point that
  (1) its parameters lie in some accepted leaf of the run (a point of Z must never fall outside
      the accepted leaves: everything else was discarded as region-empty),
  (2) its dependent coordinate x (chart: xi) lies in that leaf's x-piece, and
      eps = x - x_lin(p) lies in the leaf's sheared range E,
  (3) trS(y) at the point, with the y stored for that leaf (chart: y = Y/b), is <= the leaf's
      certified upper bound.
Points: uniformly random parameters, and random parameters inside the leaves with the largest
certified bounds (the tightest leaves).

usage: check_leaves.py direct|chart LEAVES.npy N_UNIFORM N_TOP [--workers W] [--seed S]
"""
import argparse
import random
import sys
import time
from multiprocessing import Pool

import mpmath
import numpy as np
from mpmath import mp, mpf

import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from formulas import literal, trS_literal

DOMAINS = {'direct': [(0.5, 1.75, 20), (0.125, 1.0, 14), (0.0, 1.75, 28)],
           'chart': [(0.0, 0.125, 8), (-2.0, 1.0, 24), (0.0, 3.0, 24)]}
LEAVES = None
CELLS = None
MODE = None


def to_direct(mode, p, t):
    """(p1, p2, p3, dependent) -> (a, b, c, x)"""
    if mode == 'direct':
        return p[0], p[1], p[2], t
    b, al, ga = p
    return 1 + b * al, b, b * ga, mpf(1) / 2 + b * t


def dep_scan_range(mode, p):
    """a generous scan range for the dependent coordinate, NOT derived from the verifier's
    region bounds: x in (-1, 1) (direct); xi in [-6, 3] intersected with x in (-1, 1) (chart)"""
    if mode == 'direct':
        return -1 + 1e-9, 1 - 1e-9
    b = float(p[0])
    return max(-6.0, (-1.5 + 1e-9) / b), min(3.0, (0.5 - 1e-9) / b)


def Fscaled(mode, p, t):
    a, b, c, x = to_direct(mode, p, t)
    F = literal(a, b, c, x)['F']
    return F if mode == 'direct' else F / (b * b)


def find_roots(mode, p, n=600):
    """roots of F (chart: F/b^2) found by a float sign scan of the literal F over the scan range,
    refined at 50 digits; kept if all six region inequalities hold (literal formulas)"""
    lo, hi = dep_scan_range(mode, p)
    pf = tuple(float(v) for v in p)
    ts = [lo + (hi - lo) * k / n for k in range(n + 1)]
    vals = []
    for t in ts:
        try:
            vals.append(float(Fscaled(mode, pf, t)))
        except (ZeroDivisionError, ValueError):
            vals.append(float('nan'))
    out = []
    mp.dps = 50
    for k in range(n):
        f0, f1 = vals[k], vals[k + 1]
        if not (f0 * f1 < 0):
            continue
        try:
            cand = mpmath.findroot(lambda t: Fscaled(mode, p, t), (mpf(ts[k]), mpf(ts[k + 1])),
                                   solver='anderson')
        except (ZeroDivisionError, ValueError):
            continue
        if not (mpf(ts[k]) <= cand <= mpf(ts[k + 1])):
            continue
        L = literal(*to_direct(mode, p, cand))
        if min(L['region']) >= -mpf(10) ** -40:
            out.append((cand, L))
    return out


def init(mode, leaves):
    global LEAVES, CELLS, MODE
    MODE = mode
    LEAVES = leaves
    CELLS = {}
    for r, row in enumerate(leaves):
        CELLS.setdefault((int(row[0]), int(row[1]), int(row[2])), []).append(r)


def cell_of(mode, p):
    idx = []
    for v, (lo, hi, n) in zip(p, DOMAINS[mode]):
        w = (hi - lo) / n
        k = int((float(v) - lo) // w)
        idx.append(min(max(k, 0), n - 1))
    return tuple(idx)


def check_point(p):
    """returns list of result dicts for the Z points above parameter point p"""
    mode = MODE
    res = []
    for t, L in find_roots(mode, p):
        rows = [r for r in CELLS.get(cell_of(mode, p), [])
                if all(LEAVES[r][3 + 2 * k] <= p[k] <= LEAVES[r][4 + 2 * k] for k in range(3))]
        if not rows:
            res.append(dict(p=[float(v) for v in p], t=float(t), status='NOT COVERED'))
            continue
        best = None
        for r in rows:
            row = LEAVES[r]
            Xlo, Xhi, xc, ph1, ph2, ph3, Elo, Ehi, Y1, Y2, ub = row[9:20]
            mids = [(row[3 + 2 * k] + row[4 + 2 * k]) / 2 for k in range(3)]
            eps = t - (xc + sum(mpf(ph) * (p[k] - mids[k]) for k, ph in enumerate((ph1, ph2, ph3))))
            inX = Xlo - 1e-13 <= t <= Xhi + 1e-13
            inE = Elo - 1e-13 <= eps <= Ehi + 1e-13
            b = p[0] if mode == 'chart' else 1
            y = (mpf(Y1) / b, mpf(Y2) / b)
            tr = trS_literal(L, y)
            d = dict(p=[float(v) for v in p], t=float(t), leaf=int(r), inX=bool(inX), inE=bool(inE),
                     trS=float(tr), ub=float(ub), margin=float(ub - tr))
            d['status'] = 'ok' if (inX and inE and tr <= ub + 1e-12) else 'BAD'
            if best is None or (d['status'] == 'ok' and best['status'] != 'ok'):
                best = d
        res.append(best)
    return res


def rand_in(lo, hi, rng):
    return mpf(lo) + (mpf(hi) - mpf(lo)) * mpf(rng.random())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['direct', 'chart'])
    ap.add_argument('leaves')
    ap.add_argument('n_uniform', type=int)
    ap.add_argument('n_top', type=int)
    ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--seed', type=int, default=12345)
    args = ap.parse_args()
    assert args.workers <= 16
    mp.dps = 50
    leaves = np.load(args.leaves)
    rng = random.Random(args.seed)
    dom = DOMAINS[args.mode]
    pts = []
    for _ in range(args.n_uniform):
        pts.append(tuple(rand_in(lo, hi, rng) for lo, hi, _n in dom))
    order = np.argsort(-leaves[:, 19])[:args.n_top]
    for r in order:
        row = leaves[r]
        pts.append(tuple(rand_in(row[3 + 2 * k], row[4 + 2 * k], rng) for k in range(3)))
    if args.mode == 'chart':
        pts = [p for p in pts if p[0] > 0]
    t0 = time.time()
    with Pool(args.workers, initializer=init, initargs=(args.mode, leaves)) as pool:
        results = pool.map(check_point, pts, chunksize=4)
    flat = [d for lst in results for d in lst]
    nu = sum(len(lst) for lst in results[:args.n_uniform])
    bad = [d for d in flat if d['status'] != 'ok']
    ok = [d for d in flat if d['status'] == 'ok']
    print(f"check_leaves {args.mode}: {len(leaves)} leaves, {len(pts)} parameter samples "
          f"({args.n_uniform} uniform, {args.n_top} in the top-bound leaves), {time.time() - t0:.1f}s")
    print(f"  Z points found: {len(flat)} ({nu} from uniform samples)")
    print(f"  ok: {len(ok)}   NOT covered / bad: {len(bad)}")
    if ok:
        m = min(ok, key=lambda d: d['margin'])
        mx = max(ok, key=lambda d: d['trS'])
        print(f"  largest trS(y_leaf) at a sampled Z point: {mx['trS']:.6f} (leaf bound {mx['ub']:.6f}) "
              f"at p={['%.6f' % v for v in mx['p']]}, dep={mx['t']:.6f}")
        print(f"  smallest margin (certified bound - trS): {m['margin']:.3e} (trS {m['trS']:.6f}, "
              f"bound {m['ub']:.6f})")
    for d in bad[:20]:
        print("  BAD:", d)


if __name__ == '__main__':
    main()
