"""Branch-and-bound verifier for SPEC.md.

usage: verify.py direct|chart [--thresh P/Q] [--workers N] [--budget K] [--select i:j:k,...]
                              [--log FILE] [--json FILE]

For every box B of the initial grid (and recursively its halves):
  (i)   enclose all x (xi) such that B x {x} meets the region and F = 0 (chart: F/b^2 = 0)
        [core.enclose_x: linear region bounds, interval Newton with a mean-value numerator,
        bisection in x]; then, per x-piece X, a sheared enclosure eps = x - x_lin(p) in E
        [core.shear], x_lin being the float linearisation of the root at the centre of B;
  (ii)  nothing left -> discarded;
  (iii) else, per piece: Y = argmin trS at the centre (float); rigorous upper bound of trS(Y)
        on {(p, x_lin(p) + eps): p in B, eps in E}  (min of the interval value and of the
        mean-value form in (p, eps), arb) must be < thresh  [core.trS_upper_sheared];
  (iv)  otherwise bisect B (the dimension with the largest estimated contribution to the
        mean-value remainder).
"""
import argparse
import json
import math
import os
import sys
import time
from fractions import Fraction
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core
from core import enclose_x, center_data, shear, choose_Y, trS_upper_sheared
from flint import arb, ctx

DOMAINS = {
    # (lo, hi, number of initial cells) per coordinate; all cell endpoints are dyadic
    'direct': [(0.5, 1.75, 20), (0.125, 1.0, 14), (0.0, 1.75, 28)],
    'chart': [(0.0, 0.125, 8), (-2.0, 1.0, 24), (0.0, 3.0, 24)],
}
NAMES = {'direct': ('a', 'b', 'c', 'x'), 'chart': ('b', 'al', 'ga', 'xi')}
MIN_REL = 2.0 ** -18          # give up on a box when all its sides are below this x initial cell


def initial_boxes(mode):
    dom = DOMAINS[mode]
    cells = []
    for lo, hi, n in dom:
        w = (Fraction(hi) - Fraction(lo)) / n
        assert w.denominator & (w.denominator - 1) == 0, "non-dyadic grid"
        cells.append([(float(Fraction(lo) + k * w), float(Fraction(lo) + (k + 1) * w)) for k in range(n)])
    boxes = []
    for i, c1 in enumerate(cells[0]):
        for j, c2 in enumerate(cells[1]):
            for k, c3 in enumerate(cells[2]):
                boxes.append(((i, j, k), [c1, c2, c3]))
    return boxes


def process_initial(args):
    mode, idx, B0, thresh_str, node_budget, want_leaves = args
    ctx.prec = core.PREC
    thresh = arb(Fraction(thresh_str).numerator) / Fraction(thresh_str).denominator
    assert thresh.is_exact()
    w0 = [hi - lo for lo, hi in B0]
    t0 = time.time()
    st = dict(idx=idx, accepted=0, pieces=0, disc_region=0, disc_F=0, splits=0, fail=[],
              maxub=None, maxub_box=None, hist={}, maxdepth=0, minrelw=1.0, nodes=0, leaves=[])
    stack = [(B0, 0)]
    while stack:
        B, depth = stack.pop()
        st['nodes'] += 1
        st['maxdepth'] = max(st['maxdepth'], depth)
        rel = [(hi - lo) / w for (lo, hi), w in zip(B, w0)]
        Xs, status, cd = enclose_x(mode, B)
        ok = True
        split_dim = None
        recs = []
        live = 0
        for X in Xs:
            xc, phi = cd
            if not (float(X.lower()) <= xc <= float(X.upper())):
                xc, phi = center_data(mode, B, X)
            E, contrib = shear(mode, B, X, xc, phi)
            if E is None:
                continue                      # proved: no root of F in B x X
            live += 1
            Y = choose_Y(mode, B, xc + float(E.mid()))
            ub, gmag = trS_upper_sheared(mode, B, E, xc, phi, Y)
            if ub is not None and ub < thresh:
                u = float(ub.upper())
                if arb(u) < ub:                 # report a float that is >= the certified bound
                    u = math.nextafter(u, math.inf)
                recs.append((u, [float(X.lower()), float(X.upper())],
                             [xc] + list(phi), [float(E.lower()), float(E.upper())], list(Y)))
                continue
            ok = False
            if gmag is not None and contrib is not None:
                rads = [(hi - lo) / 2 for lo, hi in B]
                score = [gmag[k] * rads[k] + gmag[3] * contrib[k] for k in range(3)]
                split_dim = max(range(3), key=lambda k: score[k])
            break
        if ok and live == 0:
            if status == 'region':
                st['disc_region'] += 1
            else:
                st['disc_F'] += 1
            continue
        if ok:
            st['accepted'] += 1
            st['pieces'] += live
            st['minrelw'] = min(st['minrelw'], max(rel))
            for (u, Xb, lin, Eb, Y) in recs:
                if want_leaves:
                    st['leaves'].append(list(idx) + [v for lh in B for v in lh] + Xb + lin + Eb + Y + [u])
                key = f"{math.floor(u * 100) / 100:.2f}"
                st['hist'][key] = st['hist'].get(key, 0) + 1
                if st['maxub'] is None or u > st['maxub']:
                    st['maxub'], st['maxub_box'] = u, dict(box=B, xpiece=Xb, xlin=lin, eps=Eb, Y=Y)
            continue
        if max(rel) < MIN_REL or st['nodes'] > node_budget:
            st['fail'].append(dict(box=B, xpieces=[[float(X.lower()), float(X.upper())] for X in Xs],
                                   depth=depth, reason='min size' if max(rel) < MIN_REL else 'budget'))
            continue
        if split_dim is None:
            split_dim = max(range(3), key=lambda k: rel[k])
        lo, hi = B[split_dim]
        mid = 0.5 * (lo + hi)
        assert lo < mid < hi and arb(mid) == (arb(lo) + arb(hi)) / 2
        B1 = list(B)
        B2 = list(B)
        B1[split_dim] = (lo, mid)
        B2[split_dim] = (mid, hi)
        stack.append((B2, depth + 1))
        stack.append((B1, depth + 1))
        st['splits'] += 1
    st['time'] = time.time() - t0
    return st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['direct', 'chart'])
    ap.add_argument('--thresh', default='3/4')
    ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--limit', type=int, default=0, help='only the first K initial boxes (testing)')
    ap.add_argument('--select', default='', help='comma list of i:j:k initial boxes (testing)')
    ap.add_argument('--log', default='')
    ap.add_argument('--json', default='')
    ap.add_argument('--budget', type=int, default=200000, help='max B&B nodes per initial box')
    ap.add_argument('--leaves', default='', help='save accepted leaves (.npy) for check_leaves.py')
    args = ap.parse_args()
    assert args.workers <= 16
    logf = open(args.log, 'a') if args.log else sys.stdout

    def log(msg):
        logf.write(msg + '\n')
        logf.flush()

    boxes = initial_boxes(args.mode)
    if args.select:
        sel = {tuple(int(t) for t in s.split(':')) for s in args.select.split(',')}
        boxes = [b for b in boxes if b[0] in sel]
    if args.limit:
        boxes = boxes[:args.limit]
    log(f"=== verify.py {args.mode}  thresh={args.thresh}  prec={core.PREC} bits  workers={args.workers}  "
        f"initial boxes={len(boxes)}  domain={DOMAINS[args.mode]}  start {time.strftime('%Y-%m-%d %H:%M:%S')}")
    t0 = time.time()
    tot = dict(accepted=0, pieces=0, disc_region=0, disc_F=0, splits=0, fails=0, cpu=0.0, nodes=0)
    maxub, maxub_box = None, None
    hist = {}
    fails = []
    maxdepth = 0
    minrelw = 1.0
    done = 0
    tasks = [(args.mode, idx, B, args.thresh, args.budget, bool(args.leaves)) for idx, B in boxes]
    leaves = []
    with Pool(args.workers) as pool:
        for st in pool.imap_unordered(process_initial, tasks, chunksize=1):
            done += 1
            for k in ('accepted', 'pieces', 'disc_region', 'disc_F', 'splits', 'nodes'):
                tot[k] += st[k]
            tot['fails'] += len(st['fail'])
            tot['cpu'] += st['time']
            fails.extend(st['fail'])
            leaves.extend(st['leaves'])
            maxdepth = max(maxdepth, st['maxdepth'])
            if st['accepted']:
                minrelw = min(minrelw, st['minrelw'])
            for k, v in st['hist'].items():
                hist[k] = hist.get(k, 0) + v
            if st['maxub'] is not None and (maxub is None or st['maxub'] > maxub):
                maxub, maxub_box = st['maxub'], st['maxub_box']
            if st['fail']:
                log(f"  FAIL in initial box {st['idx']}: {len(st['fail'])} unproved leaves, e.g. {st['fail'][0]}")
            if done % 250 == 0 or done == len(tasks):
                log(f"  [{done}/{len(tasks)}] {time.time() - t0:8.1f}s  accepted={tot['accepted']} "
                    f"disc_region={tot['disc_region']} disc_F={tot['disc_F']} fails={tot['fails']} "
                    f"maxub={maxub}")
    wall = time.time() - t0
    log(f"=== DONE {args.mode}: wall {wall:.1f}s, cpu {tot['cpu']:.1f}s")
    log(f"    leaves accepted (trS bound < {args.thresh}): {tot['accepted']}  (x-pieces {tot['pieces']})")
    log(f"    leaves discarded: region {tot['disc_region']}, F != 0 {tot['disc_F']}")
    log(f"    bisections: {tot['splits']}, max depth {maxdepth}, smallest accepted box (rel. to initial cell) {minrelw:.3g}")
    log(f"    unproved leaves: {tot['fails']}")
    log(f"    max certified upper bound of trS over accepted leaves: {maxub!r}")
    log(f"      attained on box {maxub_box}")
    log(f"    histogram of leaf upper bounds: {dict(sorted(hist.items()))}")
    if args.leaves:
        # columns: i j k  p1lo p1hi p2lo p2hi p3lo p3hi  Xlo Xhi  xc phi1 phi2 phi3  Elo Ehi  Y1 Y2  ub
        np.save(args.leaves, np.array(leaves, dtype=np.float64))
        log(f"    accepted leaves saved to {args.leaves}")
    if args.json:
        with open(args.json, 'w') as f:
            json.dump(dict(mode=args.mode, thresh=args.thresh, prec=core.PREC, workers=args.workers,
                           initial_boxes=len(boxes), domain=DOMAINS[args.mode], wall=wall, **tot,
                           maxdepth=maxdepth, minrelw=minrelw, maxub=maxub, maxub_box=maxub_box,
                           hist=dict(sorted(hist.items())), fail_list=fails), f, indent=1)


if __name__ == '__main__':
    main()
