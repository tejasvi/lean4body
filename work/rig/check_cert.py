"""Checks of a certificate file written by bb.py (run from work/rig):
    python check_cert.py CERTFILE [nsample=2000] [seed=1]
1. Structure: the initial boxes are exactly the grid of the header (mode, N, B0); inside each initial box the leaf
   paths form a complete binary tree (every internal node has both children, split in the same dimension), so the
   leaves tile the initial box; every accepted leaf has tr S upper bounds below the threshold.  Also prints the
   number of accepted leaves with one, two, ... runs of the dependent coordinate.
2. Sampling: random central configurations in the covered region (roots of the Dziobek function found in floats and
   polished in mpmath at 40 digits), each located in the tree.  Checks: the leaf is accepted (not discarded), the
   dependent coordinate lies in one of its x-runs, tr S(y) at the configuration with that run's y is below the
   certified bound, and lambda_max(G) computed independently (cc4.Gmat: cotangent formula for beta, optimal y)
   does not exceed it.  In the chart, tr S is also recomputed in the direct coordinates with y/b."""
import sys, os, math, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import mpmath as mp
import cert
import cc4

mp.mp.dps = 40
FN = sys.argv[1]
NS = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
random.seed(int(sys.argv[3]) if len(sys.argv) > 3 else 1)

# ---------- parse ----------
boxes = {}
hdr = None
cur = None
for line in open(FN):
    if line.startswith('#'):
        t = line.split()
        hdr = dict(mode=t[2], N=(int(t[4]), int(t[5]), int(t[6])), B0=float(t[8]), prec=int(t[10]), thresh=float(t[12]))
        continue
    t = line.split()
    if t[0] == 'B':
        key = tuple(float(v) for v in t[1:])
        assert key not in boxes, "initial box listed twice"
        cur = boxes[key] = {}
    elif t[0] == 'D':
        path = '' if t[1] == '-' else t[1]
        assert path not in cur
        cur[path] = None
    elif t[0] == 'A':
        path = '' if t[1] == '-' else t[1]
        assert path not in cur
        runs = []
        for tok in t[2:]:
            x0, x1, y0, y1, tu = (float(v) for v in tok.split(','))
            runs.append((x0, x1, y0, y1, tu))
        assert runs
        cur[path] = runs
MODE, N, B0, TH = hdr['mode'], hdr['N'], hdr['B0'], hdr['thresh']
RANGES = [(0.5, 1.75), (B0, 1.0), (0.0, 1.75)] if MODE == 'direct' else [(0.0, B0), (-2.0, 1.0), (0.0, 3.0)]

# ---------- 1. structure ----------
pts = []
for d in range(3):
    p = [RANGES[d][0] + (RANGES[d][1] - RANGES[d][0]) * k / N[d] for k in range(N[d] + 1)]
    assert p[0] == RANGES[d][0] and p[-1] == RANGES[d][1] and all(u < v for u, v in zip(p, p[1:]))
    pts.append(p)
grid = {(pts[0][i], pts[0][i + 1], pts[1][j], pts[1][j + 1], pts[2][k], pts[2][k + 1])
        for i in range(N[0]) for j in range(N[1]) for k in range(N[2])}
assert set(boxes) == grid, "initial boxes differ from the grid"


def tokens(path):
    assert len(path) % 2 == 0
    return [(int(path[i]), path[i + 1]) for i in range(0, len(path), 2)]


def complete(leaves, prefix):
    """leaves: set of paths (as tuples of tokens) under this node, relative to it"""
    if leaves == {()}:
        return True
    if () in leaves:
        return False               # a leaf that also has descendants
    dims = {p[0][0] for p in leaves}
    if len(dims) != 1:
        return False
    L = {p[1:] for p in leaves if p[0][1] == 'L'}
    R = {p[1:] for p in leaves if p[0][1] == 'R'}
    return bool(L) and bool(R) and complete(L, prefix) and complete(R, prefix)


sys.setrecursionlimit(10000)
nA = nD = 0
tmax = 0.0
depth = 0
nruns = {}
for key, recs in boxes.items():
    assert complete({tuple(tokens(p)) for p in recs}, ()), "incomplete tree in %r" % (key,)
    for p, runs in recs.items():
        depth = max(depth, len(p) // 2)
        if runs is None:
            nD += 1
        else:
            nA += 1
            nruns[len(runs)] = nruns.get(len(runs), 0) + 1
            for r in runs:
                assert r[0] <= r[1] and r[4] < TH
                tmax = max(tmax, r[4])
print("%s: %d initial boxes = the %dx%dx%d grid; %d accepted + %d discarded leaves; every tree complete; max depth %d;"
      " accepted leaves by number of runs %s; max certified tr S %.9f < %g" % (
          MODE, len(boxes), N[0], N[1], N[2], nA, nD, depth, dict(sorted(nruns.items())), tmax, TH), flush=True)


# ---------- 2. sampling ----------
def leaf_of(p):
    """the leaf containing the point p (3 floats) and its box; points on a face may go either way"""
    idx = []
    for d in range(3):
        k = min(int((p[d] - RANGES[d][0]) / (RANGES[d][1] - RANGES[d][0]) * N[d]), N[d] - 1)
        while pts[d][k] > p[d]:
            k -= 1
        while pts[d][k + 1] < p[d]:
            k += 1
        idx.append(k)
    key = (pts[0][idx[0]], pts[0][idx[0] + 1], pts[1][idx[1]], pts[1][idx[1] + 1], pts[2][idx[2]], pts[2][idx[2] + 1])
    recs = boxes[key]
    bx = [[key[0], key[1]], [key[2], key[3]], [key[4], key[5]]]
    path = ''
    while path not in recs:
        nxt = [q for q in recs if q.startswith(path)][0][len(path)]
        k = int(nxt)
        m = 0.5 * (bx[k][0] + bx[k][1])
        if p[k] <= m:
            path += '%dL' % k; bx[k][1] = m
        else:
            path += '%dR' % k; bx[k][0] = m
    return path, recs[path], bx


FUN = cert.direct if MODE == 'direct' else cert.chart
FP = cert.direct_P if MODE == 'direct' else cert.chart_P


def Pval(p, x):
    Q = FUN(*[mp.mpf(v) for v in p], mp.mpf(x))
    return Q['pos'] - Q['neg'], min(Q['g'])


def roots(p):
    """all x (or xi) in the scan range with P = 0 and all g >= 0"""
    lo, hi = (-1.0, 1.0) if MODE == 'direct' else (-4.0, 4.0)
    n = 400
    xs = [lo + (hi - lo) * k / n for k in range(n + 1)]
    vals = []
    for x in xs:
        try:
            v, g = Pval(p, x)
            vals.append(v if g >= 0 else None)
        except (ZeroDivisionError, ValueError):
            vals.append(None)
    out = []
    for k in range(n):
        if vals[k] is not None and vals[k + 1] is not None and vals[k] * vals[k + 1] < 0:
            r = mp.findroot(lambda x: Pval(p, x)[0], (mp.mpf(xs[k]), mp.mpf(xs[k + 1])), solver='anderson')
            if Pval(p, r)[1] >= 0:
                out.append(r)
    return out


def trS_at(Q, y):
    return sum(Q['D'][e] * ((Q['B'][e][0] + y[0] * Q['T'][e]) ** 2 + (Q['B'][e][1] + y[1] * Q['T'][e]) ** 2)
               for e in cert.EDGES)


ns = nroot = 0
worst_gap = 1.0
maxrat = 0.0
maxdiff = 0.0
while ns < NS:
    p = [random.uniform(*RANGES[d]) for d in range(3)]
    if MODE == 'chart' and p[0] == 0.0:
        continue
    ns += 1
    for x in roots(p):
        nroot += 1
        path, runs, bx = leaf_of(p)
        assert runs is not None, "a central configuration lies in a DISCARDED leaf: %r x=%s" % (p, x)
        inrun = [r for r in runs if r[0] <= x <= r[1]]
        assert inrun, "x = %s outside every run of leaf %r: %r" % (x, bx, runs)
        Q = FUN(*[mp.mpf(v) for v in p], x)
        best = min(trS_at(Q, (mp.mpf(r[2]), mp.mpf(r[3]))) for r in inrun)
        bound = min(r[4] for r in inrun)
        assert best <= bound, (p, x, best, bound)
        worst_gap = min(worst_gap, float(bound - best))
        # independent lambda_max(G) at the same configuration
        if MODE == 'direct':
            a, b, c, xx = mp.mpf(p[0]), mp.mpf(p[1]), mp.mpf(p[2]), x
        else:
            b = mp.mpf(p[0]); a = 1 + b * p[1]; c = b * p[2]; xx = mp.mpf(1) / 2 + b * x
            Qd = cert.direct(a, b, c, xx)
            for r in inrun:
                y = (mp.mpf(r[2]), mp.mpf(r[3]))
                t1 = trS_at(Q, y)
                t2 = trS_at(Qd, (y[0] / b, y[1] / b))
                maxdiff = max(maxdiff, float(abs(t1 - t2) / t1))
        t = mp.acos(xx)
        P = cc4.pos(a, b, c, t)
        G, _ = cc4.Gmat(P)
        lm = cc4.lam_max(G)
        assert lm <= best * (1 + mp.mpf(10) ** -30), (p, x, lm, best)
        maxrat = max(maxrat, float(lm / best))
print("sampling: %d parameter samples, %d central configurations; all in accepted leaves and inside an x-run;"
      " min (certified bound - tr S at the point) = %.4f;  max lambda_max(G)/tr S = %.4f%s" % (
          ns, nroot, worst_gap, maxrat, ("; chart vs direct tr S max rel diff %.1e" % maxdiff) if MODE == 'chart' else ''))
