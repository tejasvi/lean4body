"""Rigorous verification of  tr S(y) < THRESH  on closure(E) (normalized convex four-body central configurations),
with THRESH from the environment (default 1; Proposition 5.1 uses THRESH=0.75).
Mean-value (first-order jet) enclosures in flint.arb; the dependent coordinate (x = cos theta, or xi in the blow-up
chart) is enclosed by interval Newton on the Dziobek function P (F = 0 <=> P = 0).
   direct : (a, b, c) in [0.5, 1.75] x [b0, 1] x [0, 1.75]
   chart  : (b, al, ga) in [0, b0] x [-2, 1] x [0, 3]        (a = 1 + b al, c = b ga, x = 1/2 + b xi)
usage: python bb.py direct|chart N1 N2 N3 [b0inv=8] [minw=1e-5] [nproc=32] [outfile]"""
import sys, os, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from multiprocessing import Pool
from flint import arb, ctx
import cert, jet

ctx.prec = int(os.environ.get('PREC', '64'))
MODE = sys.argv[1]
N = (int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
B0 = 1.0 / int(sys.argv[5]) if len(sys.argv) > 5 else 1.0 / 8
MINW = float(sys.argv[6]) if len(sys.argv) > 6 else 1e-5
BUDGET = 400000
NPROC = int(sys.argv[7]) if len(sys.argv) > 7 else 32
OUT = sys.argv[8] if len(sys.argv) > 8 else None
THRESH = float(os.environ.get('THRESH', '1.0'))


def iv(lo, hi):
    return arb(lo).union(arb(hi))


def up(v):
    """a float >= every point of the arb v"""
    f = float(v.upper())
    return f + 4e-16 * abs(f) + 1e-300


def dn(v):
    f = float(v.lower())
    return f - 4e-16 * abs(f) - 1e-300


def mid_rad(lo, hi):
    m = 0.5 * (lo + hi)
    r = max(up(arb(hi) - arb(m)), up(arb(m) - arb(lo)))
    return m, r


def xrange_direct(A, Bb, C):
    L = [(A * A - Bb * Bb - 2 * Bb) / (2 * A), (1 - C * C - 2 * A * C) / (2 * A), (C - A) / (2 * Bb)]
    if C > 0:
        L.append((Bb - 1) / (2 * C))
    U = [(A - C) / 2, (1 - Bb) / (2 * A)]
    return max(max(dn(l) for l in L), -1.0), min(min(up(u) for u in U), 1.0)


def xrange_chart(Bb, AL, GA):
    den = 2 * (1 + Bb * AL)
    L = [(AL - 2 - Bb * (1 - AL * AL)) / den, -(2 * GA + AL + Bb * (GA * GA + 2 * AL * GA)) / den]
    U = [-(1 + AL) / den, (AL - GA) / 2]
    return max(dn(l) for l in L), min(up(u) for u in U)


if MODE == 'direct':
    FUN, FUNP, XRANGE = cert.direct, cert.direct_P, xrange_direct
    RANGES = [(0.5, 1.75), (B0, 1.0), (0.0, 1.75)]
    WSCALE = (1.0, 1.0, 1.0)
else:
    FUN, FUNP, XRANGE = cert.chart, cert.chart_P, xrange_chart
    RANGES = [(0.0, B0), (-2.0, 1.0), (0.0, 3.0)]
    WSCALE = (2.0, 1.0, 1.0)


def box4(bf, X):
    mr = [mid_rad(lo, hi) for lo, hi in list(bf) + [X]]
    mids = [m for m, r in mr]
    rads = [r for m, r in mr]
    boxes = [iv(lo, hi) for lo, hi in list(bf) + [X]]
    return mids, rads, boxes


def mvP(bf, X):
    """None if the 4-box contains no point of closure(C) with P = 0; else (enclosure of P(B, x_mid), P_x over box)"""
    mids, rads, boxes = box4(bf, X)
    gJ, PJ = FUNP(*jet.seeds(boxes))
    gm, Pm = FUNP(*[arb(m) for m in mids])
    for gj, g0 in zip(gJ, gm):
        if not isinstance(gj, jet.J):
            continue
        if jet.mv(g0, gj, rads) < 0:
            return None
    Pall = jet.mv(Pm, PJ, rads)
    if Pall > 0 or Pall < 0:
        return None
    Pd = Pm
    for i in range(3):
        Pd = Pd + PJ.d[i] * arb(0, rads[i])
    return Pd, PJ.d[3], mids[3], [float(PJ.d[i].mid()) for i in range(4)]


def xenclose(bf, stats):
    p = [iv(lo, hi) for lo, hi in bf]
    xl, xu = XRANGE(*p)
    if xl > xu:
        return None, None
    w = max(hi - lo for lo, hi in bf)
    tol = max(0.25 * w, 1e-13)
    work = [(xl, xu)]
    kept = []
    grad = None
    while work:
        X = work.pop()
        while True:
            stats['nP'] += 1
            res = mvP(bf, X)
            if res is None:
                X = None
                break
            Pd, dx, xm, g = res
            grad = g
            if dx > 0 or dx < 0:
                Nw = arb(xm) - Pd / dx
                lo, hi = max(X[0], dn(Nw)), min(X[1], up(Nw))
                if lo > hi:
                    X = None
                    break
                if hi - lo <= 0.7 * (X[1] - X[0]):
                    X = (lo, hi)
                    continue
                X = (lo, hi)
            if X[1] - X[0] <= tol:
                break
            m = 0.5 * (X[0] + X[1])
            work.append((X[0], m))
            X = (m, X[1])
        if X is not None:
            kept.append(X)
    if not kept:
        return None, None
    kept.sort()
    runs = [list(kept[0])]
    for lo, hi in kept[1:]:
        if lo <= runs[-1][1]:
            runs[-1][1] = max(runs[-1][1], hi)
        else:
            runs.append([lo, hi])
    return runs, grad


def mvcert(bf, X):
    mids, rads, boxes = box4(bf, X)
    Qf = FUN(*mids)
    y = cert.opt_y(Qf)
    y = (arb(y[0]), arb(y[1]))
    tJ = cert.trS(FUN(*jet.seeds(boxes)), y)
    tm = cert.trS(FUN(*[arb(m) for m in mids]), y)
    t = jet.mv(tm, tJ, rads)
    contrib = [float(abs(tJ.d[i]).upper()) * rads[i] for i in range(4)]
    return t, contrib, (float(y[0].mid()), float(y[1].mid()))


def process(bf0):
    stats = dict(acc=0, disc=0, fail=[], maxt=0.0, nP=0, ncert=0, minw=1.0, rec=[])
    stack = [(bf0, '')]
    while stack:
        bf, path = stack.pop()
        runs, grad = xenclose(bf, stats)
        if runs is None:
            stats['disc'] += 1
            stats['rec'].append(('D', path))
            continue
        ok = True
        worst = None
        tmax = 0.0
        cert_runs = []
        for X in runs:
            stats['ncert'] += 1
            t, contrib, yf = mvcert(bf, tuple(X))
            tu = up(t) if math.isfinite(float(t.mid())) and math.isfinite(float(t.rad())) else float('inf')
            if not (tu < THRESH):
                ok = False
                worst = contrib
                break
            tmax = max(tmax, tu)
            cert_runs.append((X[0], X[1], yf[0], yf[1], tu))
        if ok:
            stats['rec'].append(('A', path, cert_runs))
            stats['acc'] += 1
            stats['maxt'] = max(stats['maxt'], tmax)
            stats['minw'] = min(stats['minw'], max(hi - lo for lo, hi in bf))
            continue
        w = [hi - lo for lo, hi in bf]
        if max(w) < MINW or stats['acc'] + stats['disc'] > BUDGET:
            stats['fail'].append((bf, runs))
            continue
        # split: attribute the x-contribution to (a, b, c) through dx/di = -P_i / P_x; widest dimension if anything
        # is not finite; never split a dimension narrower than a quarter of the widest (aspect control)
        allowed = [i for i in range(3) if w[i] * WSCALE[i] >= 0.25 * max(w[j] * WSCALE[j] for j in range(3))]
        k = max(allowed, key=lambda i: w[i] * WSCALE[i])
        c = worst[:3]
        if all(math.isfinite(v) for v in worst) and grad is not None and grad[3] != 0 and all(math.isfinite(v) for v in grad):
            dxd = [abs(grad[i] / grad[3]) * w[i] for i in range(3)]
            sdx = sum(dxd) or 1.0
            c = [c[i] + worst[3] * dxd[i] / sdx for i in range(3)]
            k = max(allowed, key=lambda i: c[i])
        lo, hi = bf[k]
        m = 0.5 * (lo + hi)
        b1 = list(bf); b1[k] = (lo, m)
        b2 = list(bf); b2[k] = (m, hi)
        stack.append((b1, path + '%dL' % k)); stack.append((b2, path + '%dR' % k))
    return stats


def breakpoints():
    """grid breakpoints per dimension; checked to start and end exactly at the range ends and to increase, so the
    initial boxes (consecutive breakpoints, identical floats on both sides of every face) tile the full box"""
    pts = []
    for d in range(3):
        p = [RANGES[d][0] + (RANGES[d][1] - RANGES[d][0]) * k / N[d] for k in range(N[d] + 1)]
        assert p[0] == RANGES[d][0] and p[-1] == RANGES[d][1] and all(u < v for u, v in zip(p, p[1:]))
        pts.append(p)
    return pts


def grid():
    pts = breakpoints()
    return [[(pts[0][i], pts[0][i + 1]), (pts[1][j], pts[1][j + 1]), (pts[2][k], pts[2][k + 1])]
            for i in range(N[0]) for j in range(N[1]) for k in range(N[2])]


def work(bf):
    t0 = time.time()
    s = process(bf)
    s['time'] = time.time() - t0
    return bf, s


if __name__ == '__main__':
    boxes = grid()
    t0 = time.time()
    tot = dict(acc=0, disc=0, fail=[], maxt=0.0, nP=0, ncert=0, minw=1.0)
    done = 0
    slow = []
    fout = open(OUT, 'w') if OUT else None
    if fout:
        fout.write("# mode %s N %d %d %d B0 %r prec %d thresh %r\n" % (MODE, N[0], N[1], N[2], B0, ctx.prec, THRESH))
    with Pool(NPROC) as pool:
        for bf, s in pool.imap_unordered(work, boxes, chunksize=1):
            done += 1
            if fout:
                fout.write("B %s\n" % ' '.join(repr(v) for lh in bf for v in lh))
                for r in s['rec']:
                    if r[0] == 'D':
                        fout.write("D %s\n" % (r[1] or '-'))
                    else:
                        fout.write("A %s %s\n" % (r[1] or '-', ' '.join('%r,%r,%r,%r,%r' % c for c in r[2])))
            for key in ('acc', 'disc', 'nP', 'ncert'):
                tot[key] += s[key]
            tot['maxt'] = max(tot['maxt'], s['maxt'])
            tot['minw'] = min(tot['minw'], s['minw'])
            tot['fail'] += s['fail']
            slow.append((s['time'], bf, s['acc']))
            if s['fail']:
                print("FAIL in", bf, len(s['fail']), "e.g.", s['fail'][0], flush=True)
            if done % max(1, len(boxes) // 10) == 0:
                print("%d/%d  acc %d disc %d fail %d  max trS %.4f  minw %.2e  P-evals %d cert-evals %d  %.0fs" % (
                    done, len(boxes), tot['acc'], tot['disc'], len(tot['fail']), tot['maxt'], tot['minw'], tot['nP'], tot['ncert'],
                    time.time() - t0), flush=True)
    if fout:
        fout.close()
    slow.sort(reverse=True)
    print("slowest initial boxes:", [(round(t, 1), bf, n) for t, bf, n in slow[:5]])
    print("DONE %s: acc %d disc %d fail %d  max trS upper %.9f  minw %.2e  P-evals %d cert-evals %d  %.0fs" % (
        MODE, tot['acc'], tot['disc'], len(tot['fail']), tot['maxt'], tot['minw'], tot['nP'], tot['ncert'], time.time() - t0))
