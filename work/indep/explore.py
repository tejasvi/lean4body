"""Float exploration of Z (not part of the proof): shape of Z, number of x-roots per (a,b,c),
size of F_x at the roots, min_y trS and the size of grad trS (fixed y) on Z.
Usage: explore.py direct|chart N
"""
import sys
import numpy as np

import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from formulas import prims_direct, prims_chart, certificate, trS, optimal_Y, Ftil, region_scaled

np.seterr(all='ignore')


def xrange_direct(a, b, c):
    lo = np.maximum.reduce([(a * a - 2 * b - b * b) / (2 * a), (1 - 2 * a * c - c * c) / (2 * a),
                            -(1 - b) / (2 * c), (c - a) / (2 * b), -np.ones_like(a)])
    hi = np.minimum.reduce([(a - c) / 2, (1 - b) / (2 * a), np.ones_like(a)])
    return lo, hi


def xrange_chart(b, al, ga):
    u = 1 + b * al
    lo = np.maximum(al / 2 - (2 + b) / (2 * u), -ga - (b * ga * ga + al) / (2 * u))
    hi = np.minimum((al - ga) / 2, -(1 + al) / (2 * u))
    return lo, hi


def roots(fun, P1, P2, P3, lo, hi, K=60, iters=60):
    """all sign changes of fun(P1,P2,P3,x) on a K-grid of [lo,hi], refined by bisection"""
    n = len(P1)
    ts = np.linspace(0, 1, K)
    X = lo[:, None] + (hi - lo)[:, None] * ts[None, :]
    Fv = fun(np.repeat(P1[:, None], K, 1), np.repeat(P2[:, None], K, 1), np.repeat(P3[:, None], K, 1), X)
    sc = np.sign(Fv[:, :-1]) * np.sign(Fv[:, 1:]) < 0
    ii, kk = np.nonzero(sc)
    xa, xb = X[ii, kk], X[ii, kk + 1]
    fa = Fv[ii, kk]
    p1, p2, p3 = P1[ii], P2[ii], P3[ii]
    for _ in range(iters):
        xm = 0.5 * (xa + xb)
        fm = fun(p1, p2, p3, xm)
        left = np.sign(fm) == np.sign(fa)
        xa = np.where(left, xm, xa)
        fa = np.where(left, fm, fa)
        xb = np.where(left, xb, xm)
    return ii, 0.5 * (xa + xb)


def Fd(a, b, c, x):
    return Ftil(prims_direct(a, b, c, x))


def Fc(b, al, ga, xi):
    return Ftil(prims_chart(b, al, ga, xi))


def analyse(mode, N):
    if mode == 'direct':
        g1 = np.linspace(0.5, 1.75, N)
        g2 = np.linspace(0.125, 1, N)
        g3 = np.linspace(0.001, 1.75, N)
        prims, fun, xr = prims_direct, Fd, xrange_direct
        names = ('a', 'b', 'c', 'x')
    else:
        g1 = np.linspace(1e-9, 0.125, N)
        g2 = np.linspace(-2, 1, N)
        g3 = np.linspace(0.001, 3, N)
        prims, fun, xr = prims_chart, Fc, xrange_chart
        names = ('b', 'al', 'ga', 'xi')
    P1, P2, P3 = [m.ravel() for m in np.meshgrid(g1, g2, g3, indexing='ij')]
    lo, hi = xr(P1, P2, P3)
    ok = lo < hi
    print(f"{mode}: grid {N}^3 = {len(P1)}; with nonempty linear x-range: {ok.sum()}")
    P1, P2, P3, lo, hi = P1[ok], P2[ok], P3[ok], lo[ok], hi[ok]
    ii, xs = roots(fun, P1, P2, P3, lo, hi)
    p1, p2, p3 = P1[ii], P2[ii], P3[ii]
    reg = np.min(np.array(region_scaled(prims(p1, p2, p3, xs))), axis=0)
    inreg = reg >= -1e-12
    ii, xs, p1, p2, p3 = ii[inreg], xs[inreg], p1[inreg], p2[inreg], p3[inreg]
    cnt = np.bincount(ii, minlength=len(P1))
    print(f"  roots in region: {len(xs)}; grid points with 1 root: {(cnt == 1).sum()}, "
          f"2 roots: {(cnt == 2).sum()}, >2: {(cnt > 2).sum()}")
    for nm, v in zip(names, (p1, p2, p3, xs)):
        print(f"  {nm}: [{v.min():.4f}, {v.max():.4f}]")
    cert = certificate(prims(p1, p2, p3, xs))
    Y = optimal_Y(cert)
    tr = trS(cert, Y)
    print(f"  min_y trS on Z samples: max {np.nanmax(tr):.5f}, median {np.nanmedian(tr):.4f}, "
          f"min {np.nanmin(tr):.4f}")
    k = np.nanargmax(tr)
    print(f"  argmax: {names[0]}={p1[k]:.5f} {names[1]}={p2[k]:.5f} {names[2]}={p3[k]:.5f} {names[3]}={xs[k]:.5f}")
    D = np.array(cert[0])
    print(f"  min D_e on Z samples: {np.nanmin(D):.3e}; max D_e {np.nanmax(D):.3e}")
    # derivatives: F_x and grad trS at fixed Y
    h = 1e-6
    Fx = (fun(p1, p2, p3, xs + h) - fun(p1, p2, p3, xs - h)) / (2 * h)
    F1 = (fun(p1 + h, p2, p3, xs) - fun(p1 - h, p2, p3, xs)) / (2 * h)
    F2 = (fun(p1, p2 + h, p3, xs) - fun(p1, p2 - h, p3, xs)) / (2 * h)
    F3 = (fun(p1, p2, p3 + h, xs) - fun(p1, p2, p3 - h, xs)) / (2 * h)
    scale = np.sqrt(F1 ** 2 + F2 ** 2 + F3 ** 2 + Fx ** 2)
    print(f"  |F_x|/|grad F| on Z: min {np.nanmin(np.abs(Fx) / scale):.3e}, "
          f"5% quantile {np.nanquantile(np.abs(Fx) / scale, 0.05):.3e}")
    dxd = [np.abs(Fi / Fx) for Fi in (F1, F2, F3)]
    print(f"  |dx/d{names[0]}|,|dx/d{names[1]}|,|dx/d{names[2]}| 95% quantiles: "
          + ", ".join(f"{np.nanquantile(d, 0.95):.3f}" for d in dxd)
          + "; max: " + ", ".join(f"{np.nanmax(d):.3f}" for d in dxd))
    grads = []
    for dv in range(4):
        pp = [p1.copy(), p2.copy(), p3.copy(), xs.copy()]
        pm = [p1.copy(), p2.copy(), p3.copy(), xs.copy()]
        pp[dv] += h
        pm[dv] -= h
        tp = trS(certificate(prims(*pp)), Y)
        tm = trS(certificate(prims(*pm)), Y)
        grads.append((tp - tm) / (2 * h))
    for dv in range(4):
        g = np.abs(grads[dv])
        print(f"  |d trS/d{names[dv]}| (fixed Y): median {np.nanmedian(g):.3f}, 95% {np.nanquantile(g, 0.95):.3f}, "
              f"max {np.nanmax(g):.3f}")
    tot = [grads[dv] + grads[3] * (-(F1, F2, F3)[dv] / Fx) for dv in range(3)]
    for dv in range(3):
        g = np.abs(tot[dv])
        print(f"  |total d trS/d{names[dv]}| along Z: median {np.nanmedian(g):.3f}, 95% {np.nanquantile(g, 0.95):.3f}, "
              f"max {np.nanmax(g):.3f}")
    return p1, p2, p3, xs, tr


if __name__ == '__main__':
    analyse(sys.argv[1], int(sys.argv[2]))
