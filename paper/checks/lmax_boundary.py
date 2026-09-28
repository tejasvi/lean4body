"""Suprema of lambda_max(G) at fixed b, for the remarks after Conjecture 7.3 (floating point, not part of the proof).
Nelder-Mead in the direct coordinates (a, c) at fixed b, with x running over the zeros of P, started at the points
(a0, c0) listed below.  The maximizing points approach the boundary of C, where one of the masses tends to zero, so
the script prints the smallest g_i and the normalized masses m/M at the final point."""
import os, sys
import numpy as np
from scipy.optimize import minimize
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'work', 'rig'))
import cert


def xrange_(a, b, c):
    """interval of x allowed by g1, g2 > 0 and g3, ..., g6 >= 0."""
    lo = [(a * a - b * b - 2 * b) / (2 * a), (1 - c * c - 2 * a * c) / (2 * a), (c - a) / (2 * b)]
    if c > 0:
        lo.append((b - 1) / (2 * c))
    hi = [(a - c) / 2, (1 - b) / (2 * a)]
    return max(max(lo), -1.0), min(min(hi), 1.0)


def zeros(a, b, c, n=400):
    lo, hi = xrange_(a, b, c)
    if not lo < hi:
        return []
    xs = np.linspace(lo, hi, n)
    Pv = [cert.direct(a, b, c, x)['pos'] - cert.direct(a, b, c, x)['neg'] for x in xs]
    out = []
    for k in range(n - 1):
        if Pv[k] == 0 or Pv[k] * Pv[k + 1] < 0:
            xl, xh, fl = xs[k], xs[k + 1], Pv[k]
            for _ in range(80):
                xm = 0.5 * (xl + xh)
                Q = cert.direct(a, b, c, xm)
                fm = Q['pos'] - Q['neg']
                if fm * fl <= 0:
                    xh = xm
                else:
                    xl, fl = xm, fm
            out.append(0.5 * (xl + xh))
    return out


def lmax_at(a, b, c, x):
    Q = cert.direct(a, b, c, x)
    S11, S12, S22 = cert.Smat(Q, cert.opt_y(Q))
    tr = S11 + S22; det = S11 * S22 - S12 * S12
    return tr / 2 + np.sqrt(max(tr * tr / 4 - det, 0)), min(Q['g']), Q


def masses(Q):
    """normalized masses from m_i m_j proportional to A_i A_j / w_ij (Dziobek's equations)."""
    A, w = Q['A'], Q['w']
    mu = {e: abs(A[int(e[0]) - 1] * A[int(e[1]) - 1] / w[e]) for e in w}
    m1 = np.sqrt(mu['12'] * mu['13'] / mu['23'])
    m = np.array([m1, mu['12'] / m1, mu['13'] / m1, mu['14'] / m1])
    return m / m.sum()


def best(a, b, c):
    vals = [(lmax_at(a, b, c, x)[0], x) for x in zeros(a, b, c)]
    vals = [v for v in vals if np.isfinite(v[0])]
    return max(vals) if vals else (np.nan, np.nan)


if __name__ == "__main__":
    with np.errstate(all='ignore'):
        for b, a0, c0 in ((0.05, 1.02525, 0.04875), (0.1, 1.052344, 0.094694), (0.2, 1.10525, 0.17947),
                          (0.4, 1.24055, 0.31660), (0.8, 1.53743, 0.50211)):
            f = lambda p: -best(p[0], b, p[1])[0] if np.isfinite(best(p[0], b, p[1])[0]) else 1.0
            res = minimize(f, [a0, c0], method='Nelder-Mead', options=dict(xatol=1e-11, fatol=1e-13, maxiter=800))
            a, c = res.x
            v, x = best(a, b, c)
            L, gmin, Q = lmax_at(a, b, c, x)
            print("b=%.2f  start %.6f  end %.6f  (1/3-end)/b^2=%.4f  at a=%.6f c=%.6f x=%.6f  min g=%.1e  m/M=%s" % (
                b, best(a0, b, c0)[0], L, (1 / 3 - L) / b ** 2, a, c, x, gmin,
                np.array2string(masses(Q), precision=3)), flush=True)
