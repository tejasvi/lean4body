"""Check (4): floating-point survey of  tr G = min_y tr S(y)  over the normalized convex CCs (direct region b >= b0 and
blow-up chart b < b0).  Vectorized root finding of P = pos - neg in x (or xi), then tr S at the optimal y.  Not
rigorous: it shows how far tr G stays below the threshold 3/4 of bb.py, and where the certificate is tightest."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cert

rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
MODE = sys.argv[1] if len(sys.argv) > 1 else 'direct'
NS = 200000
NX = 400


def trS_opt(Q):
    D, T, B = Q['D'], Q['T'], Q['B']
    y = cert.opt_y(Q)
    return sum(D[k] * ((B[k][0] + y[0] * T[k]) ** 2 + (B[k][1] + y[1] * T[k]) ** 2) for k in cert.EDGES)


if MODE == 'direct':
    a = rng.uniform(0.5, 1.75, NS); b = rng.uniform(0.125, 1.0, NS); c = rng.uniform(0, 1.75, NS)
    lo, hi = -1.0, 1.0
    FUN = cert.direct
else:
    a = rng.uniform(0, 0.125, NS) ** 1.0; b = rng.uniform(-2, 1, NS); c = rng.uniform(0, 3, NS)
    lo, hi = -3.0, 3.0
    FUN = cert.chart
xs = np.linspace(lo, hi, NX)
res = []
with np.errstate(all='ignore'):
    Pv = []
    for x in xs:
        Q = FUN(a, b, c, np.full(NS, x))
        ok = np.all(np.array(Q['g']) >= 0, axis=0)
        Pv.append(np.where(ok, Q['pos'] - Q['neg'], np.nan))
    Pv = np.array(Pv)
    for k in range(NX - 1):
        sc = (Pv[k] * Pv[k + 1] < 0)
        idx = np.nonzero(sc)[0]
        if len(idx) == 0:
            continue
        xl = np.full(len(idx), xs[k]); xh = np.full(len(idx), xs[k + 1])
        A, Bb, C = a[idx], b[idx], c[idx]
        fl = Pv[k][idx]
        for it in range(60):
            xm = 0.5 * (xl + xh)
            Q = FUN(A, Bb, C, xm)
            fm = Q['pos'] - Q['neg']
            left = fm * fl <= 0
            xh = np.where(left, xm, xh); xl = np.where(left, xl, xm); fl = np.where(left, fl, fm)
        xm = 0.5 * (xl + xh)
        Q = FUN(A, Bb, C, xm)
        t = trS_opt(Q)
        for j in range(len(idx)):
            res.append((t[j], A[j], Bb[j], C[j], xm[j]))
res.sort(reverse=True)
t = np.array([r[0] for r in res])
print(MODE, "roots", len(res), "of", NS, "samples;  max trS_opt %.6f  quantiles 50/90/99/99.9: %s" % (
    t.max(), np.round(np.quantile(t, [0.5, 0.9, 0.99, 0.999]), 4)))
for r in res[:12]:
    print("  trS %.6f at (%.5f, %.5f, %.5f, x=%.6f)" % r)
