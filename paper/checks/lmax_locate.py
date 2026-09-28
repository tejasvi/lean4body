"""Floating-point survey of lambda_max(G) on E, for the remarks after Conjecture 7.3 (not part of the proof).
The sampler is that of fig_lmax.py, with 100000 parameter points per region and another seed.  For each region the
script prints the number of CCs found, the largest values of lambda_max(G) and tr G with the coordinates (a, b, c, x)
where they occur, and the five CCs with the largest lambda_max(G)."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'work', 'rig'))
import cert
rng = np.random.default_rng(5)


def survey(MODE, NS, NX=400):
    if MODE == 'direct':
        p1 = rng.uniform(0.5, 1.75, NS); p2 = rng.uniform(0.125, 1.0, NS); p3 = rng.uniform(0, 1.75, NS)
        lo, hi, FUN = -1.0, 1.0, cert.direct
    else:
        p1 = 0.125 * rng.uniform(0, 1, NS) ** 3; p2 = rng.uniform(-2, 1, NS); p3 = rng.uniform(0, 3, NS)
        lo, hi, FUN = -3.0, 3.0, cert.chart
    xs = np.linspace(lo, hi, NX)
    out = []
    with np.errstate(all='ignore'):
        Pv = []
        for x in xs:
            Q = FUN(p1, p2, p3, np.full(NS, x))
            ok = np.all(np.array(Q['g']) >= 0, axis=0)
            Pv.append(np.where(ok, Q['pos'] - Q['neg'], np.nan))
        Pv = np.array(Pv)
        for k in range(NX - 1):
            idx = np.nonzero(Pv[k] * Pv[k + 1] < 0)[0]
            if len(idx) == 0:
                continue
            xl = np.full(len(idx), xs[k]); xh = np.full(len(idx), xs[k + 1])
            A, B, C = p1[idx], p2[idx], p3[idx]; fl = Pv[k][idx]
            for it in range(60):
                xm = 0.5 * (xl + xh); Q = FUN(A, B, C, xm); fm = Q['pos'] - Q['neg']; left = fm * fl <= 0
                xh = np.where(left, xm, xh); xl = np.where(left, xl, xm); fl = np.where(left, fl, fm)
            X = 0.5 * (xl + xh)
            Q = FUN(A, B, C, X)
            y = cert.opt_y(Q)
            S11, S12, S22 = cert.Smat(Q, y)
            tr = S11 + S22; det = S11 * S22 - S12 ** 2
            lmax = tr / 2 + np.sqrt(np.maximum(tr * tr / 4 - det, 0))
            if MODE == 'chart':
                a, b, c, xx = 1 + A * B, A, A * C, 0.5 + A * X
            else:
                a, b, c, xx = A, B, C, X
            out.append(np.stack([a, b, c, xx, lmax, tr], axis=1))
    return np.concatenate(out)


for mode in ('direct', 'chart'):
    R = survey(mode, 100000)
    R = R[np.isfinite(R).all(axis=1)]
    i = np.argmax(R[:, 4]); j = np.argmax(R[:, 5])
    print("%s: %d CCs; max lambda_max %.6f at (a,b,c,x) = %s (tr there %.4f); max tr %.6f at %s" % (
        mode, len(R), R[i, 4], np.round(R[i, :4], 5), R[i, 5], R[j, 5], np.round(R[j, :4], 5)))
    top = R[np.argsort(-R[:, 4])[:5]]
    print("   top lambda_max rows (a,b,c,x,lmax,tr):\n", np.round(top, 5))
