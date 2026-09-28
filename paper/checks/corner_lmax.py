"""Floating-point (mpmath, 40 digits) study of lambda_max(G) near the corner b = 0 of the chart, for the remarks after
Conjecture 7.3 (not part of the proof).  For small b > 0 and (alpha, gamma) in the admissible region, point() solves
Ph(b, alpha, gamma, xi) = 0 for xi near xi0 = 3(alpha gamma - alpha - gamma) / (2 (1 + gamma)) and computes G with
the certified formulas (cert.chart).  For fixed b and alpha the script maximizes lambda_max(G) over gamma by golden
section search, first for |alpha| <= 0.45 and then for alpha near the edges alpha = +-1/2 of the admissible region,
and prints the maximum minus 1/3, also divided by b^2."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'work', 'rig'))
import mpmath as mp
import cert
mp.mp.dps = 40
cert.sq = mp.sqrt


def Ph(b, al, ga, xi):
    g, P = cert.chart_P(mp.mpf(b), mp.mpf(al), mp.mpf(ga), xi)
    return P


def point(b, al, ga):
    xi0 = 3 * (al * ga - al - ga) / (2 * (1 + ga))
    xi = mp.findroot(lambda t: Ph(b, al, ga, t), mp.mpf(xi0))
    Q = cert.chart(mp.mpf(b), mp.mpf(al), mp.mpf(ga), xi)
    y = cert.opt_y(Q)
    S11, S12, S22 = cert.Smat(Q, y)
    tr, det = S11 + S22, S11 * S22 - S12 ** 2
    lmax = tr / 2 + mp.sqrt(tr * tr / 4 - det)
    gs = Q['g']
    ok = all(v > 0 for v in gs)
    return lmax, tr, ok


def golden(f, lo, hi, it=60):
    g = (mp.sqrt(5) - 1) / 2
    a, b = mp.mpf(lo), mp.mpf(hi)
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(it):
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = f(d)
    return (a + b) / 2, max(fc, fd)


def report(b, al, lo, hi, it, afmt):
    gs, v = golden(lambda ga: point(mp.mpf(b), al, ga)[0], lo, hi, it)
    d = v - mp.mpf(1) / 3
    print(("b=%s alpha=" + afmt + ": argmax gamma=%s  max lmax - 1/3 = %s  (/b^2 = %s)") % (
        b, al, mp.nstr(gs, 10), mp.nstr(d, 6), mp.nstr(d / mp.mpf(b) ** 2, 6)), flush=True)


if __name__ == '__main__':
    for b in ('3e-2', '1e-2', '1e-3'):
        for al in (-0.45, -0.3, -0.15, 0.0, 0.15, 0.3, 0.45):
            report(b, al, 0.8, 1.2, 60, '%.2f')
    for b in ('1e-2', '1e-3'):
        for al in (-0.49, -0.499, -0.4999, 0.49, 0.499, 0.4999):
            report(b, al, 0.9, 1.1, 50, '%.4f')
