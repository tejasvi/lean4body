"""Check (3), in computer algebra: the blow-up chart formulas (cert.chart) against the direct formulas (cert.direct).
Under a = 1 + b al, c = b ga, x = 1/2 + b xi, with b > 0, the script checks for the six edges e that
    T_e(chart) = tau_e(direct) / b,   beta_e(chart) = beta_e(direct),
and that g_i = b h_i for i = 1, ..., 4 and g_i = h_i for i = 5, 6.  The square roots are replaced by independent
positive symbols (r_k = sqrt(R) for each new radicand R, rho = sqrt(rho^2), s = sqrt(1 - x^2), and
sqrt(R_34) = b rho), and psi(X, Y) by (X^{-3/2} - Y^{-3/2}) / (Y - X) written with these symbols, which is the same
function.  Every difference is then a rational function of (b, al, ga, xi) and these symbols, and sympy.cancel must
reduce it to 0.  D_e and Ph = b P are compared by chart_identity_psi.py."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'work', 'rig'))
import sympy as sp
import cert

bb, al, ga, xi = sp.symbols('b alpha gamma xi', positive=True)
rho, s = sp.symbols('rho s', positive=True)
a, c, x = 1 + bb * al, bb * ga, sp.Rational(1, 2) + bb * xi
rho2 = sp.expand(1 - 2 * ga * x + ga ** 2)
ROOTS = {sp.expand(1 - x * x): s, sp.expand(rho2): rho, sp.expand(bb ** 2 * rho2): bb * rho}
names = iter(('r_%d' % k for k in range(1, 20)))


def sq(v):
    key = sp.expand(v)
    if key not in ROOTS:
        ROOTS[key] = sp.Symbol(next(names), positive=True)
    return ROOTS[key]


def hd(X, Y):
    sx, sy = sq(X), sq(Y)
    return (1 / (X * sx) - 1 / (Y * sy)) / (Y - X)


cert.sq = sq
cert.hd = hd
Qd = cert.direct(a, bb, c, x)
Qc = cert.chart(bb, al, ga, xi)
print("square-root symbols:", {str(v): sp.factor(k) for k, v in ROOTS.items()}, flush=True)


def zero(e):
    return sp.cancel(sp.together(e)) == 0


ok = True
for e in cert.EDGES:
    for task, r in (('T' + e, zero(Qd['T'][e] - bb * Qc['T'][e])),
                    ('Bx' + e, zero(Qd['B'][e][0] - Qc['B'][e][0])),
                    ('By' + e, zero(Qd['B'][e][1] - Qc['B'][e][1]))):
        print(task, "identical:", r, flush=True)
        ok &= r
r = all(zero(Qd['g'][i] - (bb if i < 4 else 1) * Qc['g'][i]) for i in range(6))
print("g identical:", r, flush=True)
print("ALL IDENTICAL" if ok and r else "MISMATCH")
