"""Symbolic check of the b = 0 values in Proposition 7.2 (not part of the proof of Theorems A and B).
(1) h_i at xi0, the limits of D_e, T_e, beta_e, and the eigenvalues of G0 (sympy, exact).
(2) Numerical comparison with the certified chart formulas (cert.chart) at b = 1e-15 (mpmath, 50 digits).
(3) Membership check (all region functions > 0) at b = 1e-2 and 1e-3 for alpha near +-1/2 and gamma near 1, where
    corner_lmax.py maximizes lambda_max(G) over gamma."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'work', 'rig'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp

al, ga = sp.symbols('alpha gamma', real=True)
xi0 = 3 * (al * ga - al - ga) / (2 * (1 + ga))
h1 = sp.factor(2 - al + 2 * xi0)
h2 = sp.factor(2 * ga + 2 * xi0 + al)
h3 = sp.factor(al - ga - 2 * xi0)
h4 = sp.factor(-(1 + al + 2 * xi0))
print("h1 =", h1, "  h2 =", h2, "  h3 =", h3, "  h4 =", h4)
rho = sp.sqrt(ga ** 2 - ga + 1)
s3 = sp.sqrt(3)
q1 = sp.Matrix([1, 0]); q2 = sp.Matrix([sp.Rational(1, 2), s3 / 2])
A1h, A2h, A3, A4 = s3 / 4, -s3 * ga / 4, s3 / 4, -s3 / 4
# b = 0 limits of the w-hat_e (formulas of the Appendix with psi = 3/2, R_e = 1, x = 1/2)
w = {'12': sp.Rational(9, 4) * h1 * h2 * rho ** 3, '13': -sp.Rational(3, 2) * h1, '24': -sp.Rational(3, 2) * h2,
     '14': sp.Rational(3, 2) * h3, '23': sp.Rational(3, 2) * h4, '34': rho ** -3}
D = {'12': 0, '13': -w['13'] / (3 * A1h * A3), '14': -w['14'] / (3 * A1h * A4), '23': -w['23'] / (3 * A2h * A3),
     '24': -w['24'] / (3 * A2h * A4), '34': -rho ** 3 * w['34'] / (3 * A3 * A4)}
T = {'12': 0, '13': A1h * A3, '14': A1h * A4, '23': A2h * A3, '24': A2h * A4, '34': A3 * A4 * rho}
# beta_e = -sum C_{e,l} q_l with q3 = q4 = 0 at b = 0; d14 = gamma, d13 = 2, d24 - d23 = 2 gamma - 1
beta = {'12': sp.zeros(2, 1),
        '13': -((ga - 2) / (8 * A2h)) * q2, '14': -((2 - ga) / (8 * A2h)) * q2,
        '23': -((2 * ga - 1) / (8 * A1h)) * q1, '24': -((1 - 2 * ga) / (8 * A1h)) * q1,
        '34': -(rho * 2 / (8 * A1h)) * q1 - (rho * 2 / (8 * A2h)) * q2}
for e in D:
    print(e, "D =", sp.simplify(D[e]), " T =", sp.simplify(T[e]), " beta =", list(sp.simplify(beta[e])))
c0 = sp.simplify(sum(D[e] * T[e] ** 2 for e in D))
b0 = sp.simplify(sum((D[e] * T[e] * beta[e] for e in D), sp.zeros(2, 1)))
Sbb = sum((D[e] * beta[e] * beta[e].T for e in D), sp.zeros(2, 2))
G0 = sp.simplify(Sbb - b0 * b0.T / c0)
print("c0 =", c0)
print("G0 =", G0)
print("G0 at gamma = 1:", sp.simplify(G0.subs(ga, 1)))
mu = (2 * ga - 1) * (2 - ga) / (3 * ga)
print("tr G0 - 1/4 - mu =", sp.simplify(G0.trace() - sp.Rational(1, 4) - mu))
print("det G0 - mu/4   =", sp.simplify(G0.det() - mu / 4))
print("1/3 - mu =", sp.factor(sp.Rational(1, 3) - mu))
print("alpha-dependence of G0:", sp.simplify(sp.diff(G0, al)))

# (2) numeric comparison with the certified chart formulas at tiny b
import mpmath as mp
import cert
mp.mp.dps = 50
cert.sq = mp.sqrt
from corner_lmax import Ph
for (a0, g0) in ((-0.3, 1.2), (0.2, 0.8), (0.0, 1.0)):
    b = mp.mpf('1e-15')
    x0 = 3 * (a0 * g0 - a0 - g0) / (2 * (1 + g0))
    xs = mp.findroot(lambda t: Ph(b, a0, g0, t), mp.mpf(x0))
    Q = cert.chart(b, mp.mpf(a0), mp.mpf(g0), xs)
    y = cert.opt_y(Q)
    S11, S12, S22 = cert.Smat(Q, y)
    Gs = G0.subs({al: a0, ga: g0})
    print("alpha=%.2f gamma=%.2f: |G(cert) - G0| = %s" % (a0, g0, mp.nstr(max(abs(S11 - float(Gs[0, 0])), abs(S12 - float(Gs[0, 1])), abs(S22 - float(Gs[1, 1]))), 3)))

# (3) membership of points near the maximizers found by corner_lmax.py for alpha near +-1/2
from corner_lmax import point
mp.mp.dps = 40
for b in ('1e-2', '1e-3'):
    for a0, g0 in ((-0.49, 1.0005), (-0.499, 1.0005), (-0.4999, 1.0005), (0.49, 0.9995), (0.499, 0.9995), (0.4999, 0.9995)):
        lm, tr, ok = point(mp.mpf(b), a0, g0)
        Qp = cert.chart(mp.mpf(b), mp.mpf(a0), mp.mpf(g0), mp.findroot(lambda t: Ph(mp.mpf(b), a0, g0, t), mp.mpf(3 * (a0 * g0 - a0 - g0) / (2 * (1 + g0)))))
        print("b=%s alpha=%.4f gamma=%.4f  in region: %s  min g = %s" % (b, a0, g0, ok, mp.nstr(min(Qp['g']), 4)))
