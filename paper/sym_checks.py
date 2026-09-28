"""Symbolic checks of identities used in the paper (run from paper):
  1. F * prod(s_e) = -P, with F the CCR/Dziobek function in r^3 form and P = E14 E23 delta + D13 D24 (E14 + E23 - E34).
  2. the w_e formulas: with lambda = (s12 s34 - s13 s24)/delta, s_e - lambda equals the stated expressions, and
     w12 w34 - w14 w23 = -P/delta.
  3. CCR coordinates: the factorizations R12-R14 = (a+c) g3 etc. and the signed areas (sum A = 0, sum A q = 0).
  4. Heron: d A / d r_1 = r_1 (r_2^2 + r_3^2 - r_1^2) / (8 A).
  5. the blow-up chart formulas of the Appendix and the b = 0 limit of the rescaled Dziobek function."""
import sympy as sp

s12, s13, s14, s23, s24, s34 = sp.symbols('s12 s13 s14 s23 s24 s34', positive=True)
t = {k: 1 / v for k, v in zip(('12', '13', '14', '23', '24', '34'), (s12, s13, s14, s23, s24, s34))}   # r^3 = 1/s
F = (t['24'] - t['14']) * (t['13'] - t['12']) * (t['23'] - t['34']) - (t['12'] - t['14']) * (t['24'] - t['34']) * (t['13'] - t['23'])
D13, D24 = s12 - s13, s12 - s24
E14, E23, E34 = s14 - s12, s23 - s12, s34 - s12
delta = D13 + D24 + E34
P = E14 * E23 * delta + D13 * D24 * (E14 + E23 - E34)
print("1. F*prod s + P =", sp.simplify(sp.expand(F * s12 * s13 * s14 * s23 * s24 * s34 + P)))
lam = (s12 * s34 - s13 * s24) / delta
w = {'12': D13 * D24 / delta, '34': (D13 + E34) * (D24 + E34) / delta, '13': -D13 * (D13 + E34) / delta,
     '24': -D24 * (D24 + E34) / delta, '14': E14 + D13 * D24 / delta, '23': E23 + D13 * D24 / delta}
sv = {'12': s12, '13': s13, '14': s14, '23': s23, '24': s24, '34': s34}
print("2. s_e - lambda - w_e:", [sp.simplify(sv[k] - lam - w[k]) for k in w],
      " w12 w34 - w13 w24:", sp.simplify(w['12'] * w['34'] - w['13'] * w['24']),
      " w12 w34 - w14 w23 + P/delta:", sp.simplify(w['12'] * w['34'] - w['14'] * w['23'] + P / delta))

a, b, c, x = sp.symbols('a b c x', real=True)
s = sp.sqrt(1 - x ** 2)
q = [sp.Matrix([1, 0]), sp.Matrix([a * x, a * s]), sp.Matrix([-b, 0]), sp.Matrix([-c * x, -c * s])]
R = lambda i, j: sp.expand(((q[i - 1] - q[j - 1]).T * (q[i - 1] - q[j - 1]))[0])
g1 = b * b + 2 * b - a * a + 2 * a * x; g2 = c * c + 2 * a * c - 1 + 2 * a * x
g3 = a - c - 2 * x; g4 = 1 - b - 2 * a * x; g5 = a - c + 2 * b * x; g6 = 1 - b + 2 * c * x
print("3. factorizations:", [sp.simplify(e) for e in (R(1, 3) - R(1, 2) - g1, R(2, 4) - R(1, 2) - g2, R(1, 2) - R(1, 4) - (a + c) * g3,
      R(1, 2) - R(2, 3) - (1 + b) * g4, R(2, 3) - R(3, 4) - (a + c) * g5, R(1, 4) - R(3, 4) - (1 + b) * g6)])
tri = lambda i, j, k: sp.Rational(1, 2) * sp.Matrix.hstack(q[j - 1] - q[i - 1], q[k - 1] - q[i - 1]).det()
A = [tri(2, 3, 4), -tri(1, 3, 4), tri(1, 2, 4), -tri(1, 2, 3)]
print("   areas:", [sp.factor(sp.simplify(v)) for v in A], " sum:", sp.simplify(sum(A)),
      " sum A q:", sp.simplify(sum((A[i] * q[i] for i in range(4)), sp.zeros(2, 1))).T)
r1, r2, r3 = sp.symbols('r1 r2 r3', positive=True)
Ar = sp.sqrt(2 * (r1**2 * r2**2 + r2**2 * r3**2 + r3**2 * r1**2) - r1**4 - r2**4 - r3**4) / 4
print("4. Heron derivative:", sp.simplify(sp.diff(Ar, r1) - r1 * (r2**2 + r3**2 - r1**2) / (8 * Ar)))
# 5. blow-up chart a = 1 + b al, c = b ga, x = 1/2 + b xi (Section 4.4, the Appendix): h_i = g_i / b (i <= 4),
#    h5 = g5, h6 = g6; R_e = 1 + b d_e for the five long edges and R34 = b^2 rho^2; and the b = 0 limit of the
#    rescaled Dziobek function Ph = b P, which is -(9/4) rho^-3 (3 al + 3 ga - 3 al ga + 2 xi (1 + ga)).
bb, al, ga, xi = sp.symbols('b alpha gamma xi', real=True)
ch = {a: 1 + bb * al, b: bb, c: bb * ga, x: sp.Rational(1, 2) + bb * xi}
xx = ch[x]
h = [2 - al + 2 * xi + bb * (1 - al**2 + 2 * al * xi), 2 * ga + 2 * xi + al + bb * (ga**2 + 2 * al * ga + 2 * al * xi),
     al - ga - 2 * xi, -(1 + 2 * xi + al + 2 * bb * al * xi), 1 + bb * (al - ga + 1 + 2 * bb * xi),
     1 - bb + bb * ga + 2 * bb**2 * ga * xi]
gs = [g1, g2, g3, g4, g5, g6]
print("5. chart h_i:", [sp.expand(gs[i].subs(ch) / (bb if i < 4 else 1) - h[i]) for i in range(6)])
d = {(1, 2): al - 2 * xi + bb * (al**2 - 2 * al * xi), (1, 3): 2 + bb, (1, 4): ga * (bb * ga + 2 * xx),
     (2, 3): 2 * al + bb * al**2 + 2 * xx * (1 + bb * al) + bb, (2, 4): 2 * (al + ga) + bb * (al + ga)**2}
rho2 = 1 - 2 * ga * xx + ga**2
print("   R_e - 1 - b d_e:", [sp.expand(R(*e).subs(ch) - 1 - bb * d[e]) for e in d],
      " R34 - b^2 rho^2:", sp.expand(R(3, 4).subs(ch) - bb**2 * rho2))
hd0 = sp.Rational(3, 2)                                   # hd(1, 1): all long R_e = 1 at b = 0
h0 = [v.subs(bb, 0) for v in h]
rho = sp.sqrt(rho2.subs(bb, 0))
Ph0 = (h0[2] * hd0) * (h0[3] * hd0) * rho**-3 + (h0[0] * hd0) * (h0[1] * hd0) * (0 - rho**-3)
print("   Ph at b = 0 + (9/4) rho^-3 (3al + 3ga - 3al ga + 2xi(1+ga)):",
      sp.simplify(Ph0 + sp.Rational(9, 4) * rho**-3 * (3 * al + 3 * ga - 3 * al * ga + 2 * xi * (1 + ga))))
