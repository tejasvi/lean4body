"""Sanity checks required before the branch-and-bound (SPEC section 2 and task rules).

 A. sympy: closed forms of R_ij, areas, E_e (direct and chart), R34 = b^2 rho^2, the
    factorisations of w_e, and F = t^2 Ftil.
 B. 50-digit mpmath at random points: sum A_i = 0, sum A_i q_i = 0, sign pattern of A_i,
    beta R = a and tau R = 0 (rigidity matrix), scaled direct formulas (t = 1) vs literal.
 C. chart formulas (t = b) vs the literal direct formulas at random points with b > 0,
    including very small b; finiteness at b = 0 exactly.
 D. points of Z (F = 0 solved at 50 digits): D_e >= 0, min_y trS, rigidity check there.
"""
import random
import sys
import mpmath
from mpmath import mp, mpf
import sympy as sp

import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from formulas import (EDGES, literal, trS_literal, rigidity_residuals, prims_direct,
                      prims_chart, certificate, trS, optimal_Y, Ftil, F_direct, region_scaled)

mp.dps = 50
random.seed(20260925)
FAIL = []


def report(name, err, tol):
    ok = err <= tol
    print(f"  {'OK  ' if ok else 'FAIL'} {name}: max err = {mpmath.nstr(err, 5)} (tol {tol})")
    if not ok:
        FAIL.append(name)


# ------------------------------------------------------------------ A. symbolic
def symbolic_checks():
    print("A. symbolic checks (sympy)")
    a, b, c, x, s = sp.symbols('a b c x s', real=True)
    al, ga, xi = sp.symbols('alpha gamma xi', real=True)
    q = {1: (1, 0), 2: (a * x, a * s), 3: (-b, 0), 4: (-c * x, -c * s)}
    subs_s = {s ** 2: 1 - x ** 2}

    def R(i, j):
        return sp.expand((q[i][0] - q[j][0]) ** 2 + (q[i][1] - q[j][1]) ** 2).subs(subs_s)

    closed = {(1, 2): 1 + a ** 2 - 2 * a * x, (1, 3): (1 + b) ** 2, (1, 4): 1 + c ** 2 + 2 * c * x,
              (2, 3): a ** 2 + b ** 2 + 2 * a * b * x, (2, 4): (a + c) ** 2,
              (3, 4): b ** 2 + c ** 2 - 2 * b * c * x}
    ok = all(sp.expand(R(*e) - closed[e]) == 0 for e in EDGES)
    print(f"  {'OK  ' if ok else 'FAIL'} R_ij closed forms")
    if not ok:
        FAIL.append('R closed forms')

    def area(u, v, w):
        return ((v[0] - u[0]) * (w[1] - u[1]) - (v[1] - u[1]) * (w[0] - u[0])) / 2

    A = {1: area(q[2], q[3], q[4]), 2: -area(q[1], q[3], q[4]),
         3: area(q[1], q[2], q[4]), 4: -area(q[1], q[2], q[3])}
    Acl = {1: b * (a + c) * s / 2, 2: -(1 + b) * c * s / 2, 3: (a + c) * s / 2, 4: -a * (1 + b) * s / 2}
    ok = all(sp.expand(A[i] - Acl[i]) == 0 for i in A)
    ok2 = sp.expand(sum(A.values())) == 0 and all(
        sp.expand(sum(A[i] * q[i][d] for i in A)) == 0 for d in (0, 1))
    print(f"  {'OK  ' if ok else 'FAIL'} signed areas closed forms; "
          f"{'OK  ' if ok2 else 'FAIL'} sum A_i = 0 and sum A_i q_i = 0 identically")
    if not (ok and ok2):
        FAIL.append('areas')

    E = {13: 2 * b + b ** 2 - a ** 2 + 2 * a * x, 14: -(a + c) * (a - c - 2 * x),
         23: -(1 + b) * (1 - b - 2 * a * x), 24: 2 * a * c + c ** 2 - 1 + 2 * a * x}
    ok = (sp.expand(closed[(1, 3)] - closed[(1, 2)] - E[13]) == 0 and
          sp.expand(closed[(1, 4)] - closed[(1, 2)] - E[14]) == 0 and
          sp.expand(closed[(2, 3)] - closed[(1, 2)] - E[23]) == 0 and
          sp.expand(closed[(2, 4)] - closed[(1, 2)] - E[24]) == 0)
    print(f"  {'OK  ' if ok else 'FAIL'} direct E_e = R_e - R12")
    if not ok:
        FAIL.append('direct E')

    # chart
    ch = {a: 1 + b * al, c: b * ga, x: sp.Rational(1, 2) + b * xi}
    A_ = 1 + b * al
    C_ = b * ga
    X_ = sp.Rational(1, 2) + b * xi
    Ech = {13: 2 + b + A_ * (2 * xi - al), 14: -(A_ + C_) * (al - ga - 2 * xi),
           23: (1 + b) * ((1 + al) + 2 * xi * A_), 24: 2 * ga * A_ + b * ga ** 2 + al + 2 * xi * A_}
    Rch = {e: sp.expand(closed[e].subs(ch)) for e in EDGES}
    ok = (sp.expand(Rch[(1, 3)] - Rch[(1, 2)] - b * Ech[13]) == 0 and
          sp.expand(Rch[(1, 4)] - Rch[(1, 2)] - b * Ech[14]) == 0 and
          sp.expand(Rch[(2, 3)] - Rch[(1, 2)] - b * Ech[23]) == 0 and
          sp.expand(Rch[(2, 4)] - Rch[(1, 2)] - b * Ech[24]) == 0)
    rho2 = 1 - ga * (1 + 2 * b * xi) + ga ** 2
    ok_r34 = sp.expand(Rch[(3, 4)] - b ** 2 * rho2) == 0
    R12ch = 1 + b * (al - 2 * xi) * A_
    R14ch = 1 + C_ * (C_ + 2 * X_)
    R23ch = A_ ** 2 + b * (b + 2 * A_ * X_)
    ok_R = (sp.expand(Rch[(1, 2)] - R12ch) == 0 and sp.expand(Rch[(1, 4)] - R14ch) == 0 and
            sp.expand(Rch[(2, 3)] - R23ch) == 0)
    # also: sympy's own quotient of (R_e - R12) by b is a polynomial equal to E_e
    okdiv = True
    for e, Ee in ((13, (1, 3)), (14, (1, 4)), (23, (2, 3)), (24, (2, 4))):
        qt, rm = sp.div(sp.Poly(Rch[Ee] - Rch[(1, 2)], b, al, ga, xi), sp.Poly(b, b, al, ga, xi))
        okdiv &= rm.is_zero and sp.expand(qt.as_expr() - Ech[e]) == 0
    print(f"  {'OK  ' if ok else 'FAIL'} chart E_e = (R_e - R12)/b;  "
          f"{'OK  ' if ok_r34 else 'FAIL'} R34 = b^2 rho^2;  "
          f"{'OK  ' if ok_R else 'FAIL'} chart R12,R14,R23 forms;  "
          f"{'OK  ' if okdiv else 'FAIL'} sympy polynomial division by b (zero remainder)")
    if not (ok and ok_r34 and ok_R and okdiv):
        FAIL.append('chart polynomials')
    # leading-order chart region at b=0
    print("   chart E_e at b=0:", {k: sp.expand(v.subs(b, 0)) for k, v in Ech.items()})
    # A2h, A1h
    ok = (sp.expand(Acl[1].subs(ch) - b * ((A_ + C_) * s / 2)) == 0 and
          sp.expand(Acl[2].subs(ch) - b * (-(1 + b) * ga * s / 2)) == 0)
    print(f"  {'OK  ' if ok else 'FAIL'} chart A1 = b A1h, A2 = b A2h")
    if not ok:
        FAIL.append('chart areas')

    # w_e factorisations
    s12, s13, s14, s23, s24, s34 = sp.symbols('s12 s13 s14 s23 s24 s34', positive=True)
    den = s12 + s34 - s13 - s24
    lam = (s12 * s34 - s13 * s24) / den
    ok = (sp.simplify(s12 - lam - (s12 - s13) * (s12 - s24) / den) == 0 and
          sp.simplify(s13 - lam - (s12 - s13) * (s13 - s34) / den) == 0 and
          sp.simplify(s24 - lam - (s12 - s24) * (s24 - s34) / den) == 0 and
          sp.simplify(s34 - lam - (s34 - s13) * (s34 - s24) / den) == 0 and
          sp.simplify(s14 - lam - ((s14 - s12) + (s12 - s13) * (s12 - s24) / den)) == 0)
    print(f"  {'OK  ' if ok else 'FAIL'} factorisations of w12, w13, w24, w34, w14")
    if not ok:
        FAIL.append('w factorisations')

    # r^3 differences
    u, v = sp.symbols('u v', positive=True)
    ok = sp.simplify(u ** 3 - v ** 3 - (u ** 2 - v ** 2) * (u ** 2 + u * v + v ** 2) / (u + v)) == 0
    print(f"  {'OK  ' if ok else 'FAIL'} r_e^3 - r_f^3 = (R_e - R_f)(r_e^2 + r_e r_f + r_f^2)/(r_e + r_f)")
    if not ok:
        FAIL.append('cube difference')


# ------------------------------------------------------------------ helpers
def rel(u, v):
    return abs(u - v) / max(mpf(1), abs(v))


def rand(lo, hi):
    return mpf(lo) + (mpf(hi) - mpf(lo)) * mpf(random.random())


# ------------------------------------------------------------------ B. direct
def direct_checks(npts=300):
    print("B. direct formulas, 50-digit mpmath, random points (a,b,c,x) in "
          "[0.5,1.75]x[0.125,1]x[0.05,1.75]x[-0.95,0.95]")
    e_sum = e_mom = e_rig_b = e_rig_t = e_R = e_A = e_D = e_T = e_B = e_F = e_F2 = e_tr = mpf(0)
    signs_ok = True
    for _ in range(npts):
        a, b, c, x = rand(0.5, 1.75), rand(0.125, 1), rand(0.05, 1.75), rand(-0.95, 0.95)
        L = literal(a, b, c, x)
        A, q = L['A'], L['q']
        e_sum = max(e_sum, abs(sum(A.values())))
        e_mom = max(e_mom, max(abs(sum(A[i] * q[i][d] for i in A)) for d in (0, 1)))
        signs_ok &= A[1] > 0 and A[3] > 0 and A[2] < 0 and A[4] < 0
        eb, et = rigidity_residuals(L)
        e_rig_b, e_rig_t = max(e_rig_b, eb), max(e_rig_t, et)
        P = prims_direct(a, b, c, x)
        for key, e in (('R12', (1, 2)), ('R13', (1, 3)), ('R14', (1, 4)), ('R23', (2, 3)), ('R24', (2, 4)),
                       ('rho2', (3, 4))):
            e_R = max(e_R, rel(P[key], L['R'][e]))
        for key, i in (('A1h', 1), ('A2h', 2), ('A3', 3), ('A4', 4)):
            e_A = max(e_A, rel(P[key], A[i]))
        D, T, beta = certificate(P)
        for k, e in enumerate(EDGES):
            e_D = max(e_D, rel(D[k], L['D'][e]))
            e_T = max(e_T, rel(T[k], L['tau'][e]))
            e_B = max(e_B, rel(beta[k][0], L['beta'][e][0]), rel(beta[k][1], L['beta'][e][1]))
        e_F = max(e_F, rel(Ftil(P), L['F']))
        e_F2 = max(e_F2, rel(F_direct(a, b, c, x), L['F']))
        y = (rand(-3, 3), rand(-3, 3))
        e_tr = max(e_tr, rel(trS((D, T, beta), y), trS_literal(L, y)))
    tol = mpf(10) ** -40
    report("sum A_i = 0", e_sum, tol)
    report("sum A_i q_i = 0", e_mom, tol)
    print(f"  {'OK  ' if signs_ok else 'FAIL'} A1, A3 > 0 > A2, A4 at all points")
    if not signs_ok:
        FAIL.append('area signs')
    report("beta R = a  (rigidity matrix)", e_rig_b, tol)
    report("tau R = 0   (self-stress)", e_rig_t, tol)
    report("closed-form R_ij vs coordinates", e_R, tol)
    report("closed-form A_i vs triangle areas", e_A, tol)
    report("scaled (t=1) D_e vs literal", e_D, tol)
    report("scaled (t=1) tau_e vs literal", e_T, tol)
    report("scaled (t=1) beta_e vs literal", e_B, tol)
    report("factored F (t=1) vs literal F", e_F, tol)
    report("F_direct vs literal F", e_F2, tol)
    report("trS(y) scaled vs literal", e_tr, tol)


# ------------------------------------------------------------------ C. chart
def chart_checks(npts=300):
    print("C. chart formulas vs literal direct ones, 50-digit mpmath, random (b,al,ga,xi), b>0")
    for (blo, bhi, label, dps) in (('1e-6', '0.125', "b in [1e-6, 1/8], 50 digits", 50),
                                   ('1e-30', '1e-20', "b in [1e-30, 1e-20], 200 digits", 200)):
        # at b ~ 1e-25 the LITERAL formulas lose ~ 5*25 + 3*25 digits to cancellation (w12 ~ b^5 is
        # computed there as s12 - lam with s34 ~ b^-3), so that range is run at 200 digits.
        mp.dps = dps
        blo, bhi = mpf(blo), mpf(bhi)
        e_D = e_T = e_B = e_F = e_reg = e_tr = mpf(0)
        n = 0
        while n < npts:
            b = blo + (bhi - blo) * mpf(random.random())
            al, ga, xi = rand(-2, 1), rand(0.2, 3), rand(-2.5, 0.7)
            a, c, x = 1 + b * al, b * ga, mpf(1) / 2 + b * xi
            if not (-1 < x < 1):
                continue
            n += 1
            L = literal(a, b, c, x)
            P = prims_chart(b, al, ga, xi)
            D, T, beta = certificate(P)
            for k, e in enumerate(EDGES):
                e_D = max(e_D, rel(D[k], L['D'][e]))
                e_T = max(e_T, rel(T[k], L['tau'][e] / b))
                e_B = max(e_B, rel(beta[k][0], L['beta'][e][0]), rel(beta[k][1], L['beta'][e][1]))
            e_F = max(e_F, rel(Ftil(P), L['F'] / b ** 2))
            rg = region_scaled(P)
            lr = L['region']
            e_reg = max(e_reg, max(rel(rg[i], lr[i] / b) for i in range(4)),
                        max(rel(rg[i], lr[i]) for i in (4, 5)))
            Y = (rand(-2, 2), rand(-2, 2))
            e_tr = max(e_tr, rel(trS((D, T, beta), Y), trS_literal(L, (Y[0] / b, Y[1] / b))))
        tol = mpf(10) ** -35
        print(f"  [{label}]")
        report("D_e chart vs literal", e_D, tol)
        report("tau_e/b chart vs literal", e_T, tol)
        report("beta_e chart vs literal", e_B, tol)
        report("Ftil vs F/b^2", e_F, tol)
        report("scaled region quantities vs literal", e_reg, tol)
        report("trS(Y/b) chart vs literal", e_tr, tol)
    # b = 0 exactly: finite values, and continuity b -> 0
    mp.dps = 50
    print("  b = 0 exactly (limit formulas) vs b = 1e-40, 50 digits:")
    worst = mpf(0)
    for _ in range(50):
        al, ga, xi = rand(-1, 0.4), rand(0.6, 1.8), rand(-1.2, -0.4)
        c0 = certificate(prims_chart(mpf(0), al, ga, xi))
        c1 = certificate(prims_chart(mpf('1e-40'), al, ga, xi))
        for k in range(6):
            worst = max(worst, rel(c0[0][k], c1[0][k]), rel(c0[1][k], c1[1][k]),
                        rel(c0[2][k][0], c1[2][k][0]), rel(c0[2][k][1], c1[2][k][1]))
        worst = max(worst, rel(Ftil(prims_chart(mpf(0), al, ga, xi)),
                               Ftil(prims_chart(mpf('1e-40'), al, ga, xi))))
    report("continuity at b=0 of D_e, T_e, beta_e, Ftil", worst, mpf(10) ** -35)


# ------------------------------------------------------------------ D. points of Z
def z_points(ntry=4000):
    print("D. points of Z (F = 0 solved at 50 digits, region inequalities checked)")
    found = []
    worst_rig = mpf(0)
    minD = mpf(10)
    maxtr = mpf(0)
    maxtr_pt = None
    for _ in range(ntry):
        a, b, c = rand(0.5, 1.75), rand(0.125, 1), rand(0.03, 1.75)
        # bracket roots of F in x on a grid, restricted to the region
        fa, fb, fc = float(a), float(b), float(c)
        xs = [-1 + 2 * k / 400 for k in range(1, 400)]
        vals = []
        for xx in xs:
            reg = region_scaled(prims_direct(fa, fb, fc, xx))
            vals.append((xx, F_direct(fa, fb, fc, xx), min(reg)))
        for (x0, f0, g0), (x1, f1, g1) in zip(vals, vals[1:]):
            if f0 * f1 < 0 and max(g0, g1) > -0.05:
                xr = mpmath.findroot(lambda t: F_direct(a, b, c, t), (mpf(x0), mpf(x1)), solver='anderson')
                L = literal(a, b, c, xr)
                if min(L['region']) < 0:
                    continue
                eb, et = rigidity_residuals(L)
                worst_rig = max(worst_rig, eb, et)
                P = prims_direct(a, b, c, xr)
                cert = certificate(P)
                minD = min(minD, min(cert[0]))
                Y = optimal_Y(cert)
                v = trS(cert, Y)
                found.append((a, b, c, xr, v))
                if v > maxtr:
                    maxtr, maxtr_pt = v, (a, b, c, xr)
    print(f"  {len(found)} points of Z found in direct box; beta R = a, tau R = 0 residual "
          f"{mpmath.nstr(worst_rig, 3)}")
    print(f"  min D_e over them = {mpmath.nstr(minD, 6)} (expected >= 0)")
    print(f"  max_Z min_y trS = {mpmath.nstr(maxtr, 8)} at (a,b,c,x) = "
          f"{[mpmath.nstr(t, 6) for t in maxtr_pt] if maxtr_pt else None}")
    if minD < 0:
        FAIL.append('D_e negative on Z')
    return found


if __name__ == '__main__':
    symbolic_checks()
    direct_checks()
    chart_checks()
    mp.dps = 50
    z_points(int(sys.argv[1]) if len(sys.argv) > 1 else 600)
    print("\nFAILED CHECKS:", FAIL if FAIL else "none")
