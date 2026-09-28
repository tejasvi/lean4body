"""Formulas for the convex four-body certificate (SPEC.md).

Two layers:

1. `literal(a,b,c,x)`: the SPEC section-1 formulas typed in literally (positions -> distances ->
   signed areas via the generic triangle-area formula -> lam, w_e, D_e, tau_e, beta_e, F).
   Used only by the numerical checks, as the reference.

2. Scaled formulas used by the rigorous verifier.  A "scale" t multiplies the quantities
   that vanish at the collision corner:
        A1 = t*A1h, A2 = t*A2h, R_e - R12 = t*E_e (e = 13,14,23,24), r34 = t*rho,
        y = Y/t, tau_e = t*T_e, w12 = t^5*w12h, w13,w14,w23,w24 = t*w..h, w34 = t^-3*w34h,
        F = t^2 * Ftil.
   Direct coordinates use t = 1 (then the hatted quantities are the plain ones, and the formulas
   are exact algebraic rewrites of SPEC section 1).  The chart a = 1 + b al, c = b ga,
   x = 1/2 + b xi uses t = b; every hatted quantity is then an explicit expression that stays
   finite (and analytic) at b = 0.  `certificate(P)` and `Ftil(P)` are the same code for both.

All functions are generic in the number type (float, mpmath mpf, arb, AD).
"""
from ad import SQRT

EDGES = ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))


def _key(i, j):
    return (i, j) if i < j else (j, i)


def _area(u, v, w):
    """counterclockwise-signed area of the triangle u, v, w"""
    return ((v[0] - u[0]) * (w[1] - u[1]) - (v[1] - u[1]) * (w[0] - u[0])) / 2


# ----------------------------------------------------------------------------------------
# 1. literal SPEC formulas (reference for the checks)
# ----------------------------------------------------------------------------------------
def literal(a, b, c, x):
    s = SQRT(1 - x * x)
    q = {1: (1, 0), 2: (a * x, a * s), 3: (-b, 0), 4: (-c * x, -c * s)}
    R, r, sv = {}, {}, {}
    for (i, j) in EDGES:
        d0 = q[i][0] - q[j][0]
        d1 = q[i][1] - q[j][1]
        R[(i, j)] = d0 * d0 + d1 * d1
        r[(i, j)] = SQRT(R[(i, j)])
        sv[(i, j)] = 1 / r[(i, j)] ** 3
    A = {1: _area(q[2], q[3], q[4]), 2: -_area(q[1], q[3], q[4]),
         3: _area(q[1], q[2], q[4]), 4: -_area(q[1], q[2], q[3])}
    s12, s13, s14, s23, s24, s34 = (sv[e] for e in EDGES)
    lam = (s12 * s34 - s13 * s24) / (s12 + s34 - s13 - s24)
    w = {e: sv[e] - lam for e in EDGES}
    D = {e: -w[e] / (3 * sv[e] * A[e[0]] * A[e[1]]) for e in EDGES}
    tau = {e: A[e[0]] * A[e[1]] * r[e] for e in EDGES}
    beta = {}
    for e in EDGES:
        i, j = e
        bx, by = 0, 0
        for l in (1, 2, 3, 4):
            if l in e:
                continue
            k = [m for m in (1, 2, 3, 4) if m not in (i, j, l)][0]
            N = R[_key(i, k)] + R[_key(j, k)] - R[e]
            coef = r[e] * N / (8 * A[l])
            bx = bx - coef * q[l][0]
            by = by - coef * q[l][1]
        beta[e] = (bx, by)
    r3 = {e: r[e] ** 3 for e in EDGES}
    F = ((r3[(2, 4)] - r3[(1, 4)]) * (r3[(1, 3)] - r3[(1, 2)]) * (r3[(2, 3)] - r3[(3, 4)])
         - (r3[(1, 2)] - r3[(1, 4)]) * (r3[(2, 4)] - r3[(3, 4)]) * (r3[(1, 3)] - r3[(2, 3)]))
    Rg = R
    region = [Rg[(1, 3)] - Rg[(1, 2)], Rg[(2, 4)] - Rg[(1, 2)], Rg[(1, 2)] - Rg[(1, 4)],
              Rg[(1, 2)] - Rg[(2, 3)], Rg[(1, 4)] - Rg[(3, 4)], Rg[(2, 3)] - Rg[(3, 4)]]
    return dict(s=s, q=q, R=R, r=r, sv=sv, A=A, lam=lam, w=w, D=D, tau=tau, beta=beta, F=F,
                region=region)


def trS_literal(L, y):
    tot = 0
    for e in EDGES:
        bx, by = L['beta'][e]
        t = L['tau'][e]
        tot = tot + L['D'][e] * ((bx + y[0] * t) ** 2 + (by + y[1] * t) ** 2)
    return tot


def rigidity_residuals(L):
    """max |beta R - a| and max |tau R| (R = 6x8 rigidity matrix, a = (A1 I, .., A4 I))"""
    q, r, A = L['q'], L['r'], L['A']
    Rm = [[0] * 8 for _ in range(6)]
    for row, (i, j) in enumerate(EDGES):
        ex = (q[i][0] - q[j][0]) / r[(i, j)]
        ey = (q[i][1] - q[j][1]) / r[(i, j)]
        Rm[row][2 * (i - 1)], Rm[row][2 * (i - 1) + 1] = ex, ey
        Rm[row][2 * (j - 1)], Rm[row][2 * (j - 1) + 1] = -ex, -ey
    amat = [[0] * 8 for _ in range(2)]
    for m in (1, 2, 3, 4):
        amat[0][2 * (m - 1)] = A[m]
        amat[1][2 * (m - 1) + 1] = A[m]
    errb = 0
    for d in range(2):
        for col in range(8):
            v = sum(L['beta'][EDGES[row]][d] * Rm[row][col] for row in range(6))
            errb = max(errb, abs(v - amat[d][col]))
    errt = 0
    for col in range(8):
        v = sum(L['tau'][EDGES[row]] * Rm[row][col] for row in range(6))
        errt = max(errt, abs(v))
    return errb, errt


# ----------------------------------------------------------------------------------------
# 2. scaled formulas (verifier)
# ----------------------------------------------------------------------------------------
def prims_direct(a, b, c, x):
    """t = 1.  All squared distances as exact polynomials in (a,b,c,x) (s^2 = 1-x^2 used)."""
    s = SQRT((1 - x) * (1 + x))
    R12 = 1 + a * a - 2 * a * x
    R13 = (1 + b) * (1 + b)
    R14 = 1 + c * c + 2 * c * x
    R23 = a * a + b * b + 2 * a * b * x
    R24 = (a + c) * (a + c)
    rho2 = b * b + c * c - 2 * b * c * x          # = R34 (t = 1)
    E13 = 2 * b + b * b - a * a + 2 * a * x        # R13 - R12
    E14 = -(a + c) * (a - c - 2 * x)               # R14 - R12
    E23 = -(1 + b) * (1 - b - 2 * a * x)           # R23 - R12
    E24 = 2 * a * c + c * c - 1 + 2 * a * x        # R24 - R12
    A3 = (a + c) * s / 2
    A4 = -a * (1 + b) * s / 2
    A1h = b * A3                                   # = A1
    A2h = -(1 + b) * c * s / 2                     # = A2
    q = ((1, 0), (a * x, a * s), (-b, 0), (-c * x, -c * s))
    return dict(t=1, q=q, R12=R12, R13=R13, R14=R14, R23=R23, R24=R24, rho2=rho2,
                E13=E13, E14=E14, E23=E23, E24=E24, A1h=A1h, A2h=A2h, A3=A3, A4=A4)


def prims_chart(b, al, ga, xi):
    """t = b, a = 1 + b al, c = b ga, x = 1/2 + b xi.  Every entry is a polynomial in
    (b, al, ga, xi) except s = sqrt(1 - x^2) (and quantities built from it)."""
    a = 1 + b * al
    c = b * ga
    x = 0.5 + b * xi
    s = SQRT((0.5 - b * xi) * (1.5 + b * xi))
    R12 = 1 + b * (al - 2 * xi) * a
    R13 = (1 + b) * (1 + b)
    R14 = 1 + c * (c + 2 * x)
    R23 = a * a + b * (b + 2 * a * x)
    R24 = (a + c) * (a + c)
    rho2 = 1 - ga * (1 + 2 * b * xi) + ga * ga    # = R34 / b^2
    E13 = 2 + b + a * (2 * xi - al)                # (R13 - R12)/b
    E14 = -(a + c) * (al - ga - 2 * xi)            # (R14 - R12)/b
    E23 = (1 + b) * ((1 + al) + 2 * xi * a)        # (R23 - R12)/b
    E24 = 2 * ga * a + b * ga * ga + al + 2 * xi * a   # (R24 - R12)/b
    A3 = (a + c) * s / 2
    A4 = -a * (1 + b) * s / 2
    A1h = A3                                       # A1 = b*A1h
    A2h = -(1 + b) * ga * s / 2                    # A2 = b*A2h
    q = ((1, 0), (a * x, a * s), (-b, 0), (-c * x, -c * s))
    return dict(t=b, q=q, R12=R12, R13=R13, R14=R14, R23=R23, R24=R24, rho2=rho2,
                E13=E13, E14=E14, E23=E23, E24=E24, A1h=A1h, A2h=A2h, A3=A3, A4=A4)


def _dists(P):
    r12 = SQRT(P['R12'])
    r13 = SQRT(P['R13'])
    r14 = SQRT(P['R14'])
    r23 = SQRT(P['R23'])
    r24 = SQRT(P['R24'])
    rho = SQRT(P['rho2'])
    return r12, r13, r14, r23, r24, rho


def certificate(P):
    """Return (D, T, beta): lists over EDGES of D_e, T_e = tau_e/t, beta_e = (bx, by).
    trS(y) = sum_e D_e |beta_e + (t y) T_e|^2."""
    t = P['t']
    R12, R13, R14, R23, R24, rho2 = P['R12'], P['R13'], P['R14'], P['R23'], P['R24'], P['rho2']
    E13, E14, E23, E24 = P['E13'], P['E14'], P['E23'], P['E24']
    A1h, A2h, A3, A4 = P['A1h'], P['A2h'], P['A3'], P['A4']
    q1, q2, q3, q4 = P['q']
    r12, r13, r14, r23, r24, rho = _dists(P)
    c12 = r12 * R12                                  # r12^3
    s12, s13, s14 = 1 / c12, 1 / (r13 * R13), 1 / (r14 * R14)
    s23, s24 = 1 / (r23 * R23), 1 / (r24 * R24)
    sh34 = 1 / (rho * rho2)                          # t^3 * s34

    # S_e = (s12 - s_e)/t = E_e (r_e^2 + r_e r12 + r12^2) / ((r_e + r12) r_e^3 r12^3)
    def S(E, re, Re):
        return E * (Re + re * r12 + R12) / ((re + r12) * (re * Re) * c12)

    S13, S14, S23, S24 = S(E13, r13, R13), S(E14, r14, R14), S(E23, r23, R23), S(E24, r24, R24)
    t3 = t * t * t
    Dn = sh34 + t3 * (s12 - s13 - s24)               # t^3 (s12 + s34 - s13 - s24)
    w12h = S13 * S24 / Dn                            # w12 = t^5 w12h
    w13h = S13 * (t3 * s13 - sh34) / Dn              # w13 = t w13h
    w24h = S24 * (t3 * s24 - sh34) / Dn              # w24 = t w24h
    t4w = t3 * t * w12h
    w14h = t4w - S14                                 # w14 = t w14h
    w23h = t4w - S23                                 # w23 = t w23h
    w34h = (sh34 - t3 * s13) * (sh34 - t3 * s24) / Dn    # w34 = t^-3 w34h
    D = [-t3 * w12h / (3 * s12 * A1h * A2h),
         -w13h / (3 * s13 * A1h * A3),
         -w14h / (3 * s14 * A1h * A4),
         -w23h / (3 * s23 * A2h * A3),
         -w24h / (3 * s24 * A2h * A4),
         -w34h / (3 * sh34 * A3 * A4)]
    T = [t * A1h * A2h * r12, A1h * A3 * r13, A1h * A4 * r14,
         A2h * A3 * r23, A2h * A4 * r24, A3 * A4 * rho]
    R34 = t * t * rho2
    trho2 = t * rho2
    # N_{e,l} (divided by t where it vanishes at the corner)
    N12_3 = R14 + R24 - R12
    N12_4 = R13 + R23 - R12
    N13_2 = E14 - E13 + trho2
    N13_4 = R12 + R23 - R13
    N14_2 = E13 - E14 + trho2
    N14_3 = R12 + R24 - R14
    N23_1 = E24 - E23 + trho2
    N23_4 = R12 + R13 - R23
    N24_1 = E23 - E24 + trho2
    N24_3 = R12 + R14 - R24
    N34_1 = R23 + R24 - R34
    N34_2 = R13 + R14 - R34

    def comb(r, c1, p1, c2, p2):
        # -r (c1 p1 + c2 p2)
        return (-r * (c1 * p1[0] + c2 * p2[0]), -r * (c1 * p1[1] + c2 * p2[1]))

    beta = [comb(r12, N12_3 / (8 * A3), q3, N12_4 / (8 * A4), q4),
            comb(r13, N13_2 / (8 * A2h), q2, N13_4 / (8 * A4), q4),
            comb(r14, N14_2 / (8 * A2h), q2, N14_3 / (8 * A3), q3),
            comb(r23, N23_1 / (8 * A1h), q1, N23_4 / (8 * A4), q4),
            comb(r24, N24_1 / (8 * A1h), q1, N24_3 / (8 * A3), q3),
            comb(rho, N34_1 / (8 * A1h), q1, N34_2 / (8 * A2h), q2)]
    return D, T, beta


def trS(cert, Y):
    D, T, beta = cert
    tot = 0
    for De, Te, (bx, by) in zip(D, T, beta):
        u = bx + Y[0] * Te
        v = by + Y[1] * Te
        tot = tot + De * (u * u + v * v)
    return tot


def optimal_Y(cert):
    """minimiser of the quadratic Y -> trS(cert, Y) (for float / mpmath evaluation)"""
    D, T, beta = cert
    den = sum(De * Te * Te for De, Te in zip(D, T))
    nx = sum(De * Te * b[0] for De, Te, b in zip(D, T, beta))
    ny = sum(De * Te * b[1] for De, Te, b in zip(D, T, beta))
    return (-nx / den, -ny / den)


def Ftil(P):
    """F / t^2 with r_e^3 - r_f^3 = (R_e - R_f) G(r_e, r_f), G(u,v) = (u^2+uv+v^2)/(u+v).
    For t = 1 this is exactly F."""
    t = P['t']
    R12, R13, R14, R23, R24 = P['R12'], P['R13'], P['R14'], P['R23'], P['R24']
    E13, E14, E23, E24 = P['E13'], P['E14'], P['E23'], P['E24']
    r12, r13, r14, r23, r24, rho = _dists(P)

    def G(u, U, v, V):
        return (U + u * v + V) / (u + v)

    tr34c = t * t * t * rho * P['rho2']              # r34^3
    return ((E24 - E14) * G(r24, R24, r14, R14) * E13 * G(r13, R13, r12, R12)
            * (r23 * R23 - tr34c)
            + E14 * G(r12, R12, r14, R14) * (r24 * R24 - tr34c)
            * (E13 - E23) * G(r13, R13, r23, R23))


def F_direct(a, b, c, x):
    """SPEC's F in direct coordinates, as written (r^3 = R*sqrt(R), R_ij exact polynomials)."""
    R12 = 1 + a * a - 2 * a * x
    R13 = (1 + b) * (1 + b)
    R14 = 1 + c * c + 2 * c * x
    R23 = a * a + b * b + 2 * a * b * x
    R24 = (a + c) * (a + c)
    R34 = b * b + c * c - 2 * b * c * x
    c12, c13, c14 = R12 * SQRT(R12), R13 * SQRT(R13), R14 * SQRT(R14)
    c23, c24, c34 = R23 * SQRT(R23), R24 * SQRT(R24), R34 * SQRT(R34)
    return (c24 - c14) * (c13 - c12) * (c23 - c34) - (c12 - c14) * (c24 - c34) * (c13 - c23)


def region_scaled(P):
    """the six region quantities, divided by the positive factor t for the first four:
    [E13, E24, -E14, -E23, R14 - R34, R23 - R34]; the region is 'all >= 0' (for t > 0)."""
    t = P['t']
    R34 = t * t * P['rho2']
    return [P['E13'], P['E24'], -P['E14'], -P['E23'], P['R14'] - R34, P['R23'] - R34]
