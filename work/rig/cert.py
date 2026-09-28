"""Mass-free nondegeneracy certificate for convex four-body central configurations, in two charts.
Generic in the number type: works with float, mpmath.mpf, or flint.arb (uses only + - * / sqrt abs).

Direct chart (CCR coordinates, x = cos(theta)):
    q1 = (1, 0), q2 = a (x, s), q3 = (-b, 0), q4 = -c (x, s),  s = sqrt(1 - x^2).
Blow-up chart at P1 = (1, 0, 0, x = 1/2):  a = 1 + b al, c = b ga, x = 1/2 + b xi.

For a configuration, returns
    g   : the six CCR inequalities (a point is in closure(C) iff all g >= 0), as used for discarding boxes
    pos, neg : the Dziobek condition F = 0 is equivalent (in C) to pos = neg  (F * prod s_e = -(pos - neg))
    D, T, B : compliances D_e, self-stress tau_e (T), particular solution beta_e (B, 2-vectors) of  beta R = a
The certificate is  tr S(y) = sum_e D_e |beta_e + y tau_e|^2 < THRESH  (3/4 in the paper), with one vector y for a
whole box.
In the blow-up chart T and B are the rescaled tau/b and beta; pos, neg are rescaled by a positive power of b;
the g's are divided by b; S(y) is unchanged (D_e tau_e^2 and D_e beta_e^2 are invariant under the rescaling)."""


def sq(v):
    return v.sqrt() if hasattr(v, 'sqrt') else v ** 0.5


def hd(X, Y):
    """-(X^{-3/2} - Y^{-3/2}) / (X - Y) > 0"""
    sx, sy = sq(X), sq(Y)
    return (X + sx * sy + Y) / ((sx + sy) * X * sx * Y * sy)


EDGES = ('12', '13', '14', '23', '24', '34')


def direct(a, b, c, x):
    one = a * 0 + 1
    s = sq(one - x * x)
    g1 = b * b + 2 * b - a * a + 2 * a * x
    g2 = c * c + 2 * a * c - 1 + 2 * a * x
    g3 = a - c - 2 * x
    g4 = 1 - b - 2 * a * x
    g5 = a - c + 2 * b * x
    g6 = 1 - b + 2 * c * x
    R = {'12': 1 + a * a - 2 * a * x, '13': (1 + b) ** 2, '14': 1 + c * c + 2 * c * x,
         '23': a * a + b * b + 2 * a * b * x, '24': (a + c) ** 2, '34': b * b + c * c - 2 * b * c * x}
    r = {k: sq(v) for k, v in R.items()}
    sv = {k: 1 / (R[k] * r[k]) for k in R}
    d12_14 = (a + c) * g3
    d12_23 = (1 + b) * g4
    d12_34 = d12_23 + (a + c) * g5
    Dl13 = g1 * hd(R['12'], R['13'])            # s12 - s13
    Dl24 = g2 * hd(R['12'], R['24'])            # s12 - s24
    E34 = d12_34 * hd(R['34'], R['12'])         # s34 - s12
    E14 = d12_14 * hd(R['14'], R['12'])         # s14 - s12
    E23 = d12_23 * hd(R['23'], R['12'])         # s23 - s12
    den = Dl13 + E34 + Dl24
    w = {}
    w['12'] = Dl13 * Dl24 / den
    w['34'] = (E34 + Dl13) * (E34 + Dl24) / den
    w['13'] = -Dl13 * (E34 + Dl13) / den
    w['24'] = -Dl24 * (E34 + Dl24) / den
    w['14'] = E14 + w['12']
    w['23'] = E23 + w['12']
    pos = E14 * E23 * den + Dl13 * Dl24 * (E14 + E23)
    neg = Dl13 * Dl24 * E34
    A = (b * (a + c) * s / 2, -c * (1 + b) * s / 2, (a + c) * s / 2, -a * (1 + b) * s / 2)
    q = ((one, 0 * one), (a * x, a * s), (-b * one, 0 * one), (-c * x, -c * s))
    D, T, B = {}, {}, {}
    for e in EDGES:
        i, j = int(e[0]) - 1, int(e[1]) - 1
        AA = A[i] * A[j]
        D[e] = -w[e] / (3 * sv[e] * AA)          # = |w_e| / (3 s_e |A_i A_j|) on closure(C)
        T[e] = AA * r[e]
    B = _betas(R, r, A, q)
    return dict(g=(g1, g2, g3, g4, g5, g6), pos=pos, neg=neg, D=D, T=T, B=B, w=w, A=A, s=s)


# N_{e,l}: for edge e = ij and l not in e, triangle = vertices other than l, third vertex k; N = R_ik + R_jk - R_ij
def _N(R, e, l):
    i, j = e[0], e[1]
    k = [v for v in '1234' if v not in (i, j, l)][0]
    return R[''.join(sorted(i + k))] + R[''.join(sorted(j + k))] - R[e]


def _betas(R, r, A, q, Nfun=None):
    B = {}
    for e in EDGES:
        bx = by = 0
        for l in '1234':
            if l in e:
                continue
            N = _N(R, e, l) if Nfun is None else Nfun(e, l)
            coef = r[e] * N / (8 * A[int(l) - 1])
            ql = q[int(l) - 1]
            bx = bx - coef * ql[0]
            by = by - coef * ql[1]
        B[e] = (bx, by)
    return B


def chart(b, al, ga, xi):
    """blow-up chart; valid down to b = 0 (inclusive)"""
    one = b * 0 + 1
    a = one + b * al
    c = b * ga
    x = one / 2 + b * xi
    s = sq(one - x * x)
    h1 = 2 - al + 2 * xi + b * (1 - al * al + 2 * al * xi)          # g1 / b
    h2 = 2 * ga + 2 * xi + al + b * (ga * ga + 2 * al * ga + 2 * al * xi)   # g2 / b
    h3 = al - ga - 2 * xi                                            # g3 / b
    h4 = -(1 + 2 * xi + al + 2 * b * al * xi)                        # g4 / b
    h5 = 1 + b * (al - ga + 1 + 2 * b * xi)                          # g5
    h6 = 1 - b + b * ga + 2 * b * b * ga * xi                        # g6
    # R_e = 1 + b d_e for the five long edges; r34^2 = b^2 rho2
    d = {'12': al - 2 * xi + b * (al * al - 2 * al * xi), '13': 2 + b, '14': ga * (b * ga + 2 * x),
         '23': 2 * al + b * al * al + 2 * x * (1 + b * al) + b, '24': 2 * (al + ga) + b * (al + ga) ** 2}
    R = {k: 1 + b * v for k, v in d.items()}
    rho2 = 1 - 2 * ga * x + ga * ga
    rho = sq(rho2)
    r = {k: sq(v) for k, v in R.items()}
    sv = {k: 1 / (R[k] * r[k]) for k in R}
    # (s_e - s_f)/b for long edges:  (R_f - R_e)/b * hd(R_e, R_f) = (d_f - d_e) hd
    Dl13 = h1 * hd(R['12'], R['13'])            # (s12 - s13)/b ;  (R13 - R12)/b = h1
    Dl24 = h2 * hd(R['12'], R['24'])            # (s12 - s24)/b
    E14 = (a + c) * h3 * hd(R['14'], R['12'])   # (s14 - s12)/b
    E23 = (1 + b) * h4 * hd(R['23'], R['12'])     # (s23 - s12)/b ; (R12 - R23)/b = (1+b) g4/b = (1+b) h4
    # s34 = b^-3 / rho^3 ;  b^3 (s34 - s12) = 1/rho^3 - b^3 s12
    E34h = 1 / (rho2 * rho) - b ** 3 * sv['12']  # b^3 (s34 - s12)
    denh = E34h + b ** 4 * (Dl13 + Dl24)        # b^3 den
    w12h = Dl13 * Dl24 / denh                   # w12 = b^5 w12h
    wh = {}                                     # w_e = b wh_e for 13, 24, 14, 23 ;  w34 = b^-3 wh34
    wh['13'] = -Dl13 * (E34h + b ** 4 * Dl13) / denh
    wh['24'] = -Dl24 * (E34h + b ** 4 * Dl24) / denh
    wh['14'] = E14 + b ** 4 * w12h
    wh['23'] = E23 + b ** 4 * w12h
    wh34 = (E34h + b ** 4 * Dl13) * (E34h + b ** 4 * Dl24) / denh
    # Dziobek: pos - neg with every term rescaled by b^-3 * b^-3 ... : use the b-free combination
    #   pos = E14 E23 den + Dl13 Dl24 (E14 + E23),  neg = Dl13 Dl24 E34   (unscaled)
    #   with E14 = b E14h..., den = b^-3 denh, E34 = b^-3 E34h:  pos = b^-1 E14h E23h denh + b^3 (...), neg = b^-1 Dl Dl E34h
    posh = E14 * E23 * denh + b ** 4 * Dl13 * Dl24 * (E14 + E23)
    negh = Dl13 * Dl24 * E34h
    A1h = s * (a + c) / 2
    A2h = -ga * s * (1 + b) / 2
    A3 = s * (a + c) / 2
    A4 = -a * s * (1 + b) / 2
    D = {}
    D['12'] = -b ** 3 * w12h / (3 * sv['12'] * A1h * A2h)
    D['13'] = -wh['13'] / (3 * sv['13'] * A1h * A3)
    D['14'] = -wh['14'] / (3 * sv['14'] * A1h * A4)
    D['23'] = -wh['23'] / (3 * sv['23'] * A2h * A3)
    D['24'] = -wh['24'] / (3 * sv['24'] * A2h * A4)
    D['34'] = -rho2 * rho * wh34 / (3 * A3 * A4)
    T = {'12': b * A1h * A2h * r['12'], '13': A1h * A3 * r['13'], '14': A1h * A4 * r['14'],
         '23': A2h * A3 * r['23'], '24': A2h * A4 * r['24'], '34': A3 * A4 * rho}
    q = ((one, 0 * one), (a * x, a * s), (-b * one, 0 * one), (-c * x, -c * s))
    R34 = b * b * rho2

    # beta_e = - sum_{l not in e} r_e N_{e,l} / (8 A_l) q_l, with A1 = b A1h, A2 = b A2h and N's rescaled accordingly
    def coef(e, l):
        if e == '12':
            return {'3': r['12'] * (R['14'] + R['24'] - R['12']) / (8 * A3),
                    '4': r['12'] * (R['13'] + R['23'] - R['12']) / (8 * A4)}[l]
        if e == '13':
            return {'2': r['13'] * (d['14'] - d['13'] + b * rho2) / (8 * A2h),
                    '4': r['13'] * (R['12'] + R['23'] - R['13']) / (8 * A4)}[l]
        if e == '14':
            return {'2': r['14'] * (d['13'] - d['14'] + b * rho2) / (8 * A2h),
                    '3': r['14'] * (R['12'] + R['24'] - R['14']) / (8 * A3)}[l]
        if e == '23':
            return {'1': r['23'] * (d['24'] - d['23'] + b * rho2) / (8 * A1h),
                    '4': r['23'] * (R['12'] + R['13'] - R['23']) / (8 * A4)}[l]
        if e == '24':
            return {'1': r['24'] * (d['23'] - d['24'] + b * rho2) / (8 * A1h),
                    '3': r['24'] * (R['12'] + R['14'] - R['24']) / (8 * A3)}[l]
        if e == '34':
            return {'1': rho * (R['23'] + R['24'] - R34) / (8 * A1h),
                    '2': rho * (R['13'] + R['14'] - R34) / (8 * A2h)}[l]
    B = {}
    for e in EDGES:
        bx = by = 0
        for l in '1234':
            if l in e:
                continue
            cf = coef(e, l)
            ql = q[int(l) - 1]
            bx = bx - cf * ql[0]
            by = by - cf * ql[1]
        B[e] = (bx, by)
    return dict(g=(h1, h2, h3, h4, h5, h6), pos=posh, neg=negh, D=D, T=T, B=B, s=s)


def opt_y(Q):
    D, T, B = Q['D'], Q['T'], Q['B']
    Tt = sum(D[k] * T[k] * T[k] for k in EDGES)
    return (-sum(D[k] * T[k] * B[k][0] for k in EDGES) / Tt, -sum(D[k] * T[k] * B[k][1] for k in EDGES) / Tt)


def Smat(Q, y):
    D, T, B = Q['D'], Q['T'], Q['B']
    S11 = S12 = S22 = 0
    for k in EDGES:
        ux = B[k][0] + y[0] * T[k]
        uy = B[k][1] + y[1] * T[k]
        S11 = S11 + D[k] * ux * ux
        S12 = S12 + D[k] * ux * uy
        S22 = S22 + D[k] * uy * uy
    return S11, S12, S22


def trS(Q, y):
    D, T, B = Q['D'], Q['T'], Q['B']
    t = 0
    for k in EDGES:
        ux = B[k][0] + y[0] * T[k]
        uy = B[k][1] + y[1] * T[k]
        t = t + D[k] * (ux * ux + uy * uy)
    return t


def direct_P(a, b, c, x):
    """(g1..g6, P) with P = pos - neg;  F = 0  <=>  P = 0"""
    g1 = b * b + 2 * b - a * a + 2 * a * x
    g2 = c * c + 2 * a * c - 1 + 2 * a * x
    g3 = a - c - 2 * x
    g4 = 1 - b - 2 * a * x
    g5 = a - c + 2 * b * x
    g6 = 1 - b + 2 * c * x
    R12 = 1 + a * a - 2 * a * x
    R13 = (1 + b) ** 2
    R14 = 1 + c * c + 2 * c * x
    R23 = a * a + b * b + 2 * a * b * x
    R24 = (a + c) ** 2
    R34 = b * b + c * c - 2 * b * c * x
    d12_23 = (1 + b) * g4
    Dl13 = g1 * hd(R12, R13)
    Dl24 = g2 * hd(R12, R24)
    E34 = (d12_23 + (a + c) * g5) * hd(R34, R12)
    E14 = (a + c) * g3 * hd(R14, R12)
    E23 = d12_23 * hd(R23, R12)
    P = E14 * E23 * (Dl13 + E34 + Dl24) + Dl13 * Dl24 * (E14 + E23 - E34)
    return (g1, g2, g3, g4, g5, g6), P


def chart_P(b, al, ga, xi):
    """(h1..h6, Ph) with Ph = posh - negh  (a positive multiple of P)"""
    x = 0.5 + b * xi
    a = 1 + b * al
    c = b * ga
    h1 = 2 - al + 2 * xi + b * (1 - al * al + 2 * al * xi)
    h2 = 2 * ga + 2 * xi + al + b * (ga * ga + 2 * al * ga + 2 * al * xi)
    h3 = al - ga - 2 * xi
    h4 = -(1 + 2 * xi + al + 2 * b * al * xi)
    h5 = 1 + b * (al - ga + 1 + 2 * b * xi)
    h6 = 1 - b + b * ga + 2 * b * b * ga * xi
    R12 = 1 + b * (al - 2 * xi + b * (al * al - 2 * al * xi))
    R13 = 1 + b * (2 + b)
    R14 = 1 + b * (ga * (b * ga + 2 * x))
    R23 = 1 + b * (2 * al + b * al * al + 2 * x * (1 + b * al) + b)
    R24 = 1 + b * (2 * (al + ga) + b * (al + ga) ** 2)
    rho2 = 1 - 2 * ga * x + ga * ga
    rho = sq(rho2)
    Dl13 = h1 * hd(R12, R13)
    Dl24 = h2 * hd(R12, R24)
    E14 = (a + c) * h3 * hd(R14, R12)
    E23 = (1 + b) * h4 * hd(R23, R12)
    E34h = 1 / (rho2 * rho) - b ** 3 / (R12 * sq(R12))
    denh = E34h + b ** 4 * (Dl13 + Dl24)
    Ph = E14 * E23 * denh + Dl13 * Dl24 * (b ** 4 * (E14 + E23) - E34h)
    return (h1, h2, h3, h4, h5, h6), Ph
