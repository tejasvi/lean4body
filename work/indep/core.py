"""Rigorous per-box routines (python-flint arb).  Used by verify.py.

A box is a list of three (lo, hi) pairs of Python floats with dyadic values (exactly
representable, converted to arb exactly).  Mode 'direct': variables (a, b, c; x).
Mode 'chart': variables (b, al, ga; xi).  In both modes the dependent coordinate is the 4th.
"""
import math
from flint import arb, ctx

from ad import AD
from formulas import (prims_direct, prims_chart, certificate, trS, optimal_Y, Ftil, F_direct,
                      region_scaled)

PREC = 96
ctx.prec = PREC
THREE_QUARTERS = arb(3) / 4          # exact (dyadic)


def ball(lo, hi):
    """arb ball containing [lo, hi] (lo, hi exact floats)"""
    return arb.union(arb(lo), arb(hi))


def midpoint(lo, hi):
    m = (arb(lo) + arb(hi)) / 2
    assert m.is_exact()
    return m


def certainly_neg(v):
    return v.is_finite() and v < 0


def certainly_pos(v):
    return v.is_finite() and v > 0


# ------------------------------------------------------------------ function evaluation
def prims(mode, p1, p2, p3, p4):
    if mode == 'direct':
        return prims_direct(p1, p2, p3, p4)
    return prims_chart(p1, p2, p3, p4)


def Fval(mode, p1, p2, p3, p4, direct_F='literal'):
    """the CC equation used for the root enclosure: SPEC's F (direct) or F/b^2 (chart)"""
    if mode == 'direct' and direct_F == 'literal':
        return F_direct(p1, p2, p3, p4)
    return Ftil(prims(mode, p1, p2, p3, p4))


# ------------------------------------------------------------------ region: x-range and discards
def linear_xrange(mode, B):
    """Rigorous [xlo, xhi] (arb points) containing every x such that some point of B x {x}
    satisfies the six region inequalities (each is linear in x).  Returns None if it is
    certainly empty, or if one of the x-free 'pair' conditions is certainly violated on B."""
    p = [ball(lo, hi) for lo, hi in B]
    if mode == 'direct':
        a, b, c = p
        # region <=> 2ax >= a^2-2b-b^2, 2ax >= 1-2ac-c^2, 2x <= a-c, 2ax <= 1-b,
        #            2cx >= b-1, 2bx >= c-a        (a > 0, b in (0,1], c >= 0)
        pairs = [2 * b + b * b - a * c, 1 + b + b * b - a * a, a * a + a * c + c * c - 1,
                 2 * a * c + c * c - b, c * (a - c) + 1 - b, a - c, b * (1 - b) - a * (c - a)]
        lows = [(a * a - 2 * b - b * b) / (2 * a), (1 - 2 * a * c - c * c) / (2 * a),
                (c - a) / (2 * b), arb(-1)]
        if c.lower() > 0:
            lows.append(-(1 - b) / (2 * c))
        highs = [(a - c) / 2, (1 - b) / (2 * a), arb(1)]
    else:
        b, al, ga = p
        u = 1 + b * al                     # = a > 0
        # region (b > 0) <=> 2 xi u >= al u - 2 - b, 2 xi u >= -(2 ga u + b ga^2 + al),
        #                    2 xi <= al - ga, 2 xi u <= -(1 + al)     (+ two non-binding ones)
        pairs = [2 + b - ga * u, 1 + b - al * (2 + b * al), (al + ga) * u + b * ga * ga + al,
                 2 * ga * u + b * ga * ga - 1]
        lows = [al / 2 - (2 + b) / (2 * u), -ga - (b * ga * ga + al) / (2 * u)]
        highs = [(al - ga) / 2, -(1 + al) / (2 * u)]
    for P in pairs:
        if certainly_neg(P):
            return None
    xlo = None
    for L in lows:
        if not L.is_finite():
            continue
        v = L.lower()
        xlo = v if xlo is None or v > xlo else xlo
    xhi = None
    for U in highs:
        if not U.is_finite():
            continue
        v = U.upper()
        xhi = v if xhi is None or v < xhi else xhi
    if xlo is None or xhi is None:
        # cannot happen on the domains used (u = 1 + b al >= 3/4, b >= 1/8 in direct mode);
        # refuse rather than guess a range
        raise RuntimeError(f"linear x-range not finite on box {B}")
    if xlo > xhi:
        return None
    return xlo, xhi


def region_violated(mode, Bb, X):
    """True if some region inequality is certainly violated on all of B x X"""
    for g in region_scaled(prims(mode, Bb[0], Bb[1], Bb[2], X)):
        if certainly_neg(g):
            return True
    return False


# ------------------------------------------------------------------ root enclosure in x
def _F_enclosure(mode, Bb, mids, X, direct_F):
    """One AD evaluation of F over B x X (4 variables) + one point evaluation.
    Returns (enc, num, gx, xm) with
        enc  = enclosure of F on B x X  (naive AD value intersected with the mean-value form),
        num  = enclosure of F(p, xm) for p in B (mean-value form around the centre),
        gx   = enclosure of dF/dx on B x X,  xm = X.mid(),
    or None if something is not finite."""
    xm = X.mid()
    vs = [AD.var(Bb[0], 0, 4), AD.var(Bb[1], 1, 4), AD.var(Bb[2], 2, 4), AD.var(X, 3, 4)]
    Fa = Fval(mode, vs[0], vs[1], vs[2], vs[3], direct_F)
    F0 = Fval(mode, mids[0], mids[1], mids[2], xm, direct_F)
    if not (Fa.v.is_finite() and F0.is_finite() and all(g.is_finite() for g in Fa.g)):
        return None
    num = F0
    for k in range(3):
        num = num + Fa.g[k] * (Bb[k] - mids[k])
    mvf = num + Fa.g[3] * (X - xm)
    enc = mvf.intersection(Fa.v) if mvf.overlaps(Fa.v) else mvf
    return enc, num, Fa.g[3], xm


def excludes_zero(v):
    return v.is_finite() and (v > 0 or v < 0)


def center_data(mode, B, X, direct_F='literal'):
    """Float data at the centre p0 of B: a root xc of F(p0, .) near X (Newton in floats,
    derivatives from AD at points) and the slopes phi_k = -F_k/F_x of the implicit function.
    Only used to choose coordinates (shear), Y and split directions -- never to decide."""
    p0 = [0.5 * (lo + hi) for lo, hi in B]
    xlo, xhi = float(X.lower()), float(X.upper())
    x = float(X.mid())
    for _ in range(6):
        vs = [AD.var(arb(p0[0]), 0, 4), AD.var(arb(p0[1]), 1, 4), AD.var(arb(p0[2]), 2, 4),
              AD.var(arb(x), 3, 4)]
        Fa = Fval(mode, vs[0], vs[1], vs[2], vs[3], direct_F)
        try:
            f, fx = float(Fa.v.mid()), float(Fa.g[3].mid())
            xn = x - f / fx
        except (ZeroDivisionError, OverflowError, ValueError):
            break
        if not math.isfinite(xn):
            break
        x = xn
    if not (xlo <= x <= xhi):
        x = float(X.mid())
    vs = [AD.var(arb(p0[0]), 0, 4), AD.var(arb(p0[1]), 1, 4), AD.var(arb(p0[2]), 2, 4),
          AD.var(arb(x), 3, 4)]
    Fa = Fval(mode, vs[0], vs[1], vs[2], vs[3], direct_F)
    try:
        fx = float(Fa.g[3].mid())
        phi = [-float(Fa.g[k].mid()) / fx for k in range(3)]
        if not all(map(math.isfinite, phi)):
            raise ValueError
    except (ZeroDivisionError, OverflowError, ValueError):
        phi = [0.0, 0.0, 0.0]
    return x, phi


def enclose_x(mode, B, direct_F='literal', max_pieces=256):
    """Return (pieces, status).  `pieces` is a list of arb balls whose union contains every x
    such that some (p1,p2,p3) in B has (p1,p2,p3,x) in the region with F = 0.
    Method: linear region bounds; then interval Newton on F in x, the other variables ranging
    over B, with the numerator F(B, xm) enclosed by a mean-value form; bisection in x where
    F_x is not bounded away from 0 or where the contraction stalls far above the expected width."""
    xr = linear_xrange(mode, B)
    if xr is None:
        return [], 'region', None
    Bb = [ball(lo, hi) for lo, hi in B]
    mids = [midpoint(lo, hi) for lo, hi in B]
    X0 = arb.union(xr[0], xr[1])
    xc, phi = center_data(mode, B, X0, direct_F)
    exp_rad = sum(abs(ph) * (hi - lo) / 2 for ph, (lo, hi) in zip(phi, B))
    tiny = 1e-13
    stack = [X0]
    out = []
    discarded_by_F = False
    while stack:
        X = stack.pop()
        if region_violated(mode, Bb, X):
            continue
        status = None
        while True:
            r = _F_enclosure(mode, Bb, mids, X, direct_F)
            if r is None:
                status = 'nan'
                break
            enc, num, gx, xm = r
            if excludes_zero(enc):
                X = None
                break
            if not excludes_zero(gx):
                status = 'gx'
                break
            N = xm - num / gx
            if not N.is_finite():
                status = 'nan'
                break
            if not N.overlaps(X):
                X = None
                break
            Xn = N.intersection(X)
            stalled = float(Xn.rad()) > 0.7 * float(X.rad())
            X = Xn
            if stalled:
                status = 'stalled'
                break
        if X is None:
            discarded_by_F = True
            continue
        rX = float(X.rad())
        if status == 'stalled' and rX <= 3 * exp_rad + tiny:
            out.append(X)
            continue
        if rX <= 0.25 * exp_rad + tiny or len(out) + len(stack) >= max_pieces:
            out.append(X)
            continue
        m = X.mid()
        stack.append(arb.union(X.lower(), m))
        stack.append(arb.union(m, X.upper()))
    if not out:
        return [], ('F' if discarded_by_F else 'region'), None
    # merge overlapping pieces (their hull is a valid, possibly larger, enclosure)
    out.sort(key=lambda z: float(z.lower()))
    merged = [out[0]]
    for z in out[1:]:
        if z.overlaps(merged[-1]) or float(z.lower()) - float(merged[-1].upper()) <= 1e-12:
            merged[-1] = arb.union(merged[-1], z)
        else:
            merged.append(z)
    return merged, 'ok', (xc, phi)


# ------------------------------------------------------------------ sheared coordinates
def xlin(p, mids, xc, phi):
    """x_lin(p) = xc + sum_k phi_k (p_k - m_k)  (xc, phi exact floats; p arb or AD)"""
    x = arb(xc)
    for k in range(3):
        x = (p[k] - mids[k]) * arb(phi[k]) + x
    return x


def shear(mode, B, X, xc, phi, direct_F='literal'):
    """Enclosure E of eps = x - x_lin(p) over all roots (p, x), p in B, x in X, of F.
    Interval Newton on G(p, eps) = F(p, x_lin(p) + eps) in eps with a mean-value numerator.
    Returns (E, contrib): E an arb ball, or None if it is proved that B x X contains no root;
    contrib[k] ~ share of dimension k in the radius of E (floats, used only for splitting)."""
    Bb = [ball(lo, hi) for lo, hi in B]
    mids = [midpoint(lo, hi) for lo, hi in B]
    E = X - xlin(Bb, mids, xc, phi)
    rads = [(hi - lo) / 2 for lo, hi in B]
    contrib = [abs(ph) * r for ph, r in zip(phi, rads)]
    for _ in range(10):
        em = E.mid()
        vs = [AD.var(Bb[0], 0, 4), AD.var(Bb[1], 1, 4), AD.var(Bb[2], 2, 4), AD.var(E, 3, 4)]
        xv = xlin(vs, mids, xc, phi) + vs[3]
        Ga = Fval(mode, vs[0], vs[1], vs[2], xv, direct_F)
        G0 = Fval(mode, mids[0], mids[1], mids[2], arb(xc) + em, direct_F)
        if not (Ga.v.is_finite() and G0.is_finite() and all(g.is_finite() for g in Ga.g)):
            break
        num = G0
        for k in range(3):
            num = num + Ga.g[k] * (Bb[k] - mids[k])
        if excludes_zero(num + Ga.g[3] * (E - em)) or excludes_zero(Ga.v):
            return None, None
        if not excludes_zero(Ga.g[3]):
            break
        N = em - num / Ga.g[3]
        if not N.is_finite():
            break
        if not N.overlaps(E):
            return None, None
        glo = float(abs(Ga.g[3]).lower())
        contrib = [float(abs(Ga.g[k]).upper()) * rads[k] / glo for k in range(3)]
        En = N.intersection(E)
        stalled = float(En.rad()) > 0.7 * float(E.rad())
        E = En
        if stalled:
            break
    return E, contrib


# ------------------------------------------------------------------ trS bound
def choose_Y(mode, B, xc):
    """float choice of Y: the minimiser of trS at (centre of B, xc)"""
    c = [0.5 * (lo + hi) for lo, hi in B] + [xc]
    cert = certificate(prims(mode, *c))
    Y = optimal_Y(cert)
    if not all(map(math.isfinite, Y)):
        return (0.0, 0.0)
    return Y


def trS_upper_sheared(mode, B, E, xc, phi, Y):
    """Rigorous upper bound of trS(Y) (i.e. of trS(y), y = Y/t) on the sheared box
    {(p, x_lin(p) + eps): p in B, eps in E}: min of the interval enclosure and of the
    mean-value form in the variables (p1, p2, p3, eps).  Returns (ub or None, |grad| or None)."""
    Bb = [ball(lo, hi) for lo, hi in B]
    mids = [midpoint(lo, hi) for lo, hi in B]
    em = E.mid()
    Ya = (arb(Y[0]), arb(Y[1]))
    x0 = arb(xc) + em
    f0 = trS(certificate(prims(mode, mids[0], mids[1], mids[2], x0)), Ya)
    vs = [AD.var(Bb[0], 0, 4), AD.var(Bb[1], 1, 4), AD.var(Bb[2], 2, 4), AD.var(E, 3, 4)]
    xv = xlin(vs, mids, xc, phi) + vs[3]
    fa = trS(certificate(prims(mode, vs[0], vs[1], vs[2], xv)), Ya)
    ub = fa.v.upper() if fa.v.is_finite() else None
    gmag = None
    if fa.v.is_finite() and f0.is_finite() and all(g.is_finite() for g in fa.g):
        enc = f0
        for k in range(3):
            enc = enc + fa.g[k] * (Bb[k] - mids[k])
        enc = enc + fa.g[3] * (E - em)
        u2 = enc.upper()
        ub = u2 if ub is None or u2 < ub else ub
        gmag = [float(abs(g).upper()) for g in fa.g]
    return ub, gmag
