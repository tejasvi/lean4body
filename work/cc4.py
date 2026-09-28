"""Core geometry for convex 4-body central configurations in Corbera-Cors-Roberts (CCR) coordinates, in mpmath
(arbitrary precision; boundary asymptotics need it).  No side effects on import.
CCR: q1=(1,0), q2=a(cos t, sin t), q3=(-b,0), q4=-c(cos t, sin t),  pi/3 < t <= pi/2.
Dziobek: m_i m_j w_ij = sigma A_i A_j,  w_ij = s_ij - lam,  s_ij = r_ij^-3,  A_i = signed area of the triangle
without body i (A1, A3 > 0 > A2, A4).  Hessian of Phi = sum m_i m_j (1/r_ij + lam r_ij^2/2) at the CC:
    Hess/|sigma| = K - a^T a,   K = sum_e W_e (e_e.(v_i - v_j))^2,  W_e = 3 s_e |A_i A_j / w_e|,  a v = sum A_i v_i.
Nondegeneracy (and local minimality mod rotations) <=> lambda_max(G) < 1,  G = a K^+ a^T = min_{bR=a} b D b^T,
D_e = 1/W_e (bar compliance).  Explicit particular solution of bR = a: beta_e = -sum_{i not in e} (dA_i/dr_e) q_i."""
import mpmath as mp

PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def pos(a, b, c, t):
    ct, st = mp.cos(t), mp.sin(t)
    return [(mp.mpf(1), mp.mpf(0)), (a*ct, a*st), (-b, mp.mpf(0)), (-c*ct, -c*st)]


def dist(a, b, c, t):
    ct = mp.cos(t)
    r12 = mp.sqrt(a**2 - 2*a*ct + 1); r23 = mp.sqrt(a**2 + 2*a*b*ct + b**2); r13 = 1 + b
    r14 = mp.sqrt(c**2 + 2*c*ct + 1); r34 = mp.sqrt(b**2 - 2*b*c*ct + c**2); r24 = a + c
    return r12, r13, r14, r23, r24, r34


def F(a, b, c, t):
    r12, r13, r14, r23, r24, r34 = dist(a, b, c, t)
    return (r24**3-r14**3)*(r13**3-r12**3)*(r23**3-r34**3) - (r12**3-r14**3)*(r24**3-r34**3)*(r13**3-r23**3)


def inC(a, b, c, t):
    r12, r13, r14, r23, r24, r34 = dist(a, b, c, t)
    return r13 > r12 and r24 > r12 and r12 >= r14 and r12 >= r23 and r14 >= r34 and r23 >= r34


def inD(a, b, c):
    if not (0 < c <= a and 0 < b <= 1): return False
    if c >= (b**2+2*b)/a: return False
    if c <= -a + mp.sqrt(a**2+b): return False
    if a <= 1:
        return c > (-a + mp.sqrt(max(4-3*a**2, 0)))/2
    return b > (-1 + mp.sqrt(4*a**2-3))/2


def theta(a, b, c, n=64):
    """all roots of F in (pi/3, pi/2] that lie in C (normally exactly one)"""
    lo, hi = mp.pi/3, mp.pi/2
    if abs(F(a, b, c, hi)) < mp.mpf(10)**(-mp.mp.dps+8):
        return [hi]
    ts = [lo + (hi-lo)*k/n for k in range(1, n+1)]
    vs = [F(a, b, c, x) for x in ts]
    out = []
    for k in range(n-1):
        if vs[k]*vs[k+1] < 0:
            x = mp.findroot(lambda x: F(a, b, c, x), (ts[k], ts[k+1]), solver='anderson')
            if inC(a, b, c, x): out.append(x)
    if not out:          # root very close to pi/3 (near P1): bisect on (pi/3, ts[0])
        x0 = lo + (ts[0]-lo)*mp.mpf(10)**-30
        if F(a, b, c, x0)*vs[0] < 0:
            x = mp.findroot(lambda x: F(a, b, c, x), (x0, ts[0]), solver='anderson')
            if inC(a, b, c, x): out.append(x)
    return out


def tri(P, i, j, k):
    return ((P[j][0]-P[i][0])*(P[k][1]-P[i][1]) - (P[k][0]-P[i][0])*(P[j][1]-P[i][1]))/2


def areas(P):
    return [tri(P, 1, 2, 3), -tri(P, 0, 2, 3), tri(P, 0, 1, 3), -tri(P, 0, 1, 2)]


def dziobek(P):
    A = areas(P)
    r = {e: mp.sqrt((P[e[0]][0]-P[e[1]][0])**2 + (P[e[0]][1]-P[e[1]][1])**2) for e in PAIRS}
    s = {e: r[e]**-3 for e in PAIRS}
    lam = (s[(0, 1)]*s[(2, 3)] - s[(0, 2)]*s[(1, 3)])/(s[(0, 1)] + s[(2, 3)] - s[(0, 2)] - s[(1, 3)])
    w = {e: s[e]-lam for e in PAIRS}
    return A, r, s, lam, w


def masses_from(P):
    """masses (normalised m1=1) from the Dziobek relations m_i m_j w_ij = sigma A_i A_j"""
    A, r, s, lam, w = dziobek(P)
    # m_j/m_1 = (A_j/A_1) * (w_1k w_1l ... ) : use m1 m_j w_1j = sigma A1 A_j and m_j m_k w_jk = sigma A_j A_k
    # => m_j^2 = sigma A_j^2 * (A_1 A_k ... ) ; simplest: sigma from (1,2),(1,3),(2,3):
    # m1^2 = sigma A1^2 w_23/(w_12 w_13)
    def w_(i, j): return w[(min(i, j), max(i, j))]
    m = []
    for i in range(4):
        j, k = [x for x in range(4) if x != i][:2]
        m.append(mp.sqrt(abs(A[i]**2 * w_(j, k)/(w_(i, j)*w_(i, k)))))
    return [x/m[0] for x in m]


def angle(P, l, j, k):
    """angle at vertex l of the triangle (l, j, k)"""
    ux, uy = P[j][0]-P[l][0], P[j][1]-P[l][1]
    vx, vy = P[k][0]-P[l][0], P[k][1]-P[l][1]
    return mp.atan2(abs(ux*vy-uy*vx), ux*vx+uy*vy)


def Gmat(P):
    """returns (G 2x2 as mp.matrix, dict with pieces).  G = min_{bR=a} b D b^T via the covariance formula."""
    A, r, s, lam, w = dziobek(P)
    D, tau, beta = {}, {}, {}
    for e in PAIRS:
        i, j = e
        D[e] = abs(w[e])/(3*s[e]*abs(A[i]*A[j]))
        tau[e] = A[i]*A[j]*r[e]
        bx = by = mp.mpf(0)
        for l in range(4):
            if l in e: continue
            # triangle without body l contains bar e; dA_l/dr_e = sign(A_l) (r_e/2) cot(angle opposite e)
            k = [x for x in range(4) if x not in e and x != l][0]     # third vertex of that triangle
            g = angle(P, k, i, j)
            dA = mp.sign(A[l])*r[e]/2*mp.cot(g)
            bx -= dA*P[l][0]; by -= dA*P[l][1]
        beta[e] = (bx, by)
    T = sum(D[e]*tau[e]**2 for e in PAIRS)
    mx = sum(D[e]*tau[e]*beta[e][0] for e in PAIRS)/T
    my = sum(D[e]*tau[e]*beta[e][1] for e in PAIRS)/T
    G = mp.matrix(2, 2)
    for e in PAIRS:
        u = (beta[e][0]-mx*tau[e], beta[e][1]-my*tau[e])
        G[0, 0] += D[e]*u[0]*u[0]; G[0, 1] += D[e]*u[0]*u[1]; G[1, 1] += D[e]*u[1]*u[1]
    G[1, 0] = G[0, 1]
    return G, dict(A=A, r=r, s=s, lam=lam, w=w, D=D, tau=tau, beta=beta)


def lam_max(G):
    tr = G[0, 0]+G[1, 1]; det = G[0, 0]*G[1, 1]-G[0, 1]**2
    return tr/2 + mp.sqrt(max(tr**2/4-det, 0))


def rigidity(P):
    """6x8 rigidity matrix and the 2x8 load matrix a (for cross-checks)"""
    A, r, s, lam, w = dziobek(P)
    R = mp.matrix(6, 8)
    for k, (i, j) in enumerate(PAIRS):
        ex, ey = (P[i][0]-P[j][0])/r[(i, j)], (P[i][1]-P[j][1])/r[(i, j)]
        R[k, 2*i], R[k, 2*i+1], R[k, 2*j], R[k, 2*j+1] = ex, ey, -ex, -ey
    am = mp.matrix(2, 8)
    for i in range(4):
        am[0, 2*i] = A[i]; am[1, 2*i+1] = A[i]
    return R, am
