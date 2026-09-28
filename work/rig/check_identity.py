"""independent checks: (1) F * prod s_e = -P at random points of C (not only on E);
(2) at random convex CCs, the Hessian of Phi (finite differences of the potential, masses from Dziobek) is positive
definite on nontrivial motions and its smallest generalised eigenvalue against K equals 1 - lambda_max(G)"""
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); sys.path.insert(0, '.')
import mpmath as mp, random
from cc4 import pos as cpos, Gmat, lam_max, theta, inD, F, dist, masses_from, dziobek, rigidity
import cert
mp.mp.dps = 40
random.seed(7)
# (1)
worst = 0
for it in range(2000):
    a = mp.mpf(random.uniform(0.58, 1.73)); b = mp.mpf(random.uniform(0.01, 1)); c = mp.mpf(random.uniform(0.01, 1.73)); x = mp.mpf(random.uniform(-0.2, 0.6))
    g, P = cert.direct_P(a, b, c, x)
    th = mp.acos(x)
    r = dist(a, b, c, th)
    prod_s = 1
    for v in r: prod_s *= v ** -3
    val = F(a, b, c, th) * prod_s + P
    worst = max(worst, abs(val) / (abs(P) + 1e-30))
print("(1) max |F prod s + P| / |P| over 2000 random points:", mp.nstr(worst, 5))
# (2)
def U_phi(Q, m, lam):
    tot = 0
    for i in range(4):
        for j in range(i + 1, 4):
            d = mp.sqrt((Q[2*i] - Q[2*j])**2 + (Q[2*i+1] - Q[2*j+1])**2)
            tot += m[i] * m[j] * (1 / d + lam * d * d / 2)
    return tot
n = 0
while n < 12:
    a = mp.mpf(random.uniform(0.58, 1.73)); b = mp.mpf(random.uniform(0, 1)); c = mp.mpf(random.uniform(0, 1.73))
    if not inD(a, b, c): continue
    ths = theta(a, b, c)
    if not ths: continue
    th = ths[0]; n += 1
    P = cpos(a, b, c, th)
    m = masses_from(P)
    A, r, s, lam, w = dziobek(P)
    Q0 = [P[i][k] for i in range(4) for k in range(2)]
    # gradient must vanish (CC check)
    grad = [mp.diff(lambda t, k=k: U_phi([Q0[i] + (t if i == k else 0) for i in range(8)], m, lam), 0) for k in range(8)]
    H = mp.matrix(8, 8)
    for k in range(8):
        for l in range(k, 8):
            f = lambda t1, t2: U_phi([Q0[i] + (t1 if i == k else 0) + (t2 if i == l else 0) for i in range(8)], m, lam)
            H[k, l] = H[l, k] = mp.diff(f, (0, 0), (1, 1))
    # basis of the complement of trivial motions: orthonormal complement of translations + rotation
    triv = [[1, 0] * 4, [0, 1] * 4, sum([[-P[i][1], P[i][0]] for i in range(4)], [])]
    Mb = mp.matrix([[mp.mpf(v) for v in row] for row in triv]).T      # 8x3
    Qm, Rm = mp.qr(Mb, mode='full')
    Nb = Qm[:, 3:8]                                                     # 8x5 complement
    Hr = Nb.T * H * Nb
    ev = sorted(mp.eigsy(Hr)[0])
    # K and generalised eigenvalue
    G, dd = Gmat(P)
    lmax = lam_max(G)
    Rr, am = rigidity(P)
    sig = None
    # sigma from m1 m2 w12 = sigma A1 A2
    sig = m[0] * m[1] * w[(0, 1)] / (A[0] * A[1]) if isinstance(w, dict) else None
    Kd = mp.matrix(6, 6)
    for e, (i, j) in enumerate([(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]):
        Kd[e, e] = 3 * m[i] * m[j] * r[(i, j)] ** -3 if isinstance(r, dict) else 3 * m[i] * m[j] * s[e]
    K = Rr.T * Kd * Rr
    Kr = Nb.T * K * Nb
    # generalised eigenvalues of (Hr, Kr): Kr^{-1/2} Hr Kr^{-1/2}
    E, V = mp.eigsy(Kr)
    Kmh = V * mp.diag([1 / mp.sqrt(e) for e in E]) * V.T
    mu = sorted(mp.eigsy(Kmh * Hr * Kmh)[0])
    print("(2) |grad| %.1e  min eig Hess Phi (nontrivial) %s  min gen. eig mu %s  1 - lam_max(G) %s" % (
        float(max(abs(g) for g in grad)), mp.nstr(ev[0], 8), mp.nstr(mu[0], 12), mp.nstr(1 - lmax, 12)))
