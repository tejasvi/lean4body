"""Check (3), in computer algebra: D_e and Ph in the blow-up chart formulas (cert.chart) against the direct formulas
(cert.direct).  Under a = 1 + b al, c = b ga, x = 1/2 + b xi, with b > 0, the script checks that
    D_e(chart) = D_e(direct) for the six edges e,   Ph = b P.
The square roots are replaced by independent positive symbols, as in chart_identity.py.  Of the five values of psi
that occur, psi(R12,R13), psi(R12,R24), psi(R14,R12) and psi(R23,R12) are replaced by independent symbols as well
(keyed by the unordered pair of expanded arguments; psi is symmetric), and psi(R34,R12) by
(X^{-3/2} - Y^{-3/2}) / (Y - X) written with the square-root symbols.  An identity that holds with these four values
free holds in particular for their true values.  The seven comparisons run in parallel processes.
usage: python chart_identity_psi.py [TASK ...]   with TASK in D12, D13, D14, D23, D24, D34, P (default: all seven)"""
import sys, os
import multiprocessing
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'work', 'rig'))
import sympy as sp
import cert

bb, al, ga, xi = sp.symbols('b alpha gamma xi', positive=True)
rho, s = sp.symbols('rho s', positive=True)
a, c, x = 1 + bb * al, bb * ga, sp.Rational(1, 2) + bb * xi
rho2 = sp.expand(1 - 2 * ga * x + ga ** 2)
R34 = sp.expand(bb ** 2 * rho2)
ROOTS = {sp.expand(1 - x * x): s, sp.expand(rho2): rho, R34: bb * rho}
names = iter(('r_%d' % k for k in range(1, 20)))
PSI = {}
pnames = iter(('Psi_%d' % k for k in range(1, 20)))


def sq(v):
    key = sp.expand(v)
    if key not in ROOTS:
        ROOTS[key] = sp.Symbol(next(names), positive=True)
    return ROOTS[key]


def hd(X, Y):
    kx, ky = sp.expand(X), sp.expand(Y)
    if R34 in (kx, ky):
        sx, sy = sq(X), sq(Y)
        return (1 / (X * sx) - 1 / (Y * sy)) / (Y - X)
    key = frozenset((kx, ky))
    if key not in PSI:
        PSI[key] = sp.Symbol(next(pnames), positive=True)
    return PSI[key]


cert.sq = sq
cert.hd = hd
Qd = cert.direct(a, bb, c, x)
Qc = cert.chart(bb, al, ga, xi)
TASKS = ['D' + e for e in cert.EDGES] + ['P']


def zero(e):
    return sp.cancel(sp.together(e)) == 0


def run(task):
    if task == 'P':
        return zero(bb * (Qd['pos'] - Qd['neg']) - (Qc['pos'] - Qc['neg']))
    return zero(Qd['D'][task[1:]] - Qc['D'][task[1:]])


if __name__ == '__main__':
    tasks = sys.argv[1:] or TASKS
    with multiprocessing.get_context('fork').Pool(len(tasks)) as pool:
        for task, r in zip(tasks, pool.imap(run, tasks)):
            print(task, "identical (psi free, %d psi symbols):" % len(PSI), r, flush=True)
