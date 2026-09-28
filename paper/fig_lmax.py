"""Figure for the paper: the sharp constant lambda_max(G) and tr G = min_y tr S(y) at sampled normalised convex central
configurations, against the coordinate b (floating point; illustration only, not part of the proof).
Run from any directory; writes fig_lmax.pdf and fig_lmax.png next to this file."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'work', 'rig'))
import cert
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rng = np.random.default_rng(11)


def survey(MODE, NS, NX=400):
    if MODE == 'direct':
        p1 = rng.uniform(0.5, 1.75, NS); p2 = rng.uniform(0.125, 1.0, NS); p3 = rng.uniform(0, 1.75, NS)
        lo, hi, FUN = -1.0, 1.0, cert.direct
    else:
        p1 = 0.125 * rng.uniform(0, 1, NS) ** 3; p2 = rng.uniform(-2, 1, NS); p3 = rng.uniform(0, 3, NS)
        lo, hi, FUN = -3.0, 3.0, cert.chart
    xs = np.linspace(lo, hi, NX)
    out = []
    with np.errstate(all='ignore'):
        Pv = []
        for x in xs:
            Q = FUN(p1, p2, p3, np.full(NS, x))
            ok = np.all(np.array(Q['g']) >= 0, axis=0)
            Pv.append(np.where(ok, Q['pos'] - Q['neg'], np.nan))
        Pv = np.array(Pv)
        for k in range(NX - 1):
            idx = np.nonzero(Pv[k] * Pv[k + 1] < 0)[0]
            if len(idx) == 0:
                continue
            xl = np.full(len(idx), xs[k]); xh = np.full(len(idx), xs[k + 1])
            A, B, C = p1[idx], p2[idx], p3[idx]; fl = Pv[k][idx]
            for it in range(60):
                xm = 0.5 * (xl + xh); Q = FUN(A, B, C, xm); fm = Q['pos'] - Q['neg']; left = fm * fl <= 0
                xh = np.where(left, xm, xh); xl = np.where(left, xl, xm); fl = np.where(left, fl, fm)
            Q = FUN(A, B, C, 0.5 * (xl + xh))
            y = cert.opt_y(Q)
            S11, S12, S22 = cert.Smat(Q, y)
            tr = S11 + S22; det = S11 * S22 - S12 ** 2
            lmax = tr / 2 + np.sqrt(np.maximum(tr * tr / 4 - det, 0))
            bb = A if MODE == 'chart' else B
            out.append(np.stack([bb, lmax, tr], axis=1))
    return np.concatenate(out)


R = np.concatenate([survey('direct', 60000), survey('chart', 60000)])
R = R[np.isfinite(R).all(axis=1) & (R[:, 0] > 0)]
print("%d configurations; b in [%.2e, %.3f]; max lambda_max %.6f; max tr G %.6f" % (
    len(R), R[:, 0].min(), R[:, 0].max(), R[:, 1].max(), R[:, 2].max()))

INK, INK2 = '#0b0b0b', '#52514e'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2,
                     'ytick.color': INK2, 'font.family': 'serif'})
fig, ax = plt.subplots(figsize=(6.0, 3.4))
ax.scatter(R[:, 0], R[:, 2], s=1.2, color='#eb6834', alpha=0.35, linewidths=0, rasterized=True,
           label=r'$\operatorname{tr} G=\min_y \operatorname{tr} S(y)$')
ax.scatter(R[:, 0], R[:, 1], s=1.2, color='#2a78d6', alpha=0.35, linewidths=0, rasterized=True,
           label=r'$\lambda_{\max}(G)$ (sharp constant)')
for yv, txt, ls in ((1.0, 'nondegeneracy\nthreshold 1', '--'), (0.75, 'certified\nbound 3/4', '--'), (1 / 3, '1/3', ':')):
    ax.axhline(yv, color=INK2, lw=0.8, ls=ls)
    ax.text(1.015, yv, txt, color=INK2, fontsize=8, va='center', transform=ax.get_yaxis_transform())
ax.set_xscale('log')
ax.set_xlim(1e-4, 1.0)
ax.set_ylim(0, 1.08)
ax.set_xlabel(r'$b$  (distance from $q_3$ to the intersection of the diagonals, with $|q_1|=1$)')
ax.set_ylabel('value at the configuration')
ax.axvline(0.125, color=INK2, lw=0.5)
ax.text(0.128, 0.05, 'direct\nregion', color=INK2, fontsize=7)
ax.text(0.03, 0.05, 'blow-up\nchart', color=INK2, fontsize=7)
leg = ax.legend(loc='center left', bbox_to_anchor=(0.0, 0.84), frameon=False, markerscale=6, fontsize=8)
for lh in leg.legend_handles:
    lh.set_alpha(1)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.subplots_adjust(right=0.84)
fig.savefig(os.path.join(HERE, 'fig_lmax.pdf'), dpi=300)
fig.savefig(os.path.join(HERE, 'fig_lmax.png'), dpi=160)
