# Independent interval-arithmetic verification of the trS inequality (SPEC.md)

I wrote this code from `SPEC.md` alone. Development and all the runs used only files in this directory.
Python: `python` (python-flint 0.9 arb, mpmath, numpy, sympy).

**Independence note.** There was one lapse, after all the runs had finished.
- While checking a number for section 8, I text-searched my own session transcript, which is stored
  outside this directory.
- The matches did not contain the number. Among them were a few lines of numerical output about this
  problem's chart (values of lam and F/b^2). I cannot tell where those lines came from.
- I stopped there and used nothing from them. No computation or result changed afterwards.

## 1. Result

**Both claims of the SPEC are verified, with no unproved box.** I also ran the verifier a second time
with the smaller threshold 19/32 = 0.59375. Both claims hold at that threshold too.

**Direct claim.** Take any (a, b, c, x) such that
- (a, b, c) is in OMEGA_direct = [1/2, 7/4] x [1/8, 1] x [0, 7/4] and x is in [-1, 1],
- the six region inequalities hold, and
- F(a, b, c, x) = 0.

Then (a, b, c) lies in one of the accepted leaf boxes B of the run. The leaves are dyadic boxes. The
region tests showed that the rest of OMEGA_direct contains no point of the region.

For the vector y_B stored with that leaf (a pair of doubles):

    trS(y_B)(a, b, c, x) <= ub_B < 3/4        (run 1)
    trS(y_B)(a, b, c, x) <= ub_B < 19/32      (run 2)

**Chart claim.** Take any (b, al, ga, xi) with b in (0, 1/8] and (al, ga) in [-2, 1] x [0, 3] such that
(a, b, c, x) = (1 + b al, b, b ga, 1/2 + b xi) is in Z. Then (b, al, ga) lies in an accepted leaf B, and,
for the Y_B stored with that leaf,

    trS(Y_B / b) <= ub_B < 3/4 (run 1),  < 19/32 (run 2).

- The check is done in the rescaled formulas, on boxes that include b = 0.
- The SPEC says that such points have al in (-2, 1) and ga in [0, 17/6), and allows this to be taken as
  given. I did so. With it, the chart claim covers every point of Z with 0 < b <= 1/8.

Each run is a complete proof on its own. The two runs have different leaves and y's.

| run | threshold | initial cells | leaves accepted | discarded (region) | discarded (F != 0) | unproved | bisections / max depth | max certified ub of trS | wall (16 procs) | CPU |
|---|---|---|---|---|---|---|---|---|---|---|
| direct | 3/4 | 7840 | 190 384 | 15 630 | 0 | **0** | 198 174 / 13 | 0.7499990406209058 | 41.8 s | 664 s |
| chart | 3/4 | 4608 | 20 189 | 4 381 | 0 | **0** | 19 962 / 6 | 0.7499902101804914 | 8.7 s | 134 s |
| direct | 19/32 | 7840 | 804 708 | 42 166 | 0 | **0** | 839 034 / 19 | 0.5937499913476635 | 172.9 s | 2667 s |
| chart | 19/32 | 4608 | 95 972 | 6 325 | 0 | **0** | 97 689 / 10 | 0.5937499972464421 | 41.3 s | 629 s |

Notes on the table:
- **The ub column.** The reported ub is a double that is >= the certified arb bound. It is rounded up
  with `nextafter` whenever the double nearest to the bound is below it.
- **Why the maximum is so close to the threshold.** A leaf is accepted as soon as its bound falls below
  the threshold, so the largest certified bound always sits just under it. The run at 19/32 is the
  informative one.
- **Margins.** A float search puts the maximum over Z of min_y trS at about **0.5927** (section 7). So
  19/32 is about 0.001 above the true value, and 3/4 leaves a margin of about 0.157.
- **One piece per leaf.** Every accepted leaf holds exactly one x-piece (the logs report
  "x-pieces" = leaves).
- **What the discard columns count.**
  - "region": every x-piece of the box was removed by a region test alone (empty linear x-range, a
    violated pair condition, or a region quantity certainly negative on B x X).
  - "F != 0": no piece survived and at least one was removed by an F test, either in the root enclosure
    or in the shear step.
- **No F-discards.** This column is 0 in every run, so every discarded box was discarded by region tests
  alone (section 7 explains why).
- **Every leaf is accounted for.** In every run, accepted + discarded = initial cells + bisections.

## 2. Files (all in `work/indep/`)

**Code**

| file | role |
|---|---|
| `ad.py` | forward-mode automatic differentiation over arb balls; `SQRT`, which dispatches on the number type |
| `formulas.py` | literal SPEC formulas (reference only), and the scaled formulas used by the verifier (direct: t = 1; chart: t = b): `certificate`, `trS`, `Ftil`, `F_direct`, region |
| `core.py` | rigorous per-box routines: linear x-range, root enclosure, shear, trS bound |
| `verify.py` | branch-and-bound driver (multiprocessing, at most 16 workers) |
| `checks.py` | the sanity checks required before the branch-and-bound (sympy, and mpmath at 50 and 200 digits) |
| `explore.py`, `maxsearch.py` | float exploration of Z, and a float search for max min_y trS (not part of the proof) |
| `check_leaves.py` | end-to-end sampling test of the accepted leaves (not part of the proof) |
| `negative_control.py` | checks that thresholds below the true maximum fail (not part of the proof) |
| `run_all.sh`, `run_tight.sh` | the exact commands of the four runs and their sampling tests |

**Logs**

| file | content |
|---|---|
| `checks.log` | `checks.py 300`, the checks done before the branch-and-bound |
| `direct_run.log`, `.json` / `chart_run.log`, `.json` | the runs at threshold 3/4 |
| `direct_run_19_32.log`, `.json` / `chart_run_19_32.log`, `.json` | the runs at threshold 19/32 |
| `check_leaves_{direct,chart}[_19_32].log` | results of the sampling test |
| `negative_control.log`, `explore_{direct,chart}.log`, `maxsearch_{direct,chart}.log` | tests and exploration |

The [README](../../README.md) at the root of the repository lists the commands that write these logs.

What the run logs contain:
- The `.log` files report progress every 250 cells; then the totals, the maximum bound and the leaf where
  it occurs, and a histogram of the leaf bounds.
- The `.json` files hold the same totals and the (empty) list of failures.

**Accepted leaves.** Each row holds: cell index, box, x-piece, x_lin data, E, Y, ub.
- `run_all.sh` writes `direct_leaves.npy` and `chart_leaves.npy` (threshold 3/4), and `run_tight.sh`
  writes `direct_leaves_19_32.npy` (129 MB) and `chart_leaves_19_32.npy`.

## 3. Formulas

`formulas.literal` implements section 1 of the SPEC literally: positions, R, r, s_e, signed areas from the
triangle formula, lam, w, D, tau, beta, F, region. It serves only as the reference for the checks.

The verifier uses exact algebraic rewrites of the same formulas, written for a scale t.
- Direct coordinates: t = 1.
- Chart: a = 1 + b al, c = b ga, x = 1/2 + b xi, and t = b.

**Closed forms (sympy-verified).** With s = sqrt(1 - x^2):

    R12 = 1 + a^2 - 2ax    R13 = (1+b)^2       R14 = 1 + c^2 + 2cx
    R23 = a^2 + b^2 + 2abx R24 = (a+c)^2       R34 = b^2 + c^2 - 2bcx
    A1 = b(a+c)s/2   A2 = -(1+b)c s/2   A3 = (a+c)s/2   A4 = -a(1+b)s/2

**Scaled quantities, defined exactly.** R_e - R12 = t E_e (e = 13, 14, 23, 24), R34 = t^2 rho^2,
A1 = t A1h, A2 = t A2h. In the chart:

    E13 = 2 + b + a(2xi - al)       E14 = -(a+c)(al - ga - 2xi)
    E23 = (1+b)(1 + al + 2a xi)     E24 = 2a ga + b ga^2 + al + 2a xi
    rho^2 = 1 - 2 ga x + ga^2       A1h = A3 = (a+c)s/2     A2h = -(1+b) ga s/2
    R12 = 1 + b a (al - 2xi)        s^2 = (1/2 - b xi)(3/2 + b xi)

Apart from s, each of these is a polynomial in (b, al, ga, xi). sympy checked that the division by b
leaves zero remainder.

**The weights.** Write S_e = (s12 - s_e)/t. Using
r_e^3 - r_f^3 = (R_e - R_f)(R_e + r_e r_f + R_f)/(r_e + r_f):

    S_e = E_e (R_e + r_e r12 + R12) / ((r_e + r12) r_e^3 r12^3).

Let s34h = t^3 s34 = rho^-3, and Dn = s34h + t^3 (s12 - s13 - s24). Then Dn = t^3 (s12 + s34 - s13 - s24),
which is t^3 times the denominator of lam. The w_e factorise as follows (sympy-verified identities in the
s_e):

    w12 = (s12-s13)(s12-s24)/(s12+s34-s13-s24) = t^5 S13 S24/Dn       =: t^5 w12h
    w13 = t S13 (t^3 s13 - s34h)/Dn
    w24 = t S24 (t^3 s24 - s34h)/Dn
    w34 = t^-3 (s34h - t^3 s13)(s34h - t^3 s24)/Dn
    w14 = w12 - (s12 - s14) = t (t^4 w12h - S14)
    w23 = w12 - (s12 - s23) = t (t^4 w12h - S23)

**D_e, tau_e and beta_e in scaled form.**
- **D_e = -w_e/(3 s_e A_i A_j).** The powers of t cancel in every D_e except D12, which keeps a factor
  t^3: D12 = -t^3 w12h/(3 s12 A1h A2h). For example, D13 = -w13h/(3 s13 A1h A3) and
  D34 = -w34h/(3 s34h A3 A4).
- **tau_e = t T_e**, with T12 = t A1h A2h r12, T34 = A3 A4 rho, and the other T_e equal to A_i A_j r_e
  with A1, A2 replaced by A1h, A2h.
- **beta_e.** Every term with 1/A1 = 1/(t A1h) or 1/A2 = 1/(t A2h) carries a factor that vanishes at the
  corner, and that factor is divided by t.
  - For e = 13, 14, 23, 24 the factor is N_{e,l}. For example, N_{13,2} = R14 + R34 - R13
    = t (E14 - E13 + t rho^2).
  - For e = 34 the factor is r34 = t rho.

So with y = Y/t, trS(y) = sum_e D_e |beta_e + Y T_e|^2. `certificate()` is the same code for both modes.

These rewrites divide only by quantities that the SPEC's own formulas also divide by (r_e, A_i, s_e, the
lam denominator), or by the positive quantity r_e + r_f. So they agree with the SPEC wherever the SPEC's
formulas are defined.

### The rescaled equation in the chart (justification)

Factor each cube difference in F:
- r24^3 - r14^3 = t (E24 - E14) G(r24, r14)
- r13^3 - r12^3 = t E13 G(r13, r12)
- r12^3 - r14^3 = -t E14 G(r12, r14)
- r13^3 - r23^3 = t (E13 - E23) G(r13, r23)
- r34^3 = t^3 rho^3

where G(u,v) = (u^2 + uv + v^2)/(u + v). This gives, for every t,

    F = t^2 * Ftil,
    Ftil = (E24-E14) G(r24,r14) E13 G(r13,r12) (r23^3 - t^3 rho^3)
           + E14 G(r12,r14) (r24^3 - t^3 rho^3) (E13-E23) G(r13,r23).

Hence **for b > 0, F = 0 if and only if Ftil = 0** (Ftil = F/b^2). The chart branch-and-bound encloses the
roots of Ftil. Direct mode uses the SPEC's F exactly as written (r^3 = R sqrt(R), with R_ij the exact
polynomials).

**Ftil is nondegenerate at b = 0.** There every r_e is 1 except r34, all G are 3/2, and

    Ftil(0, al, ga, xi) = (9/4)(E24 E13 - E14 E23) = (9/4) [3 al + 3 ga - 3 al ga + 2 xi (1 + ga)].

So dFtil/dxi = (9/2)(1 + ga) >= 9/2 > 0, and the root xi = -3(al + ga - al ga)/(2(1 + ga)) is simple. I
divided out no factor other than b^2.

**All chart quantities are finite at b = 0 on the relevant set.**
- A1h, A3 and A4 are nonzero.
- A2h is nonzero, because the region forces ga >= 1/2 at b = 0 (and ga > 0.43 for b <= 1/8).
- Dn tends to rho^-3 > 0.
- rho^2 = (ga - 1/2)^2 + 3/4 at b = 0.

## 4. Region handling

Each region inequality is linear in x (in xi in the chart). For example:
- R13 >= R12 is 2ax >= a^2 - 2b - b^2;
- R12 >= R14 is (a+c)(a - c - 2x) >= 0, that is, 2x <= a - c.

**The x-range (`core.linear_xrange`).**
- Direct:
  - lower bounds: (a^2-2b-b^2)/(2a), (1-2ac-c^2)/(2a), (c-a)/(2b), -(1-b)/(2c) (used only if c > 0 on
    the box), and -1;
  - upper bounds: (a-c)/2, (1-b)/(2a), 1.
- Chart, from the four inequalities that are divisible by b, with u = a = 1 + b al >= 3/4:
  - xi >= al/2 - (2+b)/(2u) and xi >= -ga - (b ga^2 + al)/(2u);
  - xi <= (al - ga)/2 and xi <= -(1+al)/(2u).
- The bounds are evaluated in arb over the box. The lower end is the max of the lower endpoints of the
  lower bounds, and the upper end is the min of the upper endpoints of the upper bounds.
- Using only some of the inequalities is valid, since dropping one only enlarges the range.

**Pair conditions.** Requiring each lower bound to be <= each upper bound gives necessary conditions that
do not involve x. I derived each by hand; they are listed in `core.py`. A box is discarded if one of them
is certainly negative on the whole box.
- Direct: 2b + b^2 - ac, 1 + b + b^2 - a^2, a^2 + ac + c^2 - 1, 2ac + c^2 - b, c(a-c) + 1 - b, a - c,
  b(1-b) - a(c-a).
- Chart: 2 + b - ga u, 1 + b - al(2 + b al), (al+ga)u + b ga^2 + al, 2 ga u + b ga^2 - 1.

**Discarding x-pieces.** An x-piece X is also dropped if any of the six region quantities, evaluated in
arb on B x X, is certainly negative.
- Direct: R13-R12, R24-R12, R12-R14, R12-R23, R14-R34, R23-R34.
- Chart: the first four divided by b, and the last two.

## 5. Branch-and-bound (`verify.py`, `core.py`)

**Grids.**
- Direct: 20 x 14 x 28 cells of side 1/16, covering OMEGA_direct exactly.
- Chart: 8 x 24 x 24 cells of sides 1/64 x 1/8 x 1/8, covering [0, 1/8] x [-2, 1] x [0, 3] (b = 0 is
  included).
- All box endpoints are dyadic, hence exact in double and in arb. Every midpoint is checked to be exact.

**Precision.** arb at 96 bits. `AD` evaluates a formula on arb balls together with its gradient. This gives
rigorous enclosures of the value and of the gradient over a box.

For each box B = B1 x B2 x B3 with centre m, the steps are as follows.

**(i) Enclosing x.**
1. Take X0 from the linear x-range above.
2. Contract X by interval Newton on F in x, with p = (p1, p2, p3) ranging over B. Let x_m be the midpoint
   of X. The step uses a mean-value numerator:

       N(X) = x_m - [F(m, x_m) + sum_k dF/dp_k(B x X) (B_k - m_k)] / dF/dx(B x X),   X <- X ∩ N(X).

   **Why this is valid.** Let (p, x*) be a root with p in B and x* in X.
   - By the mean-value theorem in p at fixed x_m, F(p, x_m) = F(m, x_m) + sum_k dF/dp_k(zeta, x_m)(p_k - m_k),
     with (zeta, x_m) in B x X. So F(p, x_m) lies in the bracket.
   - By the mean-value theorem in x, 0 = F(p, x*) = F(p, x_m) + dF/dx(p, eta)(x* - x_m), with eta in X.
   - Hence x* is in N(X).
3. X is discarded if N(X) ∩ X is empty, or if the enclosure of F on B x X excludes 0. That enclosure is
   the naive AD value intersected with the mean-value form.
4. **When the contraction stops.** It stops when dF/dx(B x X) contains 0, when an evaluation is not
   finite, or when a step shrinks the radius by less than 30%.
5. **Keep or bisect.** Let r_exp = sum_k |phi_k| rad(B_k) be the expected spread of the root over B
   (phi is defined below).
   - X is kept as a piece if its radius is at most r_exp/4 + 1e-13, or at most 3 r_exp + 1e-13 after a
     stalled contraction.
   - X is also kept once 256 pieces (kept plus pending) exist.
   - Otherwise X is bisected and both halves are processed again.
   - Finally, pieces that overlap or are less than 1e-12 apart are merged into their hull.

**Shear.** Take x_c, a float root at m, and phi_k = -F_k/F_x there. These floats only choose
coordinates. Let x_lin(p) = x_c + sum_k phi_k (p_k - m_k). Every root (p, x) with p in B and x in X
satisfies eps = x - x_lin(p) in E.
- E starts from E0 = X - x_lin(B).
- E is then contracted by the same Newton step applied to G(p, eps) = F(p, x_lin(p) + eps).
- The piece is dropped if this step proves that it contains no root.

**(ii) Discarding.** B is discarded if every piece is eliminated.

**(iii) Bounding trS.** Take Y as the float minimiser of the quadratic Y -> trS at the point
(m, x_c + mid E), that is, Y = -sum D_e T_e beta_e / sum D_e T_e^2. Let g(p, eps) = trS(Y) at
(p, x_lin(p) + eps). Then

    ub = min( upper(g(B x E)),  upper( g(m, mid E) + sum_k dg/dp_k(B x E)(B_k - m_k) + dg/deps(B x E)(E - mid E) ) ).

- The first term g(B x E) is the AD value. The second is the mean-value form in the sheared variables
  (p, eps), which range over the box B x E.
- The leaf is accepted if ub < threshold. This is an arb comparison with an exact dyadic threshold.
- In the sheared variables, dg/dp_k is approximately the derivative of trS along Z, which is small. So
  the bound is far tighter than it would be in (p, x).

**(iv) Splitting.** Otherwise B is bisected in the dimension k that maximises
|dg/dp_k| rad_k + |dg/deps| contrib_k. Here contrib_k = |G_k| rad_k / |G_eps| is the share of dimension k in
rad E.

**Give-up rules.** A box is recorded as unproved in either of two cases:
- every side is below 2^-18 of its initial cell;
- the cell's budget of B&B nodes is used up (200 000 at 3/4, 300 000 at 19/32).

Neither happened.

### Rigor notes
- **Decisions.** Every decision to discard or accept is an arb comparison. Such a comparison is true only
  if it is certain.
- **Floats.** Floats (x_c, phi, Y, split choices) are converted to arb exactly. They only define
  coordinates and y, so their accuracy does not affect validity.
- **NaN handling in python-flint 0.9 (tested).**
  - 1/[-1,1], 1/[0,1], 0/[-1,1], 0/0 and sqrt([-1e-30, 1]) all return NaN, and 0*NaN is NaN.
  - The naive interval value is used as a bound only if it is finite.
  - The mean-value forms and the Newton steps are used only if the value and every gradient component
    are finite.
- **Square roots.**
  - The arguments are squared distances, 1 - x^2, and rho^2. These are >= 0 at every real point, so the
    ball is clipped to [0, inf) before sqrt. Clipping loses no real point.
  - In AD the derivative 1/(2r) of every square root is formed. In the trS formulas every square root is
    also a divisor.
  - So a finite value certifies that every denominator is != 0 on the whole box; the naive value is then
    a valid enclosure.
  - A finite gradient further certifies that every sqrt argument is > 0 on the whole box. The functions
    are then C^1 on the (convex) box, and the mean-value forms are valid.
  - The value path of F_direct does not divide by square roots. Clipping is harmless there, because the
    true value is still enclosed.
- **Pieces.** Pieces are dropped only by certain tests: a region inequality certainly violated on B x X,
  an enclosure of F excluding 0, or an empty Newton intersection. Merging pieces into their hull only
  enlarges the enclosure.

## 6. Checks done before the branch-and-bound (`checks.py` -> `checks.log`)

All checks passed. **No SPEC formula failed.** "Error" below means |u - v|/max(1, |v|). The rigidity
residuals are absolute.

**A. sympy.**
- R_ij closed forms and signed-area closed forms.
- sum A_i = 0 and sum A_i q_i = 0, identically.
- E_e in direct and chart form, including polynomial division by b with zero remainder.
- R34 = b^2 rho^2, A1 = b A1h, A2 = b A2h.
- The w_e factorisations, and the cube-difference identity.

**B. 300 random points of [0.5,1.75] x [0.125,1] x [0.05,1.75] x [-0.95,0.95], at 50 digits.** These are
not restricted to the region.

| check | max error |
|---|---|
| sum A_i = 0 | 5.3e-51 |
| sum A_i q_i = 0 | 5.3e-51 |
| **beta R = a** (R the 6x8 rigidity matrix: row ij has e_ij^T in the columns of body i and -e_ij^T in those of body j) | **5.5e-50** |
| **tau R = 0** | **2.1e-50** |
| closed-form R_ij and A_i against the coordinates and triangle areas | 7.8e-51 |
| scaled (t = 1) D, tau, beta, F, and trS(y) for random y, against literal | <= 1.1e-48 |

The sign pattern A1, A3 > 0 > A2, A4 held at every point.

**C. Chart formulas against the literal direct formulas, at random points with b > 0.**
- Sampling: 300 points per b-range, with al in [-2, 1], ga in [0.2, 3], xi in [-2.5, 0.7] and |x| < 1.
- Compared: D_e, tau_e/b, beta_e, Ftil against F/b^2, the region quantities, and trS(Y/b) for random Y.

| b range | precision | max error |
|---|---|---|
| [1e-6, 1/8] | 50 digits | 2.2e-46 |
| [1e-30, 1e-20] | 200 digits | 5.2e-156 |
| b = 0 limit against b = 1e-40 (50 points, al in [-1, 0.4], ga in [0.6, 1.8], xi in [-1.2, -0.4]) | 50 digits | 9.5e-40 |

The 200-digit run was needed because the literal formulas cancel catastrophically at tiny b.

**D. Points of Z solved at 50 digits.** With 300 tries the script finds 68 points. At these points:
- beta R = a and tau R = 0 hold to 2.1e-50.
- min D_e is 1.9e-4 >= 0.
- max min_y trS is 0.5713.

## 7. Tests of the verifier output (not part of the proof)

**End-to-end sampling (`check_leaves.py`).** Points of Z are found from the literal SPEC formulas only:
- a float sign scan of F in 600 steps over x in (-1, 1) (chart: xi in [-6, 3] with |x| < 1), which does
  not use the verifier's region bounds;
- refinement by findroot at 50 digits;
- a check of the region with the literal formulas.

Parameter samples: 4000 uniform in the domain, plus random points inside the leaves with the largest
bounds (2000 of them at 3/4, 3000 at 19/32).

Each Z point found must pass three tests:
1. it lies in an accepted leaf;
2. x is in the leaf's x-piece and x - x_lin(p) is in E (tolerance 1e-13);
3. the literal trS(y_leaf) is <= ub_leaf (tolerance 1e-12, far below the smallest margin in the table).

| leaves | Z points tested | failures | largest trS(y_leaf) | smallest margin ub - trS |
|---|---|---|---|---|
| direct, 3/4 | 2760 | 0 | 0.5907 | 0.079 |
| chart, 3/4 | 2711 | 0 | 0.5846 | 0.067 |
| direct, 19/32 | 3216 | 0 | 0.592634 (leaf bound 0.593740) | 0.00109 |
| chart, 19/32 | 3792 | 0 | 0.5847 | 0.0068 |

**Negative control (`negative_control.py`).**
- Setup: the same branch-and-bound on one direct and one chart initial cell, each containing a maximiser
  found by `explore.py`. The thresholds are below the values of min_y trS there: 9/16 and 37/64 in the
  direct cell, 9/16 in the chart cell.
- The verifier fails, as it must: 1773 unproved leaves in direct cell (2,3,8) and 1323 in chart cell
  (7,11,8), within a 4000-node budget.
- The unproved leaves really sit where min_y trS exceeds the threshold. Take the unproved leaves that
  have a root above their centre (1445 of the 1773 in the direct cell, 778 of the 1323 in the chart
  cell), and evaluate min_y trS at that root:
  - direct cell: 1445 of 1445 exceed 9/16, and 1444 of 1445 exceed 37/64;
  - chart cell: 778 of 778 exceed 9/16.

**Float exploration (`explore.py`, 60^3 parameter grids; grid points only, not a proof).**
- **One root per point.** Take every direct grid point with a nonempty region x-interval (45 652 of
  them). At each one, F has exactly one sign change on a 60-point grid of that interval. The same holds
  at all 46 295 such chart grid points.
- **Why there is always a root.** A short computation shows that, when the other inequalities hold:
  - F <= 0 on each of the four lower faces of the region (R13 = R12, R24 = R12, R14 = R34, R23 = R34);
  - F >= 0 on each of the two upper faces (R12 = R14, R12 = R23).

  For example, on R13 = R12 only the second product of F survives:
  F = -(r12^3 - r14^3)(r24^3 - r34^3)(r13^3 - r23^3) <= 0.
- **Consequence.** Above every point of the projection of the region there is a point of Z. At the lower
  end of the region x-interval a lower face is active, and at the upper end an upper face is active, so
  F changes sign across the interval.
  - Hence F can never remove a box that meets the projection. Only boxes outside it can be discarded.
  - The region tests caught all of those, which is why the "F != 0" column is 0.
- **Ranges found.** Direct: x in [0, 0.434]. Chart: al in [-1.69, 0.54], ga in [0.46, 2.70],
  xi in [-2.20, -0.40]. These agree with the SPEC's al and ga ranges.

**Maximum of min_y trS on Z (`maxsearch.py`).**
- **Method.** At each parameter point, x (xi) is the root of F in the region x-interval, found by
  bisection, and min_y trS is the explicit minimum of the quadratic.
  - Compass search with halving steps, clamped to the domain (chart: b >= 1e-12).
  - Started from the 25 best points of a 28^3 grid, and from the maximisers reported by `explore.py` and
    `check_leaves.py`.
- **Direct: about 0.59270.**
  - The largest value found is 0.5926966, at the kite a = c = 0.57735 (about 1/sqrt 3), b = 0.32995,
    x = 0. It lies on the face a = c, where the region's x-interval is the single point x = 0.
  - A second family of near-maxima, up to 0.5926915 at (a, b, c, x) = (1.7063, 0.9702, 0.5628, 0.0087),
    runs toward the face b = 1, where again x = 0.
  - F vanishes identically on both faces: for a = c (x = 0), R12 = R14 and R23 = R34; for b = 1 (x = 0),
    R12 = R23 and R14 = R34.
- **Chart: 0.58486**, at b = 1/8. This family continues into the direct domain. At the search's lower
  limit b = 1e-12, the local maxima found are 0.583333, at ga = 1 and xi = -0.75, for several values of al.
- The 60^3 exploration grid (max 0.5917) missed the kite maximum. The sampling test independently found
  0.592634.

## 8. Problems, discrepancies with the SPEC, remarks

- **No formula in the SPEC is wrong.** Every check of section 6 passed, including beta R = a and tau R = 0.
- **The maximum agrees with the SPEC's estimate.** The SPEC says "at most about 0.59". The numerical
  maximum of min_y trS on Z is 0.5927, which is below the certified 19/32 = 0.59375.
- **The literal formulas lose all precision near b = 0**, as the SPEC anticipates. They serve only as a
  reference. The chart comparison at b around 1e-25 needed 200 digits.
- **"A2h nonzero" holds only on the relevant set.** The SPEC says A1h and A2h are nonzero, but
  A2h = -(1+b) ga s/2 vanishes at ga = 0, which lies in OMEGA_chart. The region keeps ga >= 1/2 at b = 0,
  so it holds on the relevant set. Likewise OMEGA_direct contains c = 0 (where A2 = 0), but region points
  have c(2a + c) >= b, so c >= 0.035. Also, x = +-1 (where s = 0) is impossible in the region.
  - Boxes near these places are either discarded by the region tests or split until they are.
  - A division by an interval containing 0 there would only have produced NaN, which counts as
    "not proved".
- **The first root enclosure was too loose.**
  - The first version used a plain interval Newton step in x, with the numerator F(B, x_m) evaluated
    naively.
  - It was correct, but because of the dependency effect in F it stalled about 50 times wider than the
    true spread of the root over a box. It needed roughly 14 times as many leaves on a set of test cells.
  - Using the mean-value numerator, and bounding trS in sheared coordinates, fixed this.
  - Only efficiency was affected, not correctness.
- **Bugs fixed during development.** None of them affected a certified result.
  - A wrong return shape in `enclose_x` crashed the code.
  - A JSON key clash crashed a run after it finished; the run was repeated.
  - The logged double of a bound could round down; it is now rounded up (section 1).
  - `check_leaves.py` first scanned the chart's dependent coordinate over a placeholder range. It was
    replaced by `dep_scan_range` (section 7) before any of the reported tests.
- **Degenerate faces of the region.**
  - On the faces b = 1 (x = 0) and a = c (x = 0), F vanishes identically and the region's x-interval
    collapses to a point. F also vanishes on R12 = R13 = R14.
  - None of this caused failures: the linear x-range pins x, and every formula stays finite there.

## 9. Reproducing

The commands that write the logs of this directory, with their running times, are in the section on
the second implementation in the top-level `README.md`.
