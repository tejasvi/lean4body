# Independent replication: a computer-assisted inequality on convex four-body central configurations

You are asked to write, FROM THIS SPECIFICATION ONLY, an independent rigorous verifier (interval arithmetic,
outward rounding) of the inequality below, and run it.  Do not read any code outside ~/convex4/work/indep/ (in
particular not ~/convex4/work/rig/ or ~/convex4/work/*.py): the point is an implementation that shares no code and
no derivation shortcuts with the original one.  Python: ~/sundman/venv/bin/python (has python-flint 0.9 = arb,
mpmath, numpy, sympy).  Never run python with cwd = /tmp.  Use absolute paths.  The machine has 32 cores, shared with other
jobs; use at most 16 worker processes.  Log long runs to files in ~/convex4/work/indep/ and run them in the background.

## 1. Configuration and quantities (the "direct" coordinates)

Parameters a, b, c > 0 and x in (-1, 1), s = sqrt(1 - x^2).  Positions
    q1 = (1, 0),  q2 = (a x, a s),  q3 = (-b, 0),  q4 = (-c x, -c s).
(The diagonals 13 and 24 meet at the origin at angle theta = arccos x; bodies are counterclockwise 1,2,3,4.)
Squared distances R_ij = |q_i - q_j|^2, distances r_ij, and s_ij = r_ij^(-3).  Edges e in {12,13,14,23,24,34}.

Region (closure of the CCR region C):
    R13 >= R12,  R24 >= R12,  R12 >= R14,  R12 >= R23,  R14 >= R34,  R23 >= R34.
Central-configuration condition (Dziobek / CCR equation):
    F = (r24^3 - r14^3)(r13^3 - r12^3)(r23^3 - r34^3) - (r12^3 - r14^3)(r24^3 - r34^3)(r13^3 - r23^3) = 0.
Let Z = { (a,b,c,x) : region inequalities hold and F = 0 }.

Signed areas (A_i = signed area of the triangle NOT containing body i; with [u v w] the counterclockwise-signed
area of the triangle u, v, w):
    A1 = [q2 q3 q4],  A2 = -[q1 q3 q4],  A3 = [q1 q2 q4],  A4 = -[q1 q2 q3].
They satisfy sum A_i = 0 and sum A_i q_i = 0 identically, and A1, A3 > 0 > A2, A4 for these configurations.

Dziobek multiplier and weights (valid as formulas everywhere, no use of F = 0):
    lam = (s12 s34 - s13 s24) / (s12 + s34 - s13 - s24),     w_e = s_e - lam.

Certificate data, for e = ij:
    D_e   = - w_e / (3 s_e A_i A_j)                                   (scalar)
    tau_e = A_i A_j r_e                                               (scalar; a self-stress of the framework K4)
    beta_e = - sum_{l not in {i,j}}  r_e N_{e,l} / (8 A_l) * q_l      (a vector in R^2)
        where for l not in e, k is the fourth index (k not in {i, j, l}) and N_{e,l} = R_ik + R_jk - R_ij.
    trS(y) = sum_e D_e |beta_e + y tau_e|^2,   y in R^2 a constant vector (independent of the point is allowed
             per box: you choose y per box, e.g. the minimiser of trS at the box centre).

THE CLAIM TO VERIFY:  for every point of Z with (a, b, c) in the box
    OMEGA_direct = [0.5, 1.75] x [1/8, 1] x [0, 1.75],   x in [-1, 1],
there is (for the sub-box containing it) a y with trS(y) < 3/4.   (It is known numerically that min_y trS is at
most about 0.59 on Z, so 3/4 leaves room.)

## 2. The collision corner: blow-up chart for small b

As b -> 0 the points of Z accumulate at (a, b, c, x) = (1, 0, 0, 1/2) (an equilateral triangle with bodies 3 and 4
colliding), where the direct formulas degenerate (0/0).  For b in (0, 1/8] use the chart
    a = 1 + b al,   c = b ga,   x = 1/2 + b xi,     (b, al, ga) in OMEGA_chart = [0, 1/8] x [-2, 1] x [0, 3],
and xi in whatever range the region inequalities allow.  In this chart:
  * R13 - R12, R24 - R12, R12 - R14, R12 - R23 are each b times an expression that stays O(1) as b -> 0;
    R14 - R34 and R23 - R34 tend to 1 (r34 = b rho with rho^2 = 1 - 2 ga x + ga^2).
  * A1 = b A1h, A2 = b A2h with A1h, A2h = O(1) nonzero; A3, A4 = O(1).
  * s34 = b^-3 rho^-3 blows up; the other s_e are O(1).
  * trS(y) is a sum of terms D_e |beta_e + y tau_e|^2; after the substitution y = Y / b (Y the chart variable) it
    reads sum_e D_e |beta_e + Y (tau_e / b)|^2, and every D_e, tau_e / b and beta_e has a finite limit as b -> 0 (and the whole expression is
    analytic in (b, al, ga, xi) up to b = 0 on the relevant set).  You must derive these rescaled expressions
    yourself (sympy recommended: substitute, cancel the powers of b symbolically, and CHECK numerically at random
    points with b > 0 that the rescaled formulas reproduce the direct ones to high precision), then evaluate the
    rescaled formulas in interval arithmetic on boxes that contain b = 0.
  * F = 0 must likewise be replaced by an equivalent rescaled equation that stays nondegenerate at b = 0
    (divide out the appropriate power of b and any factor that is positive on the region).  Justify the
    equivalence in your write-up.
THE CLAIM IN THE CHART:  for every point of Z with b in (0, 1/8] (these satisfy al in (-2, 1), ga in [0, 17/6);
you may take that as given), trS(Y/b) < 3/4 for the Y chosen for its sub-box, with the check done in the rescaled
(b = 0 inclusive) formulas on OMEGA_chart.

## 3. What your verifier must do

A branch-and-bound over OMEGA_direct and over OMEGA_chart (3-dimensional boxes in the parameters; the fourth
coordinate x or xi is dependent and must be enclosed rigorously), with python-flint arb (all arithmetic in arb;
convert floats to arb exactly; never let a float rounding decide anything).  For each box B:
  (i)  enclose every x (or xi) such that some point of B x {x} lies in the region with F = 0 (for example an
       interval Newton / Krawczyk contraction on F in x with the other variables in B, plus bisection; points
       with any region inequality strictly violated on the whole sub-box can be discarded);
 (ii)  if nothing is left, the box is discarded;
(iii)  else choose y (or Y) and prove trS(y) < 3/4 on B x (each piece of the x-enclosure) intersected with the
       region, using a first-order Taylor / mean-value form or any other rigorous enclosure;
 (iv)  otherwise bisect B and recurse.
Any NaN, infinity, or division by an interval containing 0 must count as "not proved" (never as success).
Choose the initial grids and splitting rules yourself.  Independent choices (F instead of any rewritten form,
your own derivative scheme, your own y choice) are welcome; that is the point.

Deliverables in ~/convex4/work/indep/: the code, a README describing the method and the exact statement proved,
the run logs (boxes accepted / discarded, max certified upper bound of trS, time), and a short report of
anything that did not work or any discrepancy you found with this specification (for example if some formula
above is wrong, the numerical checks in section 2 will show it: say so rather than silently fixing it).
