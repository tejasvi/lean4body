"""Number-type helpers shared by the checks and the rigorous verifier.

* SQRT(x): square root that dispatches on the type of x (float, mpmath, arb, AD).
  For arb it first intersects the ball with [0, inf): every quantity whose square root is
  taken in this project (squared distances, 1 - x^2, rho^2) is mathematically >= 0 at every
  real point of the domain, so this is a valid enclosure.  If the resulting ball touches 0,
  later divisions produce NaN, which the verifier treats as "not proved".
* AD: forward-mode automatic differentiation whose value and partial derivatives are arb
  balls.  Evaluating a formula on AD inputs whose values are boxes gives rigorous enclosures
  of the function AND of its gradient over the whole box (every arb operation is rigorous;
  the chain rule is exact).
"""
import math
import mpmath
import numpy as np
from flint import arb


def arb_sqrt(x):
    if x.is_finite() and not (x.lower() >= 0):
        x = x.nonnegative_part()
    return x.sqrt()


class AD:
    __slots__ = ('v', 'g')

    def __init__(self, v, g):
        self.v = v
        self.g = g

    @staticmethod
    def var(v, k, n):
        """independent variable number k (of n) with value v (an arb ball)"""
        return AD(v, [arb(1) if i == k else arb(0) for i in range(n)])

    def __add__(self, o):
        if isinstance(o, AD):
            return AD(self.v + o.v, [x + y for x, y in zip(self.g, o.g)])
        return AD(self.v + o, self.g)

    __radd__ = __add__

    def __sub__(self, o):
        if isinstance(o, AD):
            return AD(self.v - o.v, [x - y for x, y in zip(self.g, o.g)])
        return AD(self.v - o, self.g)

    def __rsub__(self, o):
        return AD(o - self.v, [-x for x in self.g])

    def __neg__(self):
        return AD(-self.v, [-x for x in self.g])

    def __mul__(self, o):
        if isinstance(o, AD):
            sv, ov = self.v, o.v
            return AD(sv * ov, [sv * y + ov * x for x, y in zip(self.g, o.g)])
        return AD(self.v * o, [x * o for x in self.g])

    __rmul__ = __mul__

    def __truediv__(self, o):
        if isinstance(o, AD):
            inv = 1 / o.v
            q = self.v * inv
            return AD(q, [(x - q * y) * inv for x, y in zip(self.g, o.g)])
        return AD(self.v / o, [x / o for x in self.g])

    def __rtruediv__(self, o):
        inv = 1 / self.v
        q = o * inv
        f = -q * inv
        return AD(q, [f * x for x in self.g])

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 1
        r = self
        for _ in range(n - 1):
            r = r * self
        return r

    def sqrt(self):
        r = arb_sqrt(self.v)
        f = 1 / (2 * r)
        return AD(r, [f * x for x in self.g])


def SQRT(x):
    if isinstance(x, AD):
        return x.sqrt()
    if isinstance(x, arb):
        return arb_sqrt(x)
    if isinstance(x, float):
        return math.sqrt(x)
    if isinstance(x, np.ndarray):
        return np.sqrt(x)
    return mpmath.sqrt(x)
