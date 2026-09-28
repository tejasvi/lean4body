"""First-order jets over flint.arb: value + gradient, for mean-value (centred) enclosures.
If every input is J(box_i, e_i), the output's d[i] encloses the partial derivative over the whole box."""
from flint import arb


class J:
    __slots__ = ('v', 'd')

    def __init__(self, v, d):
        self.v = v
        self.d = d

    def __add__(self, o):
        if isinstance(o, J):
            return J(self.v + o.v, [p + q for p, q in zip(self.d, o.d)])
        return J(self.v + o, self.d)

    __radd__ = __add__

    def __neg__(self):
        return J(-self.v, [-p for p in self.d])

    def __sub__(self, o):
        if isinstance(o, J):
            return J(self.v - o.v, [p - q for p, q in zip(self.d, o.d)])
        return J(self.v - o, self.d)

    def __rsub__(self, o):
        return J(o - self.v, [-p for p in self.d])

    def __mul__(self, o):
        if isinstance(o, J):
            return J(self.v * o.v, [self.v * q + o.v * p for p, q in zip(self.d, o.d)])
        return J(self.v * o, [p * o for p in self.d])

    __rmul__ = __mul__

    def __truediv__(self, o):
        if isinstance(o, J):
            v = self.v / o.v
            return J(v, [(p - v * q) / o.v for p, q in zip(self.d, o.d)])
        return J(self.v / o, [p / o for p in self.d])

    def __rtruediv__(self, o):
        v = o / self.v
        return J(v, [-v * p / self.v for p in self.d])

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 1
        if n == 1:
            return self
        vn1 = self.v ** (n - 1)
        return J(vn1 * self.v, [n * vn1 * p for p in self.d])

    def sqrt(self):
        s = self.v.sqrt()
        return J(s, [p / (2 * s) for p in self.d])


def seeds(boxes):
    n = len(boxes)
    zero = arb(0)
    return [J(boxes[i], [arb(1) if k == i else zero for k in range(n)]) for i in range(n)]


def mv(fm, fj, rads):
    """mean-value enclosure: f(mid) + sum_i f_i(box) [-r_i, r_i], intersected with the naive enclosure fj.v"""
    out = fm
    for di, ri in zip(fj.d, rads):
        out = out + di * arb(0, ri)
    try:
        return out.intersection(fj.v)
    except Exception:
        return out
