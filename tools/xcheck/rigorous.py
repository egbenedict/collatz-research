"""Rigorous interval arithmetic over exact rationals (no floating point anywhere).

Every real constant or function value is returned as an interval [lo, hi] of Fractions
that is guaranteed to contain the true value. The series use integer fixed point with
floor for lower bounds, ceil for upper bounds, and explicit tail bounds. Results are
rounded outward to a grid of 2^-PBITS to keep numerators small.

    LN2, LN3, DELTA         enclosures of ln 2, ln 3, log2 3
    ln(r), log2(r)          r > 0 rational -> interval
    exp(r), pow2(r)         r >= 0 rational -> interval
    I                       interval type with + - * / (sign-aware)
"""
from __future__ import annotations

from fractions import Fraction as F
from math import floor, ceil

PBITS = 320
ONE = 1 << PBITS


def rdown(x: F) -> F:
    return F(floor(x * ONE), ONE)


def rup(x: F) -> F:
    return F(ceil(x * ONE), ONE)


class I:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        lo = F(lo)
        hi = lo if hi is None else F(hi)
        if lo > hi:
            raise ValueError(f"empty interval [{lo}, {hi}]")
        self.lo, self.hi = lo, hi

    @staticmethod
    def _c(x) -> "I":
        return x if isinstance(x, I) else I(x)

    def __add__(self, o):
        o = I._c(o)
        return I(rdown(self.lo + o.lo), rup(self.hi + o.hi))
    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, o):
        return self + (-I._c(o))

    def __rsub__(self, o):
        return I._c(o) - self

    def __mul__(self, o):
        o = I._c(o)
        ps = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return I(rdown(min(ps)), rup(max(ps)))
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = I._c(o)
        if o.lo <= 0 <= o.hi:
            raise ZeroDivisionError("interval division by an interval containing 0")
        return self * I(rdown(1 / o.hi), rup(1 / o.lo)) if o.lo > 0 else self * I(rdown(1 / o.hi), rup(1 / o.lo))

    def __rtruediv__(self, o):
        return I._c(o) / self

    def __pow__(self, n: int):
        assert n >= 0 and self.lo >= 0
        r = I(1)
        b = self
        while n:
            if n & 1:
                r = r * b
            b = b * b
            n >>= 1
        return r

    def __repr__(self):
        return f"I[{float(self.lo):.17g}, {float(self.hi):.17g}]"

    @property
    def width(self):
        return self.hi - self.lo


# ----------------------------------------------------------------------------
# Series in integer fixed point (scale 2^B, B > PBITS so final rounding dominates)
# ----------------------------------------------------------------------------
B = PBITS + 64
SCALE = 1 << B


def _atanh_inv(k: int) -> tuple[int, int]:
    """atanh(1/k) = sum_{n>=0} 1/((2n+1) k^{2n+1}), as fixed-point (lo, hi), k >= 2."""
    lo = hi = 0
    kp = k                     # k^{2n+1}
    n = 0
    while True:
        d = (2 * n + 1) * kp
        lo += SCALE // d
        hi += -(-SCALE // d)
        if SCALE // d == 0:
            # tail: sum_{j>n} 1/((2j+1) k^{2j+1}) < 1/(k^{2n+3}) * 1/(1 - 1/k^2) <= 2 / k^{2n+3}
            hi += 2
            return lo, hi
        kp *= k * k
        n += 1


def _fx(lo: int, hi: int) -> I:
    return I(rdown(F(lo, SCALE)), rup(F(hi, SCALE)))


_a3 = _atanh_inv(3)
_a5 = _atanh_inv(5)
LN2 = _fx(2 * _a3[0], 2 * _a3[1])                       # ln 2 = 2 atanh(1/3)
LN3 = LN2 + _fx(2 * _a5[0], 2 * _a5[1])                 # ln 3 = ln 2 + 2 atanh(1/5)
DELTA = LN3 / LN2


def _ln_mantissa(t_lo: int, t_hi: int) -> tuple[int, int]:
    """2*atanh(t) for fixed-point t in [0, 1/3): returns fixed-point (lo, hi)."""
    t2_lo = (t_lo * t_lo) >> B
    t2_hi = -((-(t_hi * t_hi)) >> B)
    p_lo, p_hi = t_lo, t_hi
    s_lo = s_hi = 0
    n = 0
    while True:
        s_lo += p_lo // (2 * n + 1)
        s_hi += -(-p_hi // (2 * n + 1))
        p_lo = (p_lo * t2_lo) >> B
        p_hi = -((-(p_hi * t2_hi)) >> B)
        n += 1
        if p_hi < 4:
            # remaining sum <= p_hi/(2n+1) * 1/(1-t^2) <= 2 p_hi (t < 1/3), plus rounding slack
            s_hi += 2 * p_hi + n + 2
            return 2 * s_lo, 2 * s_hi


def _ln_point(r: F) -> I:
    if r <= 0:
        raise ValueError("ln of non-positive")
    e = r.numerator.bit_length() - r.denominator.bit_length()
    s = r / (F(2) ** e)
    if s < 1:
        e -= 1
        s *= 2
    elif s >= 2:
        e += 1
        s /= 2
    assert 1 <= s < 2
    t = (s - 1) / (s + 1)                         # in [0, 1/3)
    t_lo, t_hi = floor(t * SCALE), ceil(t * SCALE)
    m_lo, m_hi = _ln_mantissa(t_lo, t_hi)
    return LN2 * e + _fx(m_lo, m_hi)


def ln(x) -> I:
    x = I._c(x)
    a = _ln_point(x.lo)
    b = a if x.hi == x.lo else _ln_point(x.hi)
    return I(a.lo, b.hi)


def log2(x) -> I:
    return ln(x) / LN2


def _exp_small(r: F) -> I:
    """exp(r) for 0 <= r < 1 by Taylor series in fixed point."""
    assert 0 <= r < 1
    x_lo, x_hi = floor(r * SCALE), ceil(r * SCALE)
    s_lo = s_hi = SCALE                 # term 0 = 1
    t_lo, t_hi = SCALE, SCALE
    n = 1
    while True:
        t_lo = (t_lo * x_lo >> B) // n
        t_hi = -((-(t_hi * x_hi)) >> B)
        t_hi = -(-t_hi // n)
        s_lo += t_lo
        s_hi += t_hi
        n += 1
        if t_hi < 4:
            s_hi += 2 * t_hi + n + 2       # tail <= 2 * last term (x < 1, n >= 2)
            return _fx(s_lo, s_hi)


def _exp_point(r: F) -> I:
    """exp(r) for r >= 0 via exp(r) = 2^q * exp(r - q ln 2) with 0 <= r - q ln 2 < 1."""
    if r < 0:
        return 1 / _exp_point(-r)
    q = floor(r / LN2.hi)
    # f = r - q ln2 lies in [r - q*LN2.hi, r - q*LN2.lo]
    f_lo, f_hi = r - q * LN2.hi, r - q * LN2.lo
    while f_hi >= 1:
        q += 1
        f_lo, f_hi = r - q * LN2.hi, r - q * LN2.lo
    f_lo = max(f_lo, F(0))
    a, b = _exp_small(f_lo), _exp_small(f_hi)
    two_q = F(2) ** q
    return I(rdown(a.lo * two_q), rup(b.hi * two_q))


def exp(x) -> I:
    x = I._c(x)
    a = _exp_point(x.lo)
    b = a if x.hi == x.lo else _exp_point(x.hi)
    return I(a.lo, b.hi)


def pow2(x) -> I:
    """2^x for x >= 0, via exact 2^floor(x) times exp(frac * ln 2)."""
    x = I._c(x)

    def pt(r: F, upper: bool) -> F:
        n = floor(r)
        f = r - n
        e = _exp_point(f * (LN2.hi if upper else LN2.lo)) if f else I(1)
        v = e.hi if upper else e.lo
        return rup(v * F(2) ** n) if upper else rdown(v * F(2) ** n)

    return I(pt(x.lo, False), pt(x.hi, True))


def ceil_log2_int(n: int) -> int:
    return (n - 1).bit_length()


if __name__ == "__main__":
    from decimal import Decimal, getcontext
    getcontext().prec = 110
    ln2, ln3 = Decimal(2).ln(), Decimal(3).ln()
    d = ln3 / ln2

    def inside(iv: I, x: Decimal, name: str):
        lo, hi = Decimal(iv.lo.numerator) / iv.lo.denominator, Decimal(iv.hi.numerator) / iv.hi.denominator
        slack = abs(x) * Decimal(10) ** -105          # the Decimal reference's own rounding
        ok = lo - slack <= x <= hi + slack
        print(f"{'ok ' if ok else 'BAD'} {name}: width {float(iv.width):.2e}")
        assert ok
    inside(LN2, ln2, "ln 2")
    inside(LN3, ln3, "ln 3")
    inside(DELTA, d, "log2 3")
    for r in [F(1, 7), F(355, 113), F(10**30 + 7, 3), F(2**71 + 1)]:
        inside(ln(r), Decimal(r.numerator).ln() - Decimal(r.denominator).ln(), f"ln({float(r):.4g})")
    for r in [F(0), F(1, 3), F(5, 2), F(123, 4), F(-7, 3)]:
        inside(exp(r), (Decimal(r.numerator) / r.denominator).exp(), f"exp({float(r):.4g})")
    for r in [F(71), F(7302, 100), F(931, 7)]:
        inside(pow2(r), (Decimal(r.numerator) / r.denominator * ln2).exp(), f"2^{float(r):.4g}")
    print("all enclosures contain the 110-digit Decimal values")
