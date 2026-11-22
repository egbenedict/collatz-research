"""Known facts the toolkit must reproduce. Run: python3 -m unittest -v (from tools/)."""
import unittest

from collatz import (COLLATZ, Map, canonical_cycle, check_bijection, find_cycles,
                     orbit, parity_vector, residue_for_parity, standard_step,
                     stopping_time, total_stopping_time)


class KnownFacts(unittest.TestCase):
    def test_27_standard(self):
        # 27 takes 111 steps under C(n) and peaks at 9232.
        x, steps, peak = 27, 0, 27
        while x != 1:
            x = standard_step(x)
            steps += 1
            peak = max(peak, x)
        self.assertEqual((steps, peak), (111, 9232))

    def test_27_accelerated(self):
        # Under T, odd steps are merged with the following halving: 111 - 41 odd steps = 70.
        self.assertEqual(total_stopping_time(27), 70)

    def test_stopping_time_small(self):
        self.assertEqual(stopping_time(2), 1)
        self.assertEqual(stopping_time(3), 4)   # 3 → 5 → 8 → 4 → 2
        self.assertEqual(stopping_time(7), 7)

    def test_positive_cycles_3x1(self):
        cyc = {c for c in find_cycles(20) if c[0] > 0}
        self.assertEqual(cyc, {(1, 2)})

    def test_negative_cycles_3x1(self):
        cyc = {c for c in find_cycles(20) if c[0] < 0}
        self.assertIn((-1,), cyc)
        self.assertIn((-5, -7, -10), cyc)
        self.assertTrue(any(c[0] == -17 and len(c) == 11 for c in cyc))
        self.assertEqual(len(cyc), 3)

    def test_5x1_cycles(self):
        cyc = {c for c in find_cycles(16, Map(5, 1)) if c[0] > 0}
        self.assertIn(canonical_cycle([1, 3, 8, 4, 2]), cyc)
        self.assertIn(canonical_cycle([13, 33, 83, 208, 104, 52, 26]), cyc)
        self.assertIn(canonical_cycle([17, 43, 108, 54, 27, 68, 34]), cyc)

    def test_orbit_finds_cycle(self):
        o = orbit(27)
        self.assertEqual(canonical_cycle(o.cycle), (1, 2))
        self.assertFalse(o.bailed)

    def test_bijection(self):
        for k in range(1, 11):
            self.assertTrue(check_bijection(k, COLLATZ))
            self.assertTrue(check_bijection(k, Map(5, 1)))

    def test_residue_roundtrip(self):
        for n in [1, 7, 27, 97, 871, 12345]:
            v = parity_vector(n, 12)
            self.assertEqual(residue_for_parity(v), n % (1 << 12))


if __name__ == "__main__":
    unittest.main()
