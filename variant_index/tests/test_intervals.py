import random
import unittest

from helpers import ROOT  # noqa: F401  (puts the package on the path)
from vindex.intervals import IntervalStore


def brute(starts, ends, qs, qe):
    return [i for i in range(len(starts)) if starts[i] <= qe and ends[i] >= qs]


class TestIntervalStore(unittest.TestCase):
    def test_empty_store(self):
        s = IntervalStore([], [])
        self.assertEqual(s.overlapping(1, 100), [])
        self.assertFalse(s.any_overlap(1, 100))

    def test_matches_brute_force_on_random_layered_intervals(self):
        rng = random.Random(7)
        for trial in range(60):
            n = rng.randint(1, 400)
            starts, ends = [], []
            for _ in range(n):
                a = rng.randint(1, 5000)
                # many short intervals and a few very long ones: the case that breaks naive early stopping
                length = rng.choice([0, 0, 1, 3, 10, 40, 400, 4000]) if rng.random() < 0.97 else rng.randint(2000, 6000)
                starts.append(a)
                ends.append(a + length)
            s = IntervalStore(starts, ends)
            for _ in range(200):
                qs = rng.randint(0, 5200)
                qe = qs + rng.choice([0, 0, 1, 5, 50, 500])
                self.assertEqual(s.overlapping(qs, qe), brute(starts, ends, qs, qe), (trial, qs, qe))
                self.assertEqual(s.any_overlap(qs, qe), bool(brute(starts, ends, qs, qe)))

    def test_intervals_of_every_length_class_are_found_together(self):
        # the store keeps up to 10 kb, up to 1 Mb and longer intervals apart (a long one no longer lengthens the walk of the short ones); the answer must not change
        rng = random.Random(11)
        for trial in range(30):
            starts, ends = [], []
            for _ in range(rng.randint(1, 300)):
                a = rng.randint(1, 4_000_000)
                length = rng.choice([0, 1, 50, 5_000, 10_000, 10_001, 200_000, 1_000_000, 1_000_001, 3_000_000]) if rng.random() < 0.9 else rng.randint(0, 3_500_000)
                starts.append(a)
                ends.append(a + length)
            s = IntervalStore(starts, ends)
            self.assertLessEqual(len(s.parts), 3)
            for _ in range(150):
                qs = rng.randint(0, 4_200_000)
                qe = qs + rng.choice([0, 1, 50, 5_000, 300_000])
                self.assertEqual(s.overlapping(qs, qe), brute(starts, ends, qs, qe), (trial, qs, qe))
                self.assertEqual(s.any_overlap(qs, qe), bool(brute(starts, ends, qs, qe)))

    def test_one_very_long_interval_has_a_store_of_its_own_and_is_still_found(self):
        starts = list(range(1000, 2_001_000, 20))                                   # 100,000 short intervals
        ends = [a + 5 for a in starts]
        starts.append(10)
        ends.append(5_000_000)                                                          # and one that reaches past all of them
        s = IntervalStore(starts, ends)
        self.assertEqual(sorted(len(p.starts) for p in s.parts), [1, 100_000])         # a short interval's walk back no longer passes through the long one's reach
        self.assertEqual(s.overlapping(1_999_000, 1_999_010), brute(starts, ends, 1_999_000, 1_999_010))
        self.assertEqual(s.overlapping(4_900_000, 4_900_001), [len(starts) - 1])

    def test_boundaries_are_inclusive(self):
        s = IntervalStore([10, 20], [15, 25])
        self.assertEqual(s.overlapping(15, 15), [0])
        self.assertEqual(s.overlapping(16, 19), [])
        self.assertEqual(s.overlapping(25, 30), [1])
        self.assertEqual(s.overlapping(1, 10), [0])

    def test_equal_starts_keep_input_positions(self):
        s = IntervalStore([5, 5, 5], [6, 100, 5])
        self.assertEqual(s.overlapping(50, 60), [1])
        self.assertEqual(s.overlapping(5, 5), [0, 1, 2])


if __name__ == '__main__':
    unittest.main()
