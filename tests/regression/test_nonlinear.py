"""Nonlinear simulation regression tests (nonlin = 1).

In the base deck m0dy = 0, and setmod drops every mode that is also an
equilibrium mode from the evolved set unless its index is <= m0dy. Every
n = 0 mode is an equilibrium mode here, so with m0dy = 0 nothing can be
driven by the n=1 x n=1 coupling and a "nonlinear" run is bitwise identical
to a linear one. These tests use m0dy = 21 so all 21 n = 0 modes evolve.

Checks:
  * at vanishing amplitude the nonlinear solver reproduces the linear one;
  * at small amplitude the n = 0 response is driven quadratically: doubling
    the seed doubles n = 1 and quadruples n = 0;
  * the result is independent of the number of MPI processes (this exercises
    the transposes between the radial and toroidal-family decompositions) --
    currently only approximately, see MpiInvariance;
  * a strongly nonlinear run matches stored reference values.
"""

import unittest

import numpy as np

import far3d_harness as h

NSTEPS = 40


def nl_case(amplitude, **overrides):
    case = h.base_case(nonlin=1, m0dy=21, maxstp=NSTEPS, widthi="%r,28*0.0" % amplitude)
    return case.set(**overrides)


def amplitudes(dump, k):
    """(max |n=1 part|, max |n=0 part|) of field k."""
    n0 = dump["nn"] == 0
    return np.abs(dump[k][~n0]).max(), np.abs(dump[k][n0]).max()


class SmallAmplitudeLimit(unittest.TestCase):

    def test_matches_linear_run(self):
        nl = h.run_case(nl_case(1e-140), "nl_tiny").final_dump()
        lin = h.run_case(nl_case(1e-140, nonlin=0), "nl_tiny_linear").final_dump()
        n1 = nl["nn"] != 0
        for k in h.FIELD_NAMES:
            self.assertLess(h.rel_diff(nl[k][n1], lin[k][n1]), 1e-12, k)
            # the quadratically driven n=0 part is ~1e-280, far below the n=1 part
            self.assertLess(np.abs(nl[k][~n1]).max(), 1e-200, k)


class QuadraticCoupling(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.a = h.run_case(nl_case(1e-8), "nl_amp1").final_dump()
        cls.b = h.run_case(nl_case(2e-8), "nl_amp2").final_dump()

    def test_n0_is_driven(self):
        for k in h.FIELD_NAMES:
            n1, n0 = amplitudes(self.a, k)
            self.assertGreater(n0, 1e-7 * n1, k)

    def test_n1_scales_linearly(self):
        for k in h.FIELD_NAMES:
            self.assertAlmostEqual(amplitudes(self.b, k)[0] / amplitudes(self.a, k)[0], 2.0, delta=1e-5, msg=k)

    def test_n0_scales_quadratically(self):
        for k in h.FIELD_NAMES:
            self.assertAlmostEqual(amplitudes(self.b, k)[1] / amplitudes(self.a, k)[1], 4.0, delta=1e-3, msg=k)


class MpiInvariance(unittest.TestCase):
    """The same nonlinear run on 1, 2 and 3 MPI processes. Linear runs are
    bitwise independent of the process count (test_linear.MpiInvariance), and
    the nonlinear terms are computed pointwise on overlapping radial slabs, so
    nonlinear runs should agree to roundoff too. They don't: np=2 and np=3
    agree with each other exactly but differ from np=1 by ~1e-8 relative after
    one step, growing to ~3e-6 after 20. The loose check guards against gross
    breakage; the tight one records the known discrepancy."""

    @classmethod
    def setUpClass(cls):
        cls.dumps = {n: h.run_case(nl_case(1e-4, maxstp=20), "nl_np%d" % n, nproc=n).final_dump()
                     for n in (1, 2, 3)}

    def test_multi_process_runs_agree_with_each_other(self):
        for k in h.FIELD_NAMES:
            self.assertLess(h.rel_diff(self.dumps[3][k], self.dumps[2][k]), 1e-12, k)

    def test_process_count_changes_result_only_slightly(self):
        for n in (2, 3):
            for k in h.FIELD_NAMES:
                self.assertLess(h.rel_diff(self.dumps[n][k], self.dumps[1][k]), 1e-4, "np=%d %s" % (n, k))

    @unittest.expectedFailure
    def test_process_count_does_not_change_result(self):
        # KNOWN ISSUE: single- vs multi-process nonlinear runs differ (see class docstring)
        for n in (2, 3):
            for k in h.FIELD_NAMES:
                self.assertLess(h.rel_diff(self.dumps[n][k], self.dumps[1][k]), 1e-10, "np=%d %s" % (n, k))


class StronglyNonlinear(unittest.TestCase):
    """Seed amplitude 1e-4: the n=0 response is comparable to n=1 by the end
    of the run and the n=1 part no longer scales linearly."""

    reference = h.Reference("nonlinear")

    @classmethod
    def setUpClass(cls):
        cls.sim = h.run_case(nl_case(1e-4, nprint=10), "nl_strong")
        cls.final = cls.sim.final_dump()

    @classmethod
    def tearDownClass(cls):
        cls.reference.save()

    def test_regime_is_nonlinear(self):
        n1, n0 = amplitudes(self.final, "phi")
        self.assertGreater(n0 / n1, 0.1)

    def test_fields_match_reference(self):
        for k in h.FIELD_NAMES:
            self.reference.check(self, "field." + k, self.final[k][:, ::20], rtol=1e-7, atol=0.0)

    def test_energy_history_matches_reference(self):
        # eke/eme: per-harmonic kinetic and magnetic energy every nprint steps
        for name in ("eke", "eme"):
            cols, table = h.read_table(self.sim.path("%s_%s" % (name, self.sim.numrun)))
            self.reference.check(self, name, table, rtol=1e-7)

    def test_all_finite(self):
        for k in h.FIELD_NAMES:
            self.assertTrue(np.all(np.isfinite(self.final[k])), k)


if __name__ == "__main__":
    unittest.main()
