"""Linear simulation regression tests (nonlin = 0).

The base case is the DIII-D Alfven-eigenmode case of Models/DIIID at reduced
radial resolution (mj=200): 8 dynamic modes (n=1, m=1..4 and their sine
partners), 200 steps of dt=2. After ~400 Alfven times the most unstable
eigenmode dominates and lincheck reports its growth rate and frequency.

Checks, from most to least physical:
  * the reported (gamma, omega) is a single eigenmode: every variable and
    every poloidal harmonic agree;
  * the solution is exactly linear in the seed amplitude;
  * the time integration converges at second order in dt (Crank-Nicolson);
  * the result is independent of the number of MPI processes;
  * gamma, omega and the final eigenfunction match stored reference values.
"""

import unittest

import numpy as np

import far3d_harness as h


def dominant(growth):
    """(gamma, omega) of each (m, n) harmonic, from the psi lines."""
    return {(g["m"], g["n"]): (g["gamma"], g["omega"]) for g in growth if g["var"] == "psi"}


class LinearBaseCase(unittest.TestCase):

    reference = h.Reference("linear")

    @classmethod
    def setUpClass(cls):
        cls.sim = h.run_case(h.base_case(), "lin_base")
        cls.growth = h.read_growth_rates(cls.sim.farprt)
        cls.final = cls.sim.final_dump()

    @classmethod
    def tearDownClass(cls):
        cls.reference.save()

    def test_lincheck_reported_all_modes(self):
        # 4 harmonics (m=1..4, n=1) x 6 variables
        self.assertEqual(len(self.growth), 24)

    def test_single_eigenmode(self):
        # The strong harmonics (m=1..3) agree to <1% in gamma and <2% in
        # omega; the weakest one (m=4) is still converging onto the eigenmode
        # at t=400, so it only has to agree in gamma to 3%.
        def spread(vals):
            vals = np.array(vals)
            return (vals.max() - vals.min()) / np.median(vals)
        strong = [g for g in self.growth if g["m"] <= 3]
        self.assertGreater(min(g["gamma"] for g in self.growth), 0.0, "mode should be unstable")
        self.assertLess(spread([g["gamma"] for g in strong]), 0.01)
        self.assertLess(spread([g["omega"] for g in strong]), 0.02)
        self.assertLess(spread([g["gamma"] for g in self.growth]), 0.03)

    def test_growth_rate_matches_reference(self):
        # lincheck prints 6 significant digits
        key = sorted(self.growth, key=lambda g: (g["var"], g["m"]))
        self.reference.check(self, "gamma", [g["gamma"] for g in key], rtol=5e-5)
        self.reference.check(self, "omega", [g["omega"] for g in key], rtol=5e-5)

    def test_eigenfunction_matches_reference(self):
        # final per-harmonic amplitude profile of each field, at sample radii
        for k in h.FIELD_NAMES:
            self.reference.check(self, "field." + k, self.final[k][:8, ::20], rtol=1e-8)

    def test_only_n1_harmonics_excited(self):
        n0 = self.final["nn"] == 0
        for k in h.FIELD_NAMES:
            self.assertEqual(np.abs(self.final[k][n0]).max(), 0.0, k)

    def test_boundary_conditions(self):
        # psi, phi, pr vanish at the wall; all non-(0,0) harmonics vanish on axis
        for k in ("psi", "phi", "pr"):
            self.assertEqual(np.abs(self.final[k][:, -1]).max(), 0.0, k)
        for k in h.FIELD_NAMES:
            self.assertEqual(np.abs(self.final[k][:, 0]).max(), 0.0, k)

    def test_linearity_in_seed_amplitude(self):
        doubled = h.run_case(h.base_case(widthi="2.e-140,28*0.0"), "lin_doubled").final_dump()
        for k in h.FIELD_NAMES:
            self.assertLess(h.rel_diff(doubled[k], 2.0 * self.final[k]), 1e-13, k)


class MpiInvariance(unittest.TestCase):

    def test_process_count_does_not_change_result(self):
        one = h.run_case(h.base_case(maxstp=50), "lin_np1", nproc=1).final_dump()
        for n in (2, 3):
            many = h.run_case(h.base_case(maxstp=50), "lin_np%d" % n, nproc=n).final_dump()
            for k in h.FIELD_NAMES:
                self.assertLess(h.rel_diff(many[k], one[k]), 1e-12, "np=%d %s" % (n, k))


class TimeStepConvergence(unittest.TestCase):
    """The same run to t=400 with dt = 4, 2, 1, 0.5. Crank-Nicolson is second
    order, so the change in the final fields should shrink by 4 each time dt
    halves (measured: 3.96-4.23). The lincheck growth rates are not used
    here: their finite-time estimate for the weak harmonics is not monotone
    in dt, and farprt prints only 6 digits."""

    def test_second_order_in_dt(self):
        finals = []
        for dt, nsteps in ((4.0, 100), (2.0, 200), (1.0, 400), (0.5, 800)):
            finals.append(h.run_case(h.base_case(dt0=dt, maxstp=nsteps), "lin_dt%g" % dt).final_dump())
        for k in h.FIELD_NAMES:
            diffs = [h.rel_diff(finals[i + 1][k], finals[i][k]) for i in range(3)]
            for ratio in (diffs[0] / diffs[1], diffs[1] / diffs[2]):
                self.assertGreater(ratio, 3.5, "%s: successive differences %s" % (k, diffs))
                self.assertLess(ratio, 4.5, "%s: successive differences %s" % (k, diffs))


if __name__ == "__main__":
    unittest.main()
