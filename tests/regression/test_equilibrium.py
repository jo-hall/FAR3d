"""Equilibrium conversion regression tests.

The equilibrium reaches FAR3d in two conversion steps:

  1. VMEC wout_*.nc --(BOOZ_XFORM, 'far' mode)--> woutb   (Boozer-coordinate
     Fourier coefficients of the metric, |B|, profiles)
  2. woutb + Data.txt --(far3d.x: seteq/vmec)--> profiles and metric
     coefficients on FAR3d's radial grid, saved in every fs##### dump.

Step 1 is checked by regenerating woutb and comparing it with
Models/DIIID/woutb (skipped if xbooz_xform has not been built). Step 2 is
checked against stored reference values, against resolution independence
(profiles are evaluated pointwise, so two grids must agree where they share
points), and against basic invariants.
"""

import os
import shutil
import subprocess
import unittest

import numpy as np

import far3d_harness as h

SAMPLE_J = list(range(0, 201, 20))         # radial indices stored in the reference
METRICS = ["sqg", "bmod", "grr", "grt", "gtt", "bst", "omdr", "omdt"]


def xbooz_exe():
    exe = os.environ.get("FAR3D_XBOOZ") or str(h.REPO / "BOOZ_XFORM" / "build" / "xbooz_xform")
    return exe if os.path.isfile(exe) and os.access(exe, os.X_OK) else None


class BoozXform(unittest.TestCase):

    def test_regenerated_woutb_matches_model(self):
        exe = xbooz_exe()
        if exe is None:
            self.skipTest("xbooz_xform not built (cd BOOZ_XFORM && make), or set FAR3D_XBOOZ")
        eqdir = h.REPO / "BOOZ_XFORM" / "Equilibria" / "DIIID"
        rundir = h.WORKDIR / "eq_booz_xform"
        if rundir.exists():
            shutil.rmtree(rundir)
        rundir.mkdir(parents=True)
        for f in ("in_booz.d3d_159243_00805", "wout_d3d_159243_00805.nc"):
            shutil.copy(eqdir / f, rundir / f)
        proc = subprocess.run([exe, "in_booz.d3d_159243_00805", "far"], cwd=rundir,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=600)
        self.assertEqual(proc.returncode, 0, proc.stdout[-2000:])

        new = h._records(rundir / "woutb")
        ref = h._records(h.MODEL_DIR / "woutb")
        self.assertEqual(len(new), len(ref), "number of records in woutb changed")
        self.assertEqual(new[0], ref[0], "woutb header record (nfp, lbmax, ns, lasym) changed")
        for i, (a, b) in enumerate(zip(new, ref)):
            if a == b:
                continue
            # allow last-bit differences in floating point records
            self.assertEqual(len(a), len(b), "record %d length changed" % i)
            self.assertEqual(len(a) % 8, 0, "record %d differs and is not a float64 record" % i)
            x, y = np.frombuffer(a, ">f8"), np.frombuffer(b, ">f8")
            self.assertLess(h.rel_diff(x, y), 1e-10, "record %d differs" % i)


class MappedEquilibrium(unittest.TestCase):

    reference = h.Reference("equilibrium")

    @classmethod
    def setUpClass(cls):
        # maxstp=0: set up the equilibrium, write the step-0 dump, stop
        cls.d = h.run_case(h.base_case(maxstp=0), "eq_mj200").dump(0)
        cls.i00 = list(cls.d["mmeq"]).index(0)

    @classmethod
    def tearDownClass(cls):
        cls.reference.save()

    def test_profiles_match_reference(self):
        for k in h.PROFILE_NAMES:
            self.reference.check(self, "profile." + k, self.d[k][SAMPLE_J], rtol=1e-9)

    def test_metric_coefficients_match_reference(self):
        for k in METRICS:
            self.reference.check(self, "metric." + k, self.d[k][:, SAMPLE_J], rtol=1e-9, atol=1e-14)

    def test_resonant_surfaces_match_reference(self):
        # rs(l): radius of the q = m/n surface of each mode, used to seed the
        # perturbation
        self.reference.check(self, "rs", self.d["rs"], rtol=1e-9)

    def test_normalization_on_axis(self):
        # thermal pressure, density and temperatures are normalized to their
        # on-axis values; |B| to its on-axis value (only the (0,0) harmonic
        # of |B| is 1, and only to the accuracy of the Boozer fit)
        for k in ("preq", "denseq", "teeq", "tieq"):
            self.assertAlmostEqual(self.d[k][0], 1.0, places=12, msg=k)
        self.assertAlmostEqual(self.d["bmod"][self.i00, 0], 1.0, delta=1e-3)

    def test_jacobian_and_metric_positive(self):
        self.assertGreater(self.d["sqg"][self.i00].min(), 0.0)
        self.assertGreater(self.d["grr"][self.i00].min(), 0.0)
        self.assertGreater(self.d["gtt"][self.i00, 1:].min(), 0.0)

    def test_up_down_symmetry(self):
        # a stellarator-symmetric (lasym = F) tokamak equilibrium: the
        # (m,0) and (-m,0) harmonics of the cosine-type |B| coefficient carry
        # the cos and sin parts; the sin parts must vanish
        mm = list(self.d["mmeq"])
        for m in range(1, 11):
            self.assertLess(np.abs(self.d["bmod"][mm.index(-m)]).max(), 1e-12, "m=%d" % m)


class ResolutionIndependence(unittest.TestCase):
    """Profiles and metric coefficients are evaluated pointwise from the
    splined equilibrium, so a mj=400 grid must reproduce the mj=200 values
    at the shared radii. Quantities built with radial finite differences
    (bst, omd*) differ at the discretization level only."""

    @classmethod
    def setUpClass(cls):
        cls.a = h.run_case(h.base_case(maxstp=0), "eq_res200").dump(0)
        cls.b = h.run_case(h.base_case(maxstp=0, mj=400, edge_p=390), "eq_res400").dump(0)

    def test_grids_share_points(self):
        np.testing.assert_allclose(self.b["r"][::2], self.a["r"], rtol=0, atol=1e-15)

    def test_pointwise_quantities_identical(self):
        for k in h.PROFILE_NAMES:
            self.assertLess(h.rel_diff(self.b[k][::2], self.a[k]), 1e-12, k)
        for k in ("sqg", "bmod", "grr", "grt", "gtt"):
            self.assertLess(h.rel_diff(self.b[k][:, ::2], self.a[k]), 1e-12, k)

    def test_finite_difference_quantities_converge(self):
        for k in ("bst", "omdr"):
            self.assertLess(h.rel_diff(self.b[k][:, 2:-2:2], self.a[k][:, 1:-1]), 1e-2, k)


class SafetyFactor(unittest.TestCase):

    def test_q_and_iota_consistent_with_vmec_q(self):
        d = h.run_case(h.base_case(maxstp=0, q_prof_on=0), "eq_qvmec").dump(0)
        np.testing.assert_allclose(d["qq"] * d["qqinv"], 1.0, rtol=0, atol=1e-12)

    @unittest.expectedFailure
    def test_q_and_iota_consistent_with_external_q(self):
        # KNOWN ISSUE: with q_prof_on=1, seteq replaces qq by the q profile
        # from Data.txt (src/equilibrium.f90, "if (q_prof_on .eq. 1)
        # qq=qprofile") but leaves qqinv -- which the operators use -- at the
        # VMEC value; the two disagree by up to ~2% in this case.
        d = h.run_case(h.base_case(maxstp=0, q_prof_on=1), "eq_qext").dump(0)
        np.testing.assert_allclose(d["qq"] * d["qqinv"], 1.0, rtol=0, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
