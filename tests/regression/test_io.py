"""I/O regression tests.

* The two input formats (Input_Model and legacy farin namelists) describe the
  same run and must produce byte-identical binary dumps.
* The fs##### dump written by wrdump round-trips: its header, grid and mode
  list match the input deck, and periodic dumps appear on schedule.
* A continuation run (nstres=1) resumed from a dump must reproduce an
  uninterrupted run. This currently FAILS: after resume the linear operator
  is unstable (fields grow ~10x per step even when resuming from the step-0
  dump), so the check is marked as an expected failure until that is fixed.
"""

import unittest

import numpy as np

import far3d_harness as h


class InputFormats(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        case = h.base_case(maxstp=20)
        cls.im = h.run_case(case, "io_input_model", input_format="input_model")
        cls.fa = h.run_case(case, "io_farin", input_format="farin")

    def test_final_dumps_byte_identical(self):
        a = self.im.path("fs0000z").read_bytes()
        b = self.fa.path("fs0000z").read_bytes()
        self.assertEqual(len(a), len(b))
        self.assertTrue(a == b, "Input_Model and farin runs wrote different final dumps")

    def test_initial_dumps_byte_identical(self):
        self.assertEqual(self.im.path("fs0000_00000").read_bytes(), self.fa.path("fs0000_00000").read_bytes())

    def test_farprt_echoes_input(self):
        # farprt starts with a verbatim copy of the input file
        text = self.fa.farprt.read_text()
        self.assertIn("copy of in:", text)
        self.assertIn("nam_par", text)


class DumpContents(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.case = h.base_case(maxstp=30, ndump=10)
        cls.sim = h.run_case(cls.case, "io_dumps")
        cls.final = cls.sim.final_dump()

    def test_header_matches_input(self):
        d = self.final
        self.assertEqual(d["maxstp"], 30)
        self.assertEqual(d["nstep"], 30)
        self.assertEqual(d["nonlin"], 0)
        self.assertEqual(d["dt0"], 2.0)
        self.assertEqual(d["ext_prof_name"], "Data.txt")
        self.assertEqual((d["mj"], d["lmax"], d["leqmax"]), (200, 29, 21))

    def test_mode_lists_match_input(self):
        mm = [int(x) for x in self.case["mm"].split(",")]
        mmeq = [int(x) for x in self.case["mmeq"].split(",")]
        np.testing.assert_array_equal(self.final["mm"], mm)
        np.testing.assert_array_equal(self.final["nn"], [1] * 4 + [-1] * 4 + [0] * 21)
        np.testing.assert_array_equal(self.final["mmeq"], mmeq)

    def test_grid_is_uniform(self):
        # Auto_grid_on = 1 in the base case
        np.testing.assert_allclose(self.final["r"], np.arange(201) / 200.0, rtol=0, atol=1e-15)

    def test_time_bookkeeping(self):
        self.assertEqual(self.final["dt"], 2.0)
        self.assertAlmostEqual(self.final["time"], 60.0, places=12)

    def test_periodic_dumps(self):
        # ndump=10: dumps at steps 10 and 20 (none at the last step, which
        # writes fs####z instead), numbered from 1; dump 0 is the initial state.
        steps = [self.sim.dump(i)["nstep"] for i in range(3)]
        self.assertEqual(steps, [0, 10, 20])
        self.assertFalse(self.sim.path("fs0000_00003").exists())

    def test_intermediate_dump_matches_shorter_run(self):
        # the state dumped at step 20 equals the final state of a 20-step run
        short = h.run_case(h.base_case(maxstp=20), "io_dumps_short").final_dump()
        mid = self.sim.dump(2)
        for k in h.FIELD_NAMES:
            self.assertEqual(h.rel_diff(mid[k], short[k]), 0.0, k)


class Restart(unittest.TestCase):
    """50 steps + a 50-step continuation must match 100 uninterrupted steps."""

    @classmethod
    def setUpClass(cls):
        cls.full = h.run_case(h.base_case(maxstp=100), "io_restart_full").final_dump()
        h.run_case(h.base_case(maxstp=50), "io_restart_split")
        # widthi = 0: don't add a new seed perturbation on resume (pert is
        # called again by resume and would otherwise add one)
        cont = h.base_case(maxstp=50, widthi="29*0.0").set(nstres=1, numrun="00 01", numruno="00 00 z")
        cls.resumed = h.run_case(cont, "io_restart_split", fresh=False).final_dump()

    def test_step_and_time_continue(self):
        self.assertEqual(self.resumed["nstep"], 100)
        self.assertAlmostEqual(self.resumed["time"], self.full["time"], places=10)

    def test_equilibrium_restored(self):
        for k in h.PROFILE_NAMES + ["sqg", "bmod", "grr", "eta"]:
            self.assertEqual(h.rel_diff(self.resumed[k], self.full[k]), 0.0, k)

    @unittest.expectedFailure
    def test_fields_match_uninterrupted_run(self):
        # KNOWN BUG: the resumed run diverges (~10x per step).
        for k in h.FIELD_NAMES:
            self.assertLess(h.rel_diff(self.resumed[k], self.full[k]), 1e-8, k)


if __name__ == "__main__":
    unittest.main()
