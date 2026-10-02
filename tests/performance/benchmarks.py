"""FAR3d performance benchmark definitions.

Each Benchmark is a FAR3d run that is (a) big enough that its cost is
dominated by the solver rather than process startup, and (b) checkable, so
a timing on a new machine, process count or architecture only counts if the
answer is still right.

Timing. Every run is made with timing_on=1, so FAR3d itself writes the
setup time (initialization and linstart) and, for every time step, the step
time split into linear / nonlinear / comm / gather / diag (src/timers.f90,
read with far3d_harness.read_timing). The driver runs each configuration for
`steps` time steps, after an untimed warm-up of `warmup_steps`, and
validates the timed run.

Adding a benchmark: write a function returning a Benchmark and add it to
BENCHMARKS. A nonlinear mode-count sweep, for example, would be a family of
Benchmarks built by one function with the mode list as a parameter.
"""

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List

import numpy as np

import far3d_harness as h

HERE = Path(__file__).resolve().parent
REFERENCE_DIR = HERE / "reference"


@dataclass
class Benchmark:
    name: str
    description: str
    make_case: Callable[[int], "h.Case"]     # steps -> input deck
    warmup_steps: int                         # untimed warm-up run length
    steps: int                                # timed (and validated) run length
    validate: Callable[["h.Run"], List[str]]  # returns a list of problems (empty = valid)
    parallelism: str                          # what the MPI decomposition can use, for the report
    extra: Dict = field(default_factory=dict)


def _mode_lists(families, eq_modes):
    """mm/nn lists in FAR3d's convention: every (m, n) of the dynamic families,
    then their (-m, -n) partners, then the n=0 equilibrium modes (which only
    evolve when m0dy covers them)."""
    mm, nn = [], []
    for n, ms in families:
        mm += ms
        nn += [n] * len(ms)
    for n, ms in families:
        mm += [-m for m in ms]
        nn += [-n] * len(ms)
    mm += list(eq_modes)
    nn += [0] * len(eq_modes)
    return mm, nn


# ---------------------------------------------------------------------------
# Linear: DIII-D Alfven eigenmode with a known answer
# ---------------------------------------------------------------------------

# Models/DIIID/farin (mj=1000, dt=2, 1000 steps, n=1 harmonics m=1..4) gives,
# for every variable and every harmonic (Models/DIIID/farprt0000):
KNOWN_GAMMA = 3.89148e-02
KNOWN_OMEGA = 3.31558e-02
# lincheck prints 6 significant digits; 1e-5 allows about +-4 units in the last
# one for other compilers and architectures, while rejecting unconverged runs.
KNOWN_RTOL = 1e-5

LINEAR_FAMILIES = [
    (1, [4, 3, 2, 1]),       # the Models/DIIID n=1 harmonics: carries the known answer
    (2, [9, 8, 7, 6]),       # extra toroidal families: independent of n=1 in an
    (3, [13, 12, 11, 10]),   # axisymmetric equilibrium, so the n=1 answer is unchanged
]


def linear_case(steps):
    mm, nn = _mode_lists(LINEAR_FAMILIES, range(-10, 11))
    lmax = len(mm)
    # seed the first harmonic of each family
    widthi = ["0.0"] * lmax
    for i in (0, 4, 8):
        widthi[i] = "1.e-140"
    # Models/DIIID/farin physics: linear, dt=2, no diffusivities
    return h.Case.from_input_model().set(
        nonlin=0, dt0=2.0, mj=1000, edge_p=990, maxstp=steps, m0dy=0,
        ndump=100000, nprint=100000, ndiag=100000,
        stdifp=0.0, stdifu=0.0, stdifv=0.0, stdifnf=0.0, stdifvf=0.0, stdifnalp=0.0, stdifvalp=0.0,
        lmax=lmax, mm=mm, nn=nn, widthi=",".join(widthi), gammai="%d*0.0" % lmax)


def linear_validate(run):
    problems = []
    growth = [g for g in h.read_growth_rates(run.farprt) if g["n"] == 1]
    if len(growth) != 24:
        return ["expected 24 n=1 growth-rate lines (4 harmonics x 6 variables), got %d" % len(growth)]
    for g in growth:
        for key, known in (("gamma", KNOWN_GAMMA), ("omega", KNOWN_OMEGA)):
            if abs(g[key] - known) > KNOWN_RTOL * known:
                problems.append("%s m=%d n=1 %s = %.6e, known answer %.6e" % (g["var"], g["m"], key, g[key], known))
    return problems


def linear_diiid():
    return Benchmark(
        name="linear_diiid",
        description="Linear DIII-D TAE, mj=1000, 3 toroidal families (n=1,2,3; 4 harmonics each), "
                    "1000 steps of dt=2. n=1 must reproduce the known gamma/omega of Models/DIIID.",
        make_case=linear_case,
        warmup_steps=1,
        steps=1000,
        validate=linear_validate,
        parallelism="linear solve distributed over 3 toroidal families: MPI speedup saturates at 3 ranks",
    )


# ---------------------------------------------------------------------------
# Nonlinear: the Models/DIIID Input_Model deck
# ---------------------------------------------------------------------------

NONLINEAR_REFERENCE = REFERENCE_DIR / "nonlinear_diiid.json"
NONLINEAR_RTOL = 1e-6


def nonlinear_case(steps):
    # Models/DIIID/Input_Model as is (mj=400, 29 modes, dt=0.01, nonlin=1)
    # except: m0dy=21 so the 21 n=0 modes evolve (with m0dy=0 nothing is
    # nonlinearly driven and the run is effectively linear), and a finite
    # seed so the nonlinear terms are not negligible.
    return h.Case.from_input_model().set(
        maxstp=steps, m0dy=21, widthi="1.e-4,28*0.0",
        ndump=100000, nprint=100000, ndiag=100000)


def nonlinear_signature(dump):
    """Per-mode max |field| of every field: a compact fingerprint of the state."""
    return {k: np.abs(dump[k]).max(axis=1).tolist() for k in h.FIELD_NAMES}


def nonlinear_validate(run):
    dump = run.final_dump()
    problems = []
    for k in h.FIELD_NAMES:
        if not np.all(np.isfinite(dump[k])):
            problems.append("%s contains non-finite values" % k)
    n0 = dump["nn"] == 0
    ratio = np.abs(dump["phi"][n0]).max() / np.abs(dump["phi"][~n0]).max()
    if not ratio > 1e-4:
        problems.append("n=0 part of phi not driven (n0/n1 = %.2e): nonlinear terms inactive?" % ratio)
    if not NONLINEAR_REFERENCE.exists():
        problems.append("no reference %s; generate it with --update-reference" % NONLINEAR_REFERENCE.name)
        return problems
    ref = json.loads(NONLINEAR_REFERENCE.read_text())
    if ref["steps"] != dump["nstep"]:
        problems.append("reference is for %d steps, run has %d" % (ref["steps"], dump["nstep"]))
        return problems
    sig = nonlinear_signature(dump)
    for k in h.FIELD_NAMES:
        err = h.rel_diff(sig[k], ref["signature"][k])
        if err > NONLINEAR_RTOL:
            problems.append("%s differs from reference by %.2e (tolerance %.0e)" % (k, err, NONLINEAR_RTOL))
    return problems


def nonlinear_update_reference(run):
    dump = run.final_dump()
    REFERENCE_DIR.mkdir(exist_ok=True)
    NONLINEAR_REFERENCE.write_text(json.dumps(
        {"steps": int(dump["nstep"]), "signature": nonlinear_signature(dump)}, indent=1) + "\n")


def nonlinear_diiid():
    return Benchmark(
        name="nonlinear_diiid",
        description="Nonlinear DIII-D, Models/DIIID/Input_Model (mj=400, 8 n=1 + 21 n=0 modes, dt=0.01) "
                    "with m0dy=21 and seed 1e-4, 100 steps.",
        make_case=nonlinear_case,
        warmup_steps=2,
        steps=102,
        validate=nonlinear_validate,
        parallelism="linear solve over 2 families (n=0 with 21 modes dominates); "
                    "nonlinear products over radial slabs (all ranks)",
        extra={"update_reference": nonlinear_update_reference},
    )


BENCHMARKS = {b.name: b for b in (linear_diiid(), nonlinear_diiid())}
