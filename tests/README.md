# FAR3d tests

Two layers, both run through CTest:

| label | what | runtime |
|---|---|---|
| `unit` | Fortran programs that call individual FAR3d routines and check them against exact answers | < 1 s |
| `regression` | Python scripts that run `far3d.x` on a small DIII-D case and check physics properties and stored reference values | ~25 s with `-j4` |

## Running

From `src/`, the usual build directory:

```
make test PYTHON=/opt/homebrew/bin/python3.11     # build release + run everything
make test CTEST_ARGS="-L unit"                     # just the unit tests
```

Or from any CMake build directory:

```
ctest --output-on-failure -j4
ctest -L regression -V          # verbose: one line per Python test case
ctest -R test_linear            # one module
```

The regression tests need a Python 3 with numpy. CMake uses the first
`python3` it finds and disables the regression tests, with a warning, if that
one has no numpy. Point it at the right interpreter with
`-DFAR3D_TEST_PYTHON=/path/to/python3` (or `PYTHON=` with `make`).

Configure with `-DBUILD_TESTING=OFF` to skip building the tests entirely.

Run output goes to `<build>/tests/runs/<test-name>/`, so a failing case can
be inspected (inputs, `farprt`, dumps, `stdout.log`) after the fact.

## Layout

```
tests/
  CMakeLists.txt          test registration
  unit/
    testing_mod.f90       check_* assertions, PASS/FAIL/XFAIL output, exit code
    test_block_tridiag.f90
    test_radial_grid.f90
    test_spectral_ops.f90
    test_numerics_tools.f90
  regression/
    far3d_harness.py      input decks, running far3d.x, reading dumps/farprt, reference values
    test_io.py
    test_equilibrium.py
    test_linear.py
    test_nonlinear.py
    reference/*.json      stored reference values (golden master)
```

The unit tests link against `far3d_core`, a static library of every source
file except `far3d.f90`, which the top-level `CMakeLists.txt` builds for both
`far3d.x` and the tests.

## What is tested

### Unit (Fortran)

- **`test_block_tridiag`**: `decbt`/`solbt` (the time-step linear solver)
  recover a known solution of random block-tridiagonal systems, including the
  (1,3)/(n,n−2) corner blocks, for block sizes up to 56 (the DIII-D case).
- **`test_radial_grid`**: the uniform and packed (`findr`) grids have the
  documented shape. The nonuniform finite-difference weights are exact on
  quadratics and converge at second order, and the cylindrical Laplacian
  `del2c` is exact on rᵐ.
- **`test_spectral_ops`**: `mult` (the Fourier-space product behind every
  nonlinear term), `dbydth_par` and `grdpar` match brute-force evaluation on a
  (θ,ζ) grid, for all four cos/sin parity combinations. This pins down the
  storage convention: for the positive member (m,n) of each pair, a type +1
  field holds cos at `ll(m,n)` and sin at `ll(-m,-n)`, and a type −1 field
  holds the reverse.
- **`test_numerics_tools`**: the equilibrium spline reproduces cubics and
  converges at 4th order. The plasma dispersion function `zzdisp` is checked
  against exact identities (Z(0), real and imaginary axes, asymptotics), and
  the `tools` `erf` is checked against the intrinsic.

### Regression (Python + far3d.x)

The base case is `Models/DIIID/Input_Model` at mj = 200, run as a linear case
with dt = 2 (`far3d_harness.base_case`).

- **`test_io`**
  - `Input_Model` and `farin` decks produce byte-identical dumps.
  - Dump header, mode lists and grid match the input.
  - Periodic dumps appear on schedule and equal shorter runs.
  - Continuation runs restore step count, time and equilibrium.
- **`test_equilibrium`**
  - BOOZ_XFORM regenerates `Models/DIIID/woutb` from the VMEC `.nc`. This is
    skipped unless `BOOZ_XFORM/build/xbooz_xform` exists or `FAR3D_XBOOZ` is
    set.
  - Mapped profiles, metric coefficients and resonant radii match the
    reference.
  - Axis normalization, positivity, and up-down symmetry hold.
  - mj = 400 reproduces mj = 200 at shared points.
- **`test_linear`**
  - lincheck's γ and ω describe a single eigenmode and match the reference.
  - The solution is exactly linear in the seed amplitude.
  - The fields converge at 2nd order in dt.
  - Results are bitwise identical for 1, 2 and 3 MPI processes.
- **`test_nonlinear`** (uses `m0dy = 21`; with the example's `m0dy = 0` no
  n = 0 mode evolves and a nonlinear run is identical to a linear one)
  - The tiny-amplitude limit reproduces the linear run.
  - Doubling the seed doubles n = 1 and quadruples n = 0 (quadratic coupling).
  - A strongly nonlinear run matches the reference fields and energy
    histories.

## Known issues encoded as expected failures

These checks assert the *correct* behavior and are marked
`@unittest.expectedFailure` (Python) or `check_known_bug` (Fortran). They
show as "expected failure" / `XFAIL` and don't fail the suite. When the
underlying bug is fixed they turn into "unexpected success" / `XPASS`, which
is the signal to make them regular checks.

| check | observation |
|---|---|
| `test_io.Restart.test_fields_match_uninterrupted_run` | A continuation run (`nstres=1`) is numerically unstable: fields grow ~10× per step, even when resuming from the step-0 dump and with no new perturbation (`widthi=0`). Equilibrium, η and dt are restored bitwise, so the problem is in how `resume` rebuilds the operator. |
| `test_nonlinear.MpiInvariance.test_process_count_does_not_change_result` | Nonlinear runs differ between np=1 and np>1 (~1e-8 relative after one step, ~3e-6 after 20). np=2 and np=3 agree with each other exactly, which points to a single- vs multi-process code-path difference rather than to where the grid is split. Linear runs are bitwise invariant. |
| `test_equilibrium.SafetyFactor.test_q_and_iota_consistent_with_external_q` | With `q_prof_on=1`, `qq` is replaced by the `Data.txt` profile but `qqinv` (ι, which the operators use) keeps the VMEC value; they differ by up to ~2%. |
| `test_numerics_tools` (3 × XFAIL) | `zzdisp` returns Re Z with the wrong sign for a purely real negative argument (it only reflects when x·y < 0). In FAR3d this needs `omegar < 0` with zero collisionality. |

## Updating reference values

The JSON files in `regression/reference/` were generated by this code on
macOS/arm64 with gfortran 13 and MPICH. After an intentional change to the
results, regenerate them and review the diff:

```
cd tests/regression
FAR3D_UPDATE_REFERENCE=1 FAR3D_EXE=<build>/far3d.x python3 -m unittest test_equilibrium test_linear test_nonlinear
git diff reference/
```

Tolerances are set for the same machine and compiler (roughly 1e-7 to 1e-9
on dumped values). lincheck output only has 6 significant digits, so γ and ω
use 5e-5. A different compiler or optimization level may need
looser values.

## Running a module by hand

The Python tests read their configuration from the environment that CTest
normally sets:

```
FAR3D_EXE=<build>/far3d.x FAR3D_MPIEXEC=mpiexec FAR3D_WORKDIR=/tmp/far3d-runs \
    python3 -m unittest -v test_linear.LinearBaseCase
```
