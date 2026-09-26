# FAR3d performance benchmarks

A benchmark suite for measuring how FAR3d's run time scales with MPI ranks
and OpenMP threads, and for comparing builds and architectures (CPU vs GPU,
compilers, machines). Every timed run is also checked for correctness, so a
configuration only counts if it gets the right answer.

These are **not** part of `ctest`. `ctest -j` runs tests concurrently, which
makes timings meaningless. Run the benchmarks on an otherwise idle machine.

## Running

```
cd tests/performance
python3 run_benchmarks.py --list                                   # describe the benchmarks
python3 run_benchmarks.py --exe ../../build-release/far3d.x --np 1 2 3 4 --label m2pro
python3 run_benchmarks.py --exe ../../build-release/far3d.x --np 1 --threads 1 2 4 8
python3 run_benchmarks.py --exe ../../build-gpu/far3d.x --np 1 2 4 --arch gpu --label a100
python3 compare_results.py results/m2pro.json results/a100.json
```

Or from `src/`: `make perf PYTHON=/opt/homebrew/bin/python3.11 PERF_ARGS="--np 1 2 4"`
builds the release configuration and benchmarks it.

Useful options:

- `--benchmarks` picks a subset.
- `--repeat` sets timed repetitions (default 3).
- `--mpiexec` and `--mpiexec-args` set the launcher and its options, e.g.
  `--mpiexec-args --bind-to core`. `--mpiexec-args` must come last.
- `--output` sets where results go (default `results/<label>.json`).
- `--workdir` sets where runs happen (default `tests/runs/performance/`).

A full sweep of both benchmarks over 1–4 ranks with `--repeat 2` takes about
4 minutes on an M2 Pro.

## What gets measured

For each (benchmark, ranks, threads) configuration, the driver:

1. Does one untimed warm-up run, so the first launch after a build doesn't
   pay for cold caches. A cold first run was observed to take 13 s instead
   of 1 s.
2. Times the benchmark at a short and a long step count, `--repeat` times
   each, and keeps the fastest of each.
3. Splits the time with `T(N) = T_setup + N * t_step`:
   - `t_step = (T_long - T_short) / (N_long - N_short)` is the time per step
     (`solve`).
   - `T_setup` is everything else: MPI start-up, reading and mapping the
     equilibrium, building and LU-factoring the operator (`linstart`), and
     end-of-run output.

   This avoids instrumenting the Fortran.
4. Validates the long run (below).
5. Reports speedup and parallel efficiency of `t_step` relative to the
   smallest `ranks × threads` in the sweep.

The results JSON also records the machine (CPU, GPU via `nvidia-smi`, core
count), the build (compiler and flags from the `CMakeCache.txt` next to
`far3d.x`), the MPI launcher version, and the git commit. It flags uncommitted
changes to `src/` as `+dirty-src`. Run-to-run noise on a laptop is about
5–10% per step, so use `--repeat 3` or more before trusting small
differences.

## The benchmarks

### `linear_diiid`: linear solve with a known answer

This is the `Models/DIIID/farin` case: linear, mj = 1000, dt = 2, 1000 steps,
no diffusivities. `Models/DIIID/farprt0000` gives, for every variable and all
four n = 1 harmonics:

    gamma = 3.89148e-02,   omega = 3.31558e-02

The original deck has a single toroidal family, and the linear solve is
parallelized over families, so extra MPI ranks would sit idle. The benchmark
therefore adds two independent families, n = 2 (m = 6–9) and n = 3
(m = 10–13), each seeded. In an axisymmetric equilibrium different n don't
couple linearly, so the n = 1 answer is unchanged (verified: identical to 6
digits), and the solve work can use 3 ranks.

Validation checks all 24 n = 1 growth-rate lines against the known answer to
1e-5 relative, about ±4 in the last printed digit. The run converges onto the
eigenmode between 300 and 500 steps. At 100–300 steps the validator rejects
it.

### `nonlinear_diiid`: basic nonlinear solve

This is the `Models/DIIID/Input_Model` deck as is (nonlinear, mj = 400,
dt = 0.01; 8 n = 1 modes plus 21 n = 0 modes), timed over 100 steps, with two
changes:

- `m0dy = 21`, so the n = 0 modes actually evolve. With the deck's
  `m0dy = 0`, `setmod` drops them from the evolved set, nothing is driven
  nonlinearly, and the run is bitwise identical to a linear one.
- A seed amplitude of 1e-4 instead of 1e-140, so the nonlinear terms are not
  negligible. The convolution cost is the same either way, but this makes
  the validation meaningful.

Validation requires:

- all fields finite;
- the n = 0 part of φ actually driven;
- a per-mode fingerprint (max |field| of each mode and field) matching
  `reference/nonlinear_diiid.json` to 1e-6 relative.

Different rank counts agree to about 7e-9 (see the known np-dependence in
`tests/README.md`), well inside that tolerance.

To regenerate the reference after an intentional physics change, run once on
a trusted build:

```
python3 run_benchmarks.py --exe ... --benchmarks nonlinear_diiid --np 1 --repeat 1 --update-reference
```

## Baseline results (Apple M2 Pro, gfortran 13, MPICH, commit 7a5f3ad)

```
benchmark        ranks  thr  setup [s]   step [ms] total [s]  speedup  effic.
------------------------------------------------------------------------------
linear_diiid         1    1       1.22       21.04     22.26     1.00     100%
linear_diiid         2    1       1.00       14.87     15.87     1.41      71%
linear_diiid         3    1       0.79        8.69      9.47     2.42      81%
linear_diiid         4    1       0.84        8.13      8.98     2.59      65%
nonlinear_diiid      1    1       1.80       68.74      8.81     1.00     100%
nonlinear_diiid      2    1       1.73       58.00      7.64     1.19      59%
nonlinear_diiid      3    1       1.78       57.38      7.63     1.20      40%
nonlinear_diiid      4    1       1.80       57.01      7.61     1.21      30%
nonlinear_diiid      1    4       1.10      109.67     12.29     0.68      17%
```

What these show:

- **Linear scales with the number of toroidal families, not ranks.** With 3
  families on 2 ranks the split is 2 + 1, so the speedup is 1.4× rather than
  2×. At 3 ranks each rank has one family (2.4×). A 4th rank has no family to
  solve.
- **Nonlinear barely scales (≈1.2×).** Its mode set has two families, and the
  n = 0 family (21 modes, 147×147 blocks) is solved on one rank. That serial
  solve bounds the speedup (Amdahl), however well the radially decomposed
  nonlinear products parallelize.
- **OpenMP threads make it slower.** The `!$OMP PARALLEL DO` loops in the
  block solver and the convolutions are too fine-grained for this problem
  size: 4 threads cost 1.6× more per step than 1.

## Adding benchmarks

Benchmarks are defined in `benchmarks.py`. Each one is a `Benchmark` with:

- a function that builds the input deck for a given step count;
- short and long step counts;
- a `validate(run)` function that returns a list of problems.

Add it to `BENCHMARKS`. The planned nonlinear mode-count sweep fits this
directly: one function that builds a `Benchmark` from a mode list (e.g.
n = 0..N families with M harmonics each), called for increasing N, with each
case's reference generated via its own `update_reference`. `_mode_lists`
already writes mode lists in FAR3d's (m,n)/(−m,−n) convention.

Reference fingerprints are specific to a deck. If you change a benchmark's
input, regenerate its reference.
