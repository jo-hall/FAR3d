FAR3d version 2.0

(This is a parallel version that requires MPI)

The earlier version of FAR3d has now been moved to the PreviousFar3d branch

Both CPU and GPU versions are included in this code. 

(to be filled in)

## Timing output

Set `timing_on=1` to have FAR3d write per-time-step wall-clock timings to
`timing_<numrun>` (e.g. `timing_0000`). It is off by default.

- `farin`: add `timing_on=1,` to the `&nam_par` namelist.
- `Input_Model`: append these two lines at the end of the file (the field is
  optional, so existing Input_Model files still work):

  ```
  !!!!!!!!!!! timing_on: 1 = write per-step timing file timing_<numrun>
  1
  ```

The file has a header (ranks, OpenMP threads, cpu/gpu build, run size, and
the setup time split into initialization and `linstart`), then one line per
time step, then a summary line with the per-column totals and the time after
the loop (`finalize_s`). Lines starting with `#` are comments, so it loads
directly with e.g. `numpy.loadtxt`. Columns, in seconds:

| column | what it covers |
|---|---|
| `step`, `time` | step number and simulation time |
| `total` | wall time of the whole step |
| `linear` | r.h.s. build and block-tridiagonal solves (`solbt`) |
| `nonlinear` | explicit nonlinear terms (zero for linear runs) |
| `comm` | MPI layout transposes (`trnsfr`) and radial halo exchanges, including time spent waiting for other ranks there |
| `gather` | gathers to rank 0 (`trnsfr0`, `trnsfr0e`) |
| `diag` | diagnostics and output in the loop (`energy`, `hifreq`, `wrdump`) |
| `other` | the rest of the step |

Timers nest and are exclusive (time is charged to the innermost category),
so on one rank the categories add up to `total`. With several ranks, each
column is the maximum over ranks taken separately, so `total` need not equal
their sum. The timers are in `src/timers.f90`.
