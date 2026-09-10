# DIIID example case

A working linear-stability example for `far3d.x`, verified on macOS
(MacPorts MPICH + gfortran 13). Run from this directory:

```
mpirun -np 1 <path-to>/far3d.x
```

`woutb` was produced from the VMEC equilibrium at
`BOOZ_XFORM/Equilibria/DIIID/wout_d3d_159243_00805.nc` via
`BOOZ_XFORM/Sources/in_booz.d3d_159243_00805`, using this repo's own
`xbooz_xform` (`BOOZ_XFORM/build/xbooz_xform` after `cmake .. && make`
in `BOOZ_XFORM/`):

```
xbooz_xform in_booz.d3d_159243_00805 far
```

`Data.txt` is the experimental profile file (same case as
`FAR3D_PSFC`'s `Tokamaks/DIIID`).

## About `farin`

`src/far3d.f90` reads a fixed-format header line followed by two raw
Fortran namelists (`nam_par`, `nam_arr`) — **not** the `Input_Model`
format described in the bundled `User_guide.pdf` (that guide documents
a newer input scheme used by other FAR3d forks, not what's actually in
`src/`). This `farin` was hand-built by cross-referencing every
`nam_par`/`nam_arr` variable against `src/dfault.f90`'s defaults and
the physical parameters of the equivalent DIIID case.

One non-obvious point if you build another case by hand: `mm`/`mmeq`
must list the n=0 family's poloidal modes **symmetrically** (both +m
and -m), not just m>=0. `initialize.f90`'s `mmlims` sets
`mmstart(-n) = -mmend(n)` for every family including n=0, and
`mult_mod.f90`/`equilibrium.f90` then index `ll(-m,0)`/`lleq(-m,0)` for
the full symmetric range — so the global `mmin` (`minval(mm)`) must
reach at least `-maxval(mm for n=0)`, or you'll hit an out-of-bounds
array index at runtime. That's why this case's `lmax`/`leqmax` (29/21)
are larger than the 19/11 mode counts you'd see in an `Input_Model`
for the same physical case elsewhere.
