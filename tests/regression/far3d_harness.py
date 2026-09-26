"""Helpers shared by the FAR3d regression tests.

Everything here is driven by environment variables that CTest sets (see
tests/CMakeLists.txt); they can also be set by hand to run a test module
directly:

    FAR3D_EXE       path to far3d.x                        (required)
    FAR3D_MPIEXEC   MPI launcher                           (default: mpirun)
    FAR3D_REPO      repository root                        (default: ../.. from here)
    FAR3D_WORKDIR   scratch directory for the runs         (default: <repo>/tests/runs)
    FAR3D_UPDATE_REFERENCE=1   rewrite reference/*.json instead of comparing

The pieces:
  * Case      -- a FAR3d input deck, parsed from an Input_Model file, that can
                 be modified by variable name and written out in either input
                 format (Input_Model or legacy farin namelists).
  * run_case  -- run far3d.x on a Case in a fresh directory.
  * read_dump -- read a binary fs##### dump (records written by wrdump in
                 src/output_mod.f90) into numpy arrays.
  * read_growth_rates, read_table -- parse farprt and the per-mode energy files.
  * Reference -- load/compare/update JSON reference values.
"""

import json
import os
import re
import shutil
import struct
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = Path(os.environ.get("FAR3D_REPO", HERE.parent.parent)).resolve()
WORKDIR = Path(os.environ.get("FAR3D_WORKDIR", REPO / "tests" / "runs")).resolve()
REFERENCE_DIR = HERE / "reference"
MODEL_DIR = REPO / "Models" / "DIIID"
UPDATE_REFERENCE = os.environ.get("FAR3D_UPDATE_REFERENCE", "") not in ("", "0")


def far3d_exe():
    exe = os.environ.get("FAR3D_EXE")
    if not exe:
        raise RuntimeError("FAR3D_EXE is not set; run the tests through ctest or set it by hand")
    return exe


# ---------------------------------------------------------------------------
# Numbers as Fortran prints them
# ---------------------------------------------------------------------------

_FORTRAN_EXP = re.compile(r"^([+-]?\d*\.?\d*)([+-]\d{3})$")


def fortran_float(text):
    """Parse a Fortran-formatted real, including the E-less 3-digit exponent
    form that Ew.d edit descriptors produce for |exponent| > 99
    (e.g. '-1.37143175-143')."""
    text = text.strip().replace("D", "E").replace("d", "e")
    m = _FORTRAN_EXP.match(text)
    if m:
        return float(m.group(1) + "E" + m.group(2))
    return float(text)


# ---------------------------------------------------------------------------
# Input decks
# ---------------------------------------------------------------------------

def _namelist_members(name):
    """Variable names in namelist `name` as declared in src/far3d.f90, so the
    farin writer always matches what the code actually reads."""
    lines = (REPO / "src" / "far3d.f90").read_text().splitlines()
    for i, line in enumerate(lines):
        if re.search(r"namelist\s*/\s*%s\s*/" % name, line, re.IGNORECASE):
            text = line.split("/", 2)[2]
            while text.rstrip().endswith("&"):
                i += 1
                text = text.rstrip()[:-1] + lines[i]
            return [v.strip().lower() for v in text.split(",") if v.strip()]
    raise RuntimeError("namelist %s not found in src/far3d.f90" % name)


class Case:
    """A FAR3d input deck.

    Parsed from an Input_Model file ("!!!!!!!!!!! name: description" comment
    line followed by a value line). Values are kept as the raw Fortran text so
    that rewriting a deck does not perturb any number."""

    HEADER = ("nstres", "numrun", "numruno", "eq_name")

    def __init__(self, entries):
        self.entries = entries          # list of [name, comment, value]

    @classmethod
    def from_input_model(cls, path=MODEL_DIR / "Input_Model"):
        lines = Path(path).read_text().splitlines()
        entries = []
        for i in range(0, len(lines) - 1, 2):
            comment, value = lines[i], lines[i + 1]
            m = re.match(r"^!+\s*([A-Za-z0-9_]+)", comment)
            if not m:
                raise ValueError("unexpected Input_Model line %d: %r" % (i + 1, comment))
            entries.append([m.group(1).lower(), comment, value.strip()])
        return cls(entries)

    def copy(self):
        return Case([list(e) for e in self.entries])

    def _entry(self, name):
        for e in self.entries:
            if e[0] == name.lower():
                return e
        raise KeyError(name)

    def __getitem__(self, name):
        return self._entry(name)[2]

    def set(self, **values):
        """Set values by name. Python values are converted to Fortran text:
        bool -> .true./.false., float -> repr, list -> comma separated."""
        for name, v in values.items():
            self._entry(name)[2] = _to_fortran(v)
        return self

    # -- writers ------------------------------------------------------------

    def write_input_model(self, path):
        with open(path, "w") as f:
            for name, comment, value in self.entries:
                f.write(comment + "\n" + value + "\n")

    def write_farin(self, path):
        nam_par = _namelist_members("nam_par")
        nam_arr = _namelist_members("nam_arr")
        numrun = self["numrun"].split()
        numruno = self["numruno"].split()
        # read by far3d.f90 with format (i1,2x,2a2,2x,2a2,a1,2x,a40)
        header = "%1d  %2s%2s  %2s%2s%1s  %s" % (int(self["nstres"]), numrun[0], numrun[1],
                                                 numruno[0], numruno[1], numruno[2], self["eq_name"])
        with open(path, "w") as f:
            f.write(header + "\n")
            for group, members in (("nam_par", nam_par), ("nam_arr", nam_arr)):
                f.write("&%s\n" % group)
                for name in members:
                    value = self[name]
                    if not _is_fortran_literal(value):
                        value = "'%s'" % value
                    f.write(" %s=%s,\n" % (name, value))
                f.write("/\n")


def _to_fortran(v):
    if isinstance(v, bool):
        return ".true." if v else ".false."
    if isinstance(v, (list, tuple)):
        return ",".join(_to_fortran(x) for x in v)
    if isinstance(v, float):
        return repr(v)
    return str(v)


def _is_fortran_literal(value):
    """True for numbers, logicals and comma lists of them (incl. r*v repeats)."""
    for item in value.split(","):
        item = item.strip()
        if "*" in item:
            item = item.split("*", 1)[1]
        if item.lower() in (".true.", ".false.", "t", "f"):
            continue
        try:
            fortran_float(item)
        except ValueError:
            return False
    return True


def base_case(**overrides):
    """The small DIII-D case used by the regression tests: the Models/DIIID
    Input_Model deck, cut down to 200 radial points so a run takes seconds,
    as a linear run with the time step of Models/DIIID/farin."""
    case = Case.from_input_model()
    case.set(mj=200, edge_p=190, dt0=2.0, nonlin=0, maxstp=200,
             ndump=100000, nprint=50, ndiag=100000)
    case.set(**overrides)
    return case


# ---------------------------------------------------------------------------
# Running far3d.x
# ---------------------------------------------------------------------------

class Run:
    def __init__(self, directory, numrun, stdout):
        self.dir = Path(directory)
        self.numrun = numrun            # e.g. "0000"
        self.stdout = stdout

    def path(self, name):
        return self.dir / name

    @property
    def farprt(self):
        return self.dir / ("farprt" + self.numrun)

    def final_dump(self):
        return read_dump(self.dir / ("fs" + self.numrun + "z"))

    def dump(self, index):
        return read_dump(self.dir / ("fs%s_%05d" % (self.numrun, index)))


def run_case(case, name, nproc=1, input_format="input_model", fresh=True, timeout=600):
    """Run far3d.x on `case` in WORKDIR/name and return a Run.

    fresh=False reuses an existing directory (for continuation runs, which
    read the previous run's dump from the working directory)."""
    rundir = WORKDIR / name
    if fresh and rundir.exists():
        shutil.rmtree(rundir)
    rundir.mkdir(parents=True, exist_ok=True)
    for f in ("woutb", "Data.txt"):
        if not (rundir / f).exists():
            shutil.copy(MODEL_DIR / f, rundir / f)
    for f in ("Input_Model", "farin"):
        if (rundir / f).exists():
            (rundir / f).unlink()

    if input_format == "input_model":
        case.write_input_model(rundir / "Input_Model")
    elif input_format == "farin":
        case.write_farin(rundir / "farin")
    else:
        raise ValueError(input_format)

    mpiexec = os.environ.get("FAR3D_MPIEXEC") or "mpirun"
    cmd = [mpiexec, "-n", str(nproc), far3d_exe()]
    proc = subprocess.run(cmd, cwd=rundir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          text=True, timeout=timeout)
    (rundir / "stdout.log").write_text(proc.stdout)
    if proc.returncode != 0 or "Simulation DONE" not in proc.stdout:
        raise RuntimeError("far3d.x failed in %s (exit %d):\n%s" % (rundir, proc.returncode, proc.stdout[-3000:]))
    return Run(rundir, "".join(case["numrun"].split()), proc.stdout)


# ---------------------------------------------------------------------------
# Binary dumps (fs#####), written by wrdump in src/output_mod.f90
# ---------------------------------------------------------------------------

def _records(path):
    """Fortran sequential unformatted, big-endian, 4-byte record markers."""
    data = Path(path).read_bytes()
    pos, out = 0, []
    while pos < len(data):
        (n,) = struct.unpack(">i", data[pos:pos + 4])
        out.append(data[pos + 4:pos + 4 + n])
        (n2,) = struct.unpack(">i", data[pos + 4 + n:pos + 8 + n])
        if n2 != n:
            raise ValueError("corrupt record markers in %s" % path)
        pos += n + 8
    return out


class _Cursor:
    def __init__(self, buf):
        self.buf, self.pos = buf, 0

    def take(self, fmt, count=1):
        size = struct.calcsize(">" + fmt) * count
        vals = struct.unpack(">" + fmt * count, self.buf[self.pos:self.pos + size])
        self.pos += size
        return vals if count > 1 else vals[0]

    def ints(self, n):
        v = np.frombuffer(self.buf, ">i4", n, self.pos).astype(np.int64)
        self.pos += 4 * n
        return v

    def reals(self, n):
        v = np.frombuffer(self.buf, ">f8", n, self.pos).astype(np.float64)
        self.pos += 8 * n
        return v

    def chars(self, n):
        v = self.buf[self.pos:self.pos + n].decode("ascii", "replace")
        self.pos += n
        return v

    def done(self):
        return self.pos == len(self.buf)


PROFILE_NAMES = ["qq", "qqinv", "preq", "feq", "cureq", "denseq", "teeq", "nfeq", "vfova", "vzt_eq",
                 "tieq", "vth_eq"]
METRIC_NAMES = ["sqgi", "sqg", "bst", "grr", "grt", "gtt", "grroj", "grtoj", "gttoj", "jbgrr", "jbgrt",
                "jbgtt", "lplrr", "lplrt", "lplrz", "lpltt", "lpltz", "lplzz", "lplr", "lplt", "lplz",
                "djroj", "djtoj", "djzoj", "omdr", "omdt", "omdz", "bmod", "dbsjtoj", "dbsjzoj"]
FIELD_NAMES = ["psi", "phi", "pr", "nf", "vprlf", "vthprlf"]


def read_dump(path):
    """Read an fs##### dump. Returns a dict of header scalars and numpy arrays:
    radial profiles have shape (mj+1,), equilibrium Fourier coefficients
    (leqmax, mj+1), and perturbed fields (lmax, mj+1), indexed [l-1, j]."""
    recs = _records(path)
    d = {}

    c = _Cursor(recs[0])
    ihist = c.take("i")
    d["numrun"], d["numruno"], d["numruns"] = c.chars(6), c.chars(6), c.chars(6)
    c.chars(8 * ihist)
    for k in ("maxstp", "nstep", "ndump", "nprint", "lplots", "itime"):
        d[k] = c.take("i")
    d["dt0"] = c.take("d")
    for k in ("nstep1", "nonlin", "m0dy", "nocpl", "noprevol_on"):
        d[k] = c.take("i")
    for k in ("adens", "bdens", "lca0", "lca1", "lca2", "lca3", "omcy", "dpres"):
        d[k] = c.take("d")
    d["ext_prof"], d["epflr_on"] = c.take("i"), c.take("i")
    d["r_epflr"] = c.take("d")
    d["alpha_on"], d["iflr_on"] = c.take("i"), c.take("i")
    d["iflr"] = c.take("d")
    d["twofl_on"], d["ieldamp_on"] = c.take("i"), c.take("i")
    c.take("d", 10)
    d["ext_prof_name"] = c.chars(40).strip()
    d["trapped_on"] = c.take("i")
    c.take("d", 2)
    d["b_par_on"] = c.take("i")
    assert c.done(), "unexpected header record length in %s" % path

    c = _Cursor(recs[1])
    mj, lmax, leqmax, _ = c.ints(4)
    d.update(mj=int(mj), lmax=int(lmax), leqmax=int(leqmax))
    nj = mj + 1

    c = _Cursor(recs[2])
    d["r"] = c.reals(nj)
    d["mm"], d["nn"] = c.ints(lmax), c.ints(lmax)
    d["mmeq"], d["nneq"] = c.ints(leqmax), c.ints(leqmax)
    d["rinv"] = c.reals(nj)
    for k in ("dc1m", "dc1p", "dc2m", "dc2p", "del2cm", "del2cp"):
        d[k] = c.reals(mj)
    assert c.done()

    i = 4
    c = _Cursor(recs[i]); i += 1
    for k in PROFILE_NAMES:
        d[k] = c.reals(nj)
    if d["alpha_on"] == 1:
        c = _Cursor(recs[i]); i += 1
        d["nalpeq"], d["valphaova"] = c.reals(nj), c.reals(nj)
    c = _Cursor(recs[i]); i += 1
    if d["ext_prof"] == 1:
        d["valfven"], d["vtherm_elecp"] = c.reals(nj), c.reals(nj)
    else:
        d["vtherm_elc"] = c.reals(nj)
    c = _Cursor(recs[i]); i += 1
    d["ndevice"] = c.chars(16)
    d["eps"], d["bet0"], d["bet0_f"] = c.take("d", 3)
    d["rs"] = c.reals(lmax)
    c = _Cursor(recs[i]); i += 1
    for k in METRIC_NAMES:
        d[k] = c.reals(nj * leqmax).reshape(leqmax, nj)
    assert c.done()
    if d["ieldamp_on"] == 1:
        i += 1
    if d["trapped_on"] == 1:
        i += 1
    c = _Cursor(recs[i]); i += 1
    d["eta"] = c.reals(nj)
    d["dt"], d["time"] = c.take("d", 2)
    names = FIELD_NAMES + (["nalp", "vprlalp"] if d["alpha_on"] == 1 else [])
    for k in names:
        c = _Cursor(recs[i]); i += 1
        d[k] = c.reals(nj * lmax).reshape(lmax, nj)
        assert c.done()
    return d


def mode_index(dump, m, n):
    """0-based row of mode (m, n) in the field arrays of a dump."""
    hits = np.nonzero((dump["mm"] == m) & (dump["nn"] == n))[0]
    if len(hits) != 1:
        raise KeyError((m, n))
    return int(hits[0])


# ---------------------------------------------------------------------------
# Text outputs
# ---------------------------------------------------------------------------

_GROWTH = re.compile(r"^\s*(\w+)\s*:\s*m=\s*(-?\d+)\s+n=\s*(-?\d+)\s+gam=\s*(\S+)\s+om_r=\s*(\S+)")


def read_growth_rates(farprt):
    """Growth rate / frequency lines printed by lincheck, as a list of
    dicts {var, m, n, gamma, omega}."""
    out = []
    for line in Path(farprt).read_text().splitlines():
        m = _GROWTH.match(line)
        if m:
            out.append(dict(var=m.group(1), m=int(m.group(2)), n=int(m.group(3)),
                            gamma=fortran_float(m.group(4)), omega=fortran_float(m.group(5))))
    return out


def read_table(path):
    """Tab-separated time history files (eke_####, eme_####, hifrq_####, ...):
    returns (column names, 2-D array)."""
    lines = Path(path).read_text().splitlines()
    names = [s.strip() for s in lines[0].split("\t")]
    rows = [[fortran_float(x) for x in line.split()] for line in lines[1:] if line.strip()]
    return names, np.array(rows)


# ---------------------------------------------------------------------------
# Reference values
# ---------------------------------------------------------------------------

def rel_diff(a, b):
    """max|a-b| / max|b| over arrays (0 if both are zero)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    scale = np.max(np.abs(b)) if b.size else 0.0
    err = np.max(np.abs(a - b)) if a.size else 0.0
    return 0.0 if scale == 0.0 and err == 0.0 else err / (scale if scale else 1.0)


class Reference:
    """JSON file of reference values for one test module.

    check(testcase, key, value, rtol) compares `value` (number or nested lists)
    with the stored one; with FAR3D_UPDATE_REFERENCE=1 it records `value`
    instead. Call save() from tearDownClass."""

    def __init__(self, name):
        self.path = REFERENCE_DIR / (name + ".json")
        self.data = json.loads(self.path.read_text()) if self.path.exists() else {}
        self.dirty = False

    def check(self, testcase, key, value, rtol, atol=0.0):
        value = _jsonable(value)
        if UPDATE_REFERENCE:
            self.data[key] = value
            self.dirty = True
            return
        if key not in self.data:
            testcase.fail("no reference value for %r in %s; run with FAR3D_UPDATE_REFERENCE=1" %
                          (key, self.path.name))
        ref = np.asarray(self.data[key], float)
        val = np.asarray(value, float)
        testcase.assertEqual(ref.shape, val.shape, "shape of %r changed" % key)
        err = np.max(np.abs(val - ref)) if ref.size else 0.0
        scale = np.max(np.abs(ref)) if ref.size else 0.0
        testcase.assertLessEqual(err, atol + rtol * scale,
                                 "%s: max abs diff %.3e exceeds %.1e + %.1e*%.3e" % (key, err, atol, rtol, scale))

    def save(self):
        if self.dirty:
            self.path.write_text(json.dumps(self.data, indent=1, sort_keys=True) + "\n")


def _jsonable(v):
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    return v
