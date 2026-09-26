#!/usr/bin/env python3
"""Run the FAR3d performance benchmarks and record the timings.

Examples
--------
  # strong scaling over MPI ranks on this machine
  python3 run_benchmarks.py --exe ../../build-release/far3d.x --np 1 2 3 4 --label m1-cpu

  # OpenMP threads instead
  python3 run_benchmarks.py --exe ../../build-release/far3d.x --np 1 --threads 1 2 4 8

  # a GPU build, one rank per GPU
  python3 run_benchmarks.py --exe ../../build-gpu/far3d.x --np 1 2 4 --label a100 --arch gpu

  # compare result files (e.g. CPU vs GPU, or two machines)
  python3 compare_results.py results/m1-cpu.json results/a100.json

Each (benchmark, ranks, threads) configuration is run at the benchmark's
short and long step counts, --repeat times each after one untimed warm-up,
and the fastest time of each is used (the usual choice for benchmarks: noise
only ever adds time). From the two step counts the driver separates the
fixed setup cost from the cost per time step; see benchmarks.py. Every long
run is validated, and a configuration with a wrong answer is reported as
INVALID rather than timed.

Timings are only meaningful on an otherwise idle machine, which is why these
benchmarks are not part of ctest (ctest -j runs tests concurrently).
"""

import argparse
import datetime
import json
import os
import platform
import re
import socket
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "regression"))

import far3d_harness as h        # noqa: E402
from benchmarks import BENCHMARKS  # noqa: E402


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------

def _run(cmd, **kw):
    try:
        return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                              timeout=30, **kw).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return ""


def cpu_model():
    if sys.platform == "darwin":
        return _run(["sysctl", "-n", "machdep.cpu.brand_string"])
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.lower().startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor()


def gpu_model():
    out = _run(["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"])
    return ", ".join(sorted(set(out.splitlines()))) if out and "not found" not in out.lower() else ""


def build_info(exe):
    """Compiler and flags from the CMakeCache.txt next to far3d.x, if any."""
    info = {}
    cache = Path(exe).resolve().parent / "CMakeCache.txt"
    if cache.exists():
        for line in cache.read_text().splitlines():
            m = re.match(r"^(CMAKE_BUILD_TYPE|CMAKE_Fortran_COMPILER|CMAKE_Fortran_FLAGS|"
                         r"CMAKE_Fortran_FLAGS_RELEASE|CMAKE_Fortran_COMPILER_ID|"
                         r"CMAKE_Fortran_COMPILER_VERSION):[A-Z]+=(.*)$", line)
            if m:
                info[m.group(1)] = m.group(2)
    compiler = info.get("CMAKE_Fortran_COMPILER")
    if compiler:
        info["compiler_version"] = _run([compiler, "--version"]).splitlines()[0:1]
    return info


def metadata(args):
    commit = _run(["git", "-C", str(h.REPO), "rev-parse", "--short", "HEAD"])
    dirty = _run(["git", "-C", str(h.REPO), "status", "--porcelain", "--untracked-files=no", "src"])
    return {
        "label": args.label,
        "arch": args.arch,
        "date": datetime.datetime.now().isoformat(timespec="seconds"),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "cpu": cpu_model(),
        "cpu_count": os.cpu_count(),
        "gpu": gpu_model(),
        "exe": str(Path(args.exe).resolve()),
        "git_commit": commit + ("+dirty-src" if dirty else ""),
        "build": build_info(args.exe),
        "mpiexec": args.mpiexec,
        "mpiexec_version": _run([args.mpiexec, "--version"]).splitlines()[0:2],
        "mpiexec_args": args.mpiexec_args,
        "repeat": args.repeat,
    }


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------

def time_config(bench, nproc, threads, args, log):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    workdir = Path(args.workdir) / bench.name
    tag = "np%d_t%d" % (nproc, threads)

    def launch(steps, name):
        return h.run_case(bench.make_case(steps), "%s_%s" % (tag, name), nproc=nproc, workdir=workdir,
                          exe=args.exe, mpiexec_args=args.mpiexec_args, env=env, timeout=args.timeout)

    launch(bench.short_steps, "warmup")
    t_short, t_long, long_run = [], [], None
    for i in range(args.repeat):
        t_short.append(launch(bench.short_steps, "short").elapsed)
        long_run = launch(bench.long_steps, "long")
        t_long.append(long_run.elapsed)
        log("    repeat %d: %d steps %.2fs, %d steps %.2fs" %
            (i + 1, bench.short_steps, t_short[-1], bench.long_steps, t_long[-1]))

    if args.update_reference and "update_reference" in bench.extra:
        bench.extra["update_reference"](long_run)
        log("    reference updated from this run")
    problems = bench.validate(long_run)

    ts, tl = min(t_short), min(t_long)
    t_step = (tl - ts) / (bench.long_steps - bench.short_steps)
    return {
        "benchmark": bench.name,
        "nproc": nproc,
        "threads": threads,
        "short_steps": bench.short_steps,
        "long_steps": bench.long_steps,
        "t_short": t_short,
        "t_long": t_long,
        "time_total": tl,
        "time_per_step": t_step,
        "time_setup": ts - bench.short_steps * t_step,
        "valid": not problems,
        "problems": problems,
    }


def add_speedups(results):
    """Speedup and parallel efficiency relative to the smallest resource count
    of the same benchmark (over nproc*threads)."""
    by_bench = {}
    for r in results:
        by_bench.setdefault(r["benchmark"], []).append(r)
    for rs in by_bench.values():
        base = min(rs, key=lambda r: r["nproc"] * r["threads"])
        base_units = base["nproc"] * base["threads"]
        for r in rs:
            units = r["nproc"] * r["threads"] / base_units
            for key in ("time_per_step", "time_total"):
                r["speedup_" + key[5:]] = base[key] / r[key]
            r["efficiency_per_step"] = r["speedup_per_step"] / units


def print_table(results, out=sys.stdout):
    hdr = "%-16s %5s %4s %10s %11s %9s %8s %7s  %s" % (
        "benchmark", "ranks", "thr", "setup [s]", "step [ms]", "total [s]", "speedup", "effic.", "valid")
    out.write(hdr + "\n" + "-" * len(hdr) + "\n")
    for r in results:
        out.write("%-16s %5d %4d %10.2f %11.2f %9.2f %8.2f %7.0f%%  %s\n" % (
            r["benchmark"], r["nproc"], r["threads"], r["time_setup"], 1e3 * r["time_per_step"],
            r["time_total"], r["speedup_per_step"], 100 * r["efficiency_per_step"],
            "yes" if r["valid"] else "INVALID"))
    for r in results:
        for p in r["problems"]:
            out.write("  INVALID %s np=%d thr=%d: %s\n" % (r["benchmark"], r["nproc"], r["threads"], p))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exe", default=os.environ.get("FAR3D_EXE"), help="far3d.x to benchmark (or FAR3D_EXE)")
    p.add_argument("--mpiexec", default=os.environ.get("FAR3D_MPIEXEC") or "mpiexec")
    p.add_argument("--mpiexec-args", nargs=argparse.REMAINDER, default=[],
                   help="extra launcher options, e.g. --bind-to core (must be the last option)")
    p.add_argument("--benchmarks", nargs="+", default=list(BENCHMARKS), choices=list(BENCHMARKS))
    p.add_argument("--np", nargs="+", type=int, default=[1], help="MPI rank counts")
    p.add_argument("--threads", nargs="+", type=int, default=[1], help="OMP_NUM_THREADS values")
    p.add_argument("--repeat", type=int, default=3)
    p.add_argument("--label", default=socket.gethostname().split(".")[0])
    p.add_argument("--arch", default="cpu", help="free-form architecture tag (cpu, gpu, ...)")
    p.add_argument("--output", help="results JSON (default: results/<label>.json)")
    p.add_argument("--workdir", default=str(h.REPO / "tests" / "runs" / "performance"))
    p.add_argument("--timeout", type=int, default=3600, help="per-run timeout [s]")
    p.add_argument("--update-reference", action="store_true",
                   help="rewrite the nonlinear reference from this run (do this once, on a trusted build)")
    p.add_argument("--list", action="store_true", help="describe the benchmarks and exit")
    args = p.parse_args(argv)

    if args.list:
        for b in BENCHMARKS.values():
            print("%s\n  %s\n  steps: %d (short), %d (long)\n  parallelism: %s\n" %
                  (b.name, b.description, b.short_steps, b.long_steps, b.parallelism))
        return 0
    if not args.exe:
        p.error("--exe (or FAR3D_EXE) is required")
    os.environ["FAR3D_MPIEXEC"] = args.mpiexec

    meta = metadata(args)
    print("FAR3d benchmarks: %s (%s), %s, commit %s" % (meta["label"], meta["arch"], meta["cpu"], meta["git_commit"]))
    results = []
    for name in args.benchmarks:
        bench = BENCHMARKS[name]
        print("\n%s: %s\n  parallelism: %s" % (name, bench.description, bench.parallelism))
        for threads in args.threads:
            for nproc in args.np:
                print("  ranks=%d threads=%d" % (nproc, threads))
                results.append(time_config(bench, nproc, threads, args, lambda s: print(s, flush=True)))
    add_speedups(results)

    print()
    print_table(results)
    output = Path(args.output) if args.output else HERE / "results" / ("%s.json" % args.label)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps({"metadata": meta, "results": results}, indent=1) + "\n")
    print("\nwrote %s" % output)
    return 0 if all(r["valid"] for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
