#!/usr/bin/env python3
"""Compare FAR3d benchmark result files side by side.

  python3 compare_results.py results/m1-cpu.json results/a100.json [...]

For every (benchmark, ranks, threads) configuration present in any file,
prints the time per step and setup time from each file, plus the per-step
speed of each file relative to the first one (>1 = faster than the first).
"""

import json
import sys
from pathlib import Path


def main(paths):
    if len(paths) < 1:
        print(__doc__)
        return 2
    runs = [json.loads(Path(p).read_text()) for p in paths]
    labels = ["%s(%s)" % (r["metadata"]["label"], r["metadata"]["arch"]) for r in runs]

    for label, run in zip(labels, runs):
        m = run["metadata"]
        gpu = (", " + m["gpu"]) if m.get("gpu") else ""
        compiler = " ".join(m["build"].get("compiler_version", [])) or "compiler unknown"
        print("%-20s %s%s | %s | commit %s | %s" % (label, m["cpu"], gpu, compiler, m["git_commit"], m["date"]))
    print()

    keys = sorted({(r["benchmark"], r["nproc"], r["threads"]) for run in runs for r in run["results"]})
    col = max(14, max(len(l) for l in labels))
    head = "%-16s %5s %4s " % ("benchmark", "ranks", "thr")
    head += " ".join("%*s" % (col, l[:col]) for l in labels)
    print("time per step [ms] (relative speed vs. %s)" % labels[0])
    print(head)
    print("-" * len(head))
    for key in keys:
        cells, base = [], None
        for i, run in enumerate(runs):
            match = [r for r in run["results"] if (r["benchmark"], r["nproc"], r["threads"]) == key]
            if not match:
                cells.append("%*s" % (col, "-"))
                continue
            r = match[0]
            t = 1e3 * r["time_per_step"]
            if i == 0:
                base = t
            rel = "%.2fx" % (base / t) if base else "n/a"
            flag = "" if r["valid"] else " INVALID"
            cells.append("%*s" % (col, "%.2f (%s)%s" % (t, rel, flag)))
        print("%-16s %5d %4d " % key + " ".join(cells))

    print("\nsetup time [s]")
    print(head)
    print("-" * len(head))
    for key in keys:
        cells = []
        for run in runs:
            match = [r for r in run["results"] if (r["benchmark"], r["nproc"], r["threads"]) == key]
            cells.append("%*s" % (col, "%.2f" % match[0]["time_setup"] if match else "-"))
        print("%-16s %5d %4d " % key + " ".join(cells))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
