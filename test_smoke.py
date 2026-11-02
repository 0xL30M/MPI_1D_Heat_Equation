#!/usr/bin/env python3
"""Local correctness smoke tests."""
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
ENV = os.environ.copy()
ENV.setdefault("OMP_NUM_THREADS", "4")
FULL = os.environ.get("FULL", "0") == "1"


def run(binary, args, cwd, mpi=False):
    launcher = shlex.split(os.environ.get("MPIEXEC", "mpirun -np 4")) if mpi else []
    result = subprocess.run(launcher + [str(ROOT / binary)] + list(map(str, args)),
                            cwd=cwd, env=ENV, text=True, capture_output=True,
                            timeout=float(os.environ.get("TEST_TIMEOUT", "300")))
    if result.returncode:
        raise RuntimeError(f"{binary} exited {result.returncode}\n{result.stderr}")
    return result.stdout


def compare(actual, expected):
    a, b = actual.split(), expected.split()
    if not a or len(a) != len(b):
        raise ValueError(f"output length mismatch: {len(a)} vs {len(b)} tokens")
    for index, (x, y) in enumerate(zip(a, b)):
        x, y = float(x), float(y)
        if not (math.isfinite(x) and math.isfinite(y)
                and math.isclose(x, y, rel_tol=5e-5, abs_tol=1e-9)):
            raise ValueError(f"numeric mismatch at token {index}: {x} vs {y}")


def main():
    with tempfile.TemporaryDirectory(prefix="mcp-smoke-") as temporary:
        work = Path(temporary)
        cases = ((1, 10000), (2, 10001), (3, 50000), (4, 200000))
        for index, points in (cases if FULL else cases[:1]):
            fixture = ROOT / "tests" / str(index)
            output = work / "heat-output"
            output.unlink(missing_ok=True)
            run("ex01.app", [points, fixture / "input", output], work, mpi=True)
            actual = output.read_text().splitlines()
            expected = (fixture / "output").read_text().splitlines()
            if len(actual) != points + 4 or len(expected) != points + 4:
                raise ValueError("invalid heat output line count")
            # The last line is elapsed time, which varies between executions.
            compare("\n".join(actual[:-1]), "\n".join(expected[:-1]))
            print(f"PASS: heat fixture {index} ({points} points)")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        sys.exit(1)
