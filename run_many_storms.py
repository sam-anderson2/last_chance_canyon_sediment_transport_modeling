#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_many_storms
===============

Purpose
-------
Dispatch all hyetographs in ``./working_climate/*.txt`` to multiple Python
processes, each calling ``driver_expt1a.run_simulation`` to produce its
own ``run_<stormname>/`` output directory. This script focuses on *robustness*
for long parallel batches on memory-constrained machines (e.g., Windows):

- Limits BLAS/NumPy threading **inside each worker** to reduce RAM spikes.
- Uses a conservative default for parallel workers (~one third of available cores).
- Detects NumPy ``ArrayMemoryError`` / "Unable to allocate ..." failures in
  parallel workers and **retries the same storm serially** in the main process.
  This avoids losing nearly-complete batches due to end-of-run fragmentation.

Typical usage
-------------
$ python run_many_storms.py
$ python run_many_storms.py --sim_max_t 25000
$ python run_many_storms.py --max_workers 4

---------------------------------------------------------------------------
HOW TO RUN FROM THE COMMAND LINE (step-by-step)
---------------------------------------------------------------------------

0) Prereqs
   - You have a working Python environment with your model dependencies.
   - Your rainfall hyetographs are in:  <folder with this script>/working_climate/*.txt
     (two-column TXT files; each file name becomes the storm tag, e.g., 3hr_100yr.txt)

1) Open a terminal
   - Windows (PowerShell):
       PS C:\path\to\experiment> cd C:\path\to\experiment
       PS C:\path\to\experiment> python --version
   - macOS/Linux:
       $ cd /path/to/experiment
       $ python3 --version

2) Run ALL storms with defaults
   - Windows:
       PS> python run_many_storms.py
   - macOS/Linux:
       $ python3 run_many_storms.py

   What happens:
   - Script scans ./working_climate for *.txt storms
   - Chooses default max_workers ≈ (CPU cores)/3
   - Launches storms in parallel; each storm creates its own run_<stormname>/ folder
   - Any OOM or similar memory error in a worker triggers **serial retry** at the end

3) Run ALL storms with a specific simulation time (e.g., 25,000 s)
   - PS> python run_many_storms.py --sim_max_t 25000
   - $ python3 run_many_storms.py --sim_max_t 25000

4) Override the parallel worker count (e.g., force 4 workers)
   - PS> python run_many_storms.py --max_workers 4
   - Combine with sim time:
       PS> python run_many_storms.py --sim_max_t 25000 --max_workers 4

5) Run ONLY specific storms (name or partial name match)
   - Exact names (without .txt):
       PS> python run_many_storms.py --storms 15min_2yr 3hr_1000yr
   - Partial substrings are allowed:
       PS> python run_many_storms.py --storms 15min 1000yr
     (This runs any file whose name contains "15min" OR "1000yr")

   You can COMBINE with other flags:
       PS> python run_many_storms.py --storms 3hr_2yr --sim_max_t 10000 --max_workers 2

6) Interpreting console output
   - Look for lines like:
       "✔️  [3hr_2yr] complete  →  C:\...\run_3hr_2yr (3/12)"
     which confirm per-storm success and show the output directory.
   - If you see:
       "❌  [<name>] failed in parallel worker ... ↳ will retry serially"
     the script will automatically re-run those storms one-by-one at the end.

7) Where do outputs go?
   - Each storm writes to:  <script folder>/run_<stormname>/
   - This script pins outputs to the script location (see `cwd=str(ROOT)`).

8) Troubleshooting tips
   - "❌ No *.txt hyetographs found": Ensure your files are in ./working_climate
   - OOM (Out-of-Memory): The script already limits math library threads and caps
     workers. If needed, try fewer workers:  --max_workers 1 or 2
   - Paths: Run the script from the folder where it lives, or ensure relative
     paths resolve (the script anchors ROOT to its own folder).

"""

# --- Set BLAS/NumPy threading to 1 **before** any heavy imports ---------------
# (This helps avoid per-process oversubscription and reduces memory spikes.)
# We set these env vars *before* importing any heavy numeric libs. This matters
# because libraries like NumPy/MKL/BLAS snapshot thread counts at import time.
import os as _os
_os.environ.setdefault("OMP_NUM_THREADS", "1")
_os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
_os.environ.setdefault("MKL_NUM_THREADS", "1")
_os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
_os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
_os.environ.setdefault("BLIS_NUM_THREADS", "1")
del _os  # keep global namespace tidy

# Standard library imports (safe after env tuning)
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import argparse
import os
import multiprocessing as mp
import traceback
import sys  # for flushing and worker-side traceback printing

# EDIT: pin all paths/outputs to this script's folder
# Using the file's parent ensures relative paths are deterministic regardless
# of the shell's current working directory when the script is launched.
ROOT = Path(__file__).parent.resolve()

# Project import: heavy work is done here
# The driver should expose: run_simulation(rainfall_file, tag, sim_max_t, cwd)
from driver_expt1a import run_simulation


# ------------------------------------------------------------------ #
def _oom_like(exc: BaseException | str) -> bool:
    """
    Heuristic: does an exception (or message) look like a NumPy/OS out-of-memory?

    We match common NumPy/Windows/Linux messages:
    - "ArrayMemoryError"
    - "Unable to allocate"
    - "bad allocation" (rare, but seen on some MKL builds)
    - "std::bad_alloc"

    Returns
    -------
    bool
        True if the error string smells like an out-of-memory condition.
        This helps us decide whether to attempt a safer *serial* retry.
    """
    msg = str(exc)
    needles = ("ArrayMemoryError", "Unable to allocate", "bad allocation", "std::bad_alloc")
    return any(k in msg for k in needles)


# ------------------------------------------------------------------ #
def _run_storm(storm_path: Path, sim_max_t: float) -> tuple[str, Path]:
    """
    Worker executed in a separate process.

    Parameters
    ----------
    storm_path : pathlib.Path
        Hyetograph file (two-column TXT).
    sim_max_t : float
        Simulation runtime in seconds.

    Returns
    -------
    tuple[str, Path]
        (storm tag, output directory)

    Notes
    -----
    * This runs inside a worker process when called via the pool.
    * We print start/finish messages and flush stdout to keep logs in order
      even when multiple workers are active.
    """
    tag = storm_path.stem  # e.g. "6hr_25yr"  → becomes run_<tag>/ output folder name

    try:
        print(f"→ Starting storm: {storm_path}")
        sys.stdout.flush()
        out_dir = run_simulation(
            rainfall_file=str(storm_path),
            tag=tag,
            sim_max_t=sim_max_t,
            cwd=str(ROOT),          # EDIT: pin outputs to this script's folder
        )
        print(f"✓ Finished storm: {storm_path}")
        sys.stdout.flush()
        return tag, out_dir
    except Exception as e:
        # IMPORTANT: Do not swallow exceptions here. We print details for
        # diagnostics, flush, and re-raise so the main process can decide how
        # to handle (e.g., mark for serial retry).
        print(f"\n⚠️  Storm {storm_path} failed in worker process!")
        print(f"    Exception type: {type(e).__name__}")
        print(f"    Message: {e}")
        print("    --- FULL TRACEBACK (worker) ---")
        traceback.print_exc(limit=8)
        sys.stdout.flush()
        sys.stderr.flush()
        raise


# ------------------------------------------------------------------ #
def _run_storm_serial_retry(storm_path: Path, sim_max_t: float) -> tuple[str, Path]:
    """
    Serial (in-process) retry for a single storm.

    This mirrors `_run_storm` but avoids a worker process. Useful if the
    parallel run died due to memory pressure—running serially often succeeds.
    """
    tag = storm_path.stem
    out_dir = run_simulation(
        rainfall_file=str(storm_path),
        tag=tag,
        sim_max_t=sim_max_t,
        cwd=str(ROOT),              # EDIT: pin outputs to this script's folder
    )
    return tag, out_dir


# ------------------------------------------------------------------ #
if __name__ == "__main__":
    # ----------------------- CLI / argument parsing ------------------------
    # The CLI lets you set a global sim time, cap workers, and filter which
    # storms to run by name. Storm names are derived from *.txt filenames.
    ap = argparse.ArgumentParser(
        description="Run every storm in ./working_climate/ on multiple processes (memory-safe)."
    )
    ap.add_argument(
        "--sim_max_t",
        type=float,
        default=4500.0,
        help="Total simulation time [s] applied to EVERY storm file",
    )
    ap.add_argument(
        "--max_workers",
        type=int,
        default=None,
        help="Optional cap for parallel workers (default ≈ one third of CPU cores).",
    )
    # NEW: optional filter to run specific storms by (partial) name
    ap.add_argument(
        "--storms",
        nargs="+",
        default=None,
        help="Optional: names or partial names of specific storms to run (e.g. --storms 15min_2yr 3hr_1000yr). "
             "Default: run all storms.",
    )
    args = ap.parse_args()

    # ----------------------- Discover storm files --------------------------
    # We anchor to the script folder so relative paths are stable even if
    # the user launches from another directory.
    storm_dir = ROOT / "working_climate"

    # Collect *.txt hyetographs; filenames are used as tags.
    storm_files = sorted(storm_dir.glob("*.txt"))
    if not storm_files:
        # Fail fast with a clear message if no inputs are present.
        raise SystemExit("❌  No *.txt hyetographs found in ./working_climate/")

    # ----------------------- Optional filtering by name --------------------
    # If --storms is provided, we keep only files whose names contain any of
    # the provided substrings. This allows exact names or partial matches.
    if args.storms:
        orig = storm_files
        storm_files = [f for f in storm_files if any(s in f.name for s in args.storms)]
        if not storm_files:
            raise SystemExit(f"❌  No rainfall files matched {args.storms}")
        print(f"🔍 Running only selected storms: {[f.name for f in storm_files]}")
    else:
        print(f"🌀  Running all storms in {storm_dir}")

    # ----------------------- Choose worker count safely --------------------
    # Default is ≈ one third of CPU cores to reduce RAM contention from
    # multiple Landlab/NumPy instances. We also never exceed the number of
    # storms to avoid idle workers.
    cpu = os.cpu_count() or 1
    default_workers = max(1, cpu // 3)  # ~one third of cores
    n_workers = args.max_workers if (args.max_workers and args.max_workers > 0) else default_workers
    n_workers = min(n_workers, len(storm_files))  # don't launch more workers than storms

    print(
        f"🌀  {len(storm_files)} storms detected  •  sim_max_t = {args.sim_max_t:.0f} s  •  "
        f"max_workers = {n_workers}\n"
    )

    # ----------------------- Process context (Windows) ---------------------
    # Windows uses "spawn" by default. We request it explicitly for clarity,
    # and to make behavior consistent across platforms when needed.
    ctx = mp.get_context("spawn")

    # ----------------------- Submit parallel jobs --------------------------
    # We submit all storms to a ProcessPoolExecutor and consume futures as they
    # complete. Any failure triggers a record in `failed` for later serial retry.
    print("▶️  Launching parallel runs...\n")
    futures = {}
    failed = []  # collect names to retry serially later
    try:
        with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as pool:
            for sf in storm_files:
                # Each future runs one storm in a worker process.
                fut = pool.submit(_run_storm, sf, args.sim_max_t)
                futures[fut] = sf

            completed = 0
            for fut in as_completed(futures):
                storm_path = futures[fut]
                name = storm_path.stem
                try:
                    # If the worker succeeded, unpack (tag, out_dir).
                    tag, out_dir = fut.result()
                    completed += 1
                    print(f"✔️  [{tag}] complete  →  {out_dir}   ({completed}/{len(storm_files)})")
                except Exception as e:
                    # We got a failure from a worker. Log details, classify
                    # whether it looks like an OOM, and mark the storm to be
                    # retried serially at the end.
                    print(f"❌  [{name}] failed in parallel worker ({completed+1}/{len(storm_files)})")
                    if _oom_like(e):
                        print("    ↳ looks like OOM; will retry serially")
                    else:
                        print("    ↳ non-OOM error; will retry serially")
                    print("    ↳ Exception type:", type(e).__name__)
                    print("    ↳ Message:", str(e))
                    tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
                    print("    ↳ Traceback (last frames):")
                    print(tb.splitlines()[-10:] if tb else "(no traceback)")
                    failed.append(name)
                finally:
                    # Ensure output is flushed promptly so logs are readable
                    # during long multi-hour runs.
                    sys.stdout.flush()
                    sys.stderr.flush()

    except KeyboardInterrupt:
        # Graceful Ctrl-C handling: cancel outstanding work and exit with the
        # conventional SIGINT code (130).
        print("\n⛔  Interrupted by user. Cancelling outstanding jobs...")
        for fut in futures:
            fut.cancel()
        raise SystemExit(130)  # conventional exit code for SIGINT

    # ----------------------- Serial retry phase ----------------------------
    # If any storms failed in parallel, we attempt a serial (in-process) run
    # for each. This often recovers from transient memory pressure.
    if failed:
        print("\n▶️  Re-running ONLY failed storms (serial, safer)…\n")
    else:
        print("\n✅  All storms finished in parallel.")
        raise SystemExit(0)

    # Serial retry for failed ones
    still_failed = []
    for name in failed:
        storm_path = storm_dir / f"{name}.txt"
        try:
            out_dir = run_simulation(
                rainfall_file=str(storm_path),
                tag=name,
                sim_max_t=args.sim_max_t,
                cwd=str(ROOT),          # EDIT: pin outputs to this script's folder
            )
            print(f"✔️  [{name}] serial retry complete  →  {out_dir}")
        except Exception as e:
            # Deep diagnostics for tricky failures: we print tracebacks and
            # also try to probe likely output folders to help you debug.
            print("\n[DEBUG] --- Serial retry exception diagnostics ---")
            print(f"Storm: {name}")
            print(f"Exception type: {type(e).__name__}")
            print(f"Message: {e}")
            tb = "".join(traceback.format_exception(type(e), e, e.__traceback__))
            print("[TRACEBACK START]")
            print(tb)
            print("[TRACEBACK END]")

            probe_path = None
            if isinstance(e, FileNotFoundError):
                # Some FileNotFoundError instances have a .filename attribute.
                # We attempt to surface it if available to guide debugging.
                try:
                    probe_path = Path(e.filename) if getattr(e, "filename", None) else None
                except Exception:
                    probe_path = None

            # Heuristic guesses for where partial outputs might live.
            guesses = []
            guesses.append(Path.cwd() / f"run_{name}")
            if probe_path is not None:
                guesses.append(Path(probe_path).parent)

            seen = set()
            for gpath in guesses:
                if gpath and gpath not in seen:
                    seen.add(gpath)
                    print(f"[DEBUG] Inspecting: {gpath}")
                    if gpath.exists():
                        try:
                            files = sorted(p.name for p in gpath.iterdir())
                            print(f"[DEBUG] Exists. Files present ({len(files)}): {files}")
                        except Exception as ee:
                            print(f"[DEBUG] Exists but cannot list files: {ee}")
                    else:
                        print("[DEBUG] Path does not exist.")

            if probe_path is not None:
                print(f"[DEBUG] Missing file path from exception: {probe_path}")
                print(f"[DEBUG] Parent exists? {probe_path.parent.exists()}")
                print(f"[DEBUG] File exists now? {probe_path.exists()}")

            print("[DEBUG] --- End diagnostics ---\n")
            still_failed.append(name)

    # ----------------------- Final status & exit code ----------------------
    if still_failed:
        print("\n❌  These storms still failed after serial retry:")
        for nm in still_failed:
            print(f"   - {nm}")
        # Non-zero exit signals partial failure to calling scripts/CI.
        raise SystemExit(1)

    print("\n✅  All storms finished (including serial retries).")
