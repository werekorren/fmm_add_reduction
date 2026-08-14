"""parallel_run_3x3.py - run `3x3.py` many times in
parallel, varying only the random seed, and classify each run.

Classification categories:
1. failed_to_converge: never enabled ternary regularisation.
2. failed_to_ternarize: ternary enabled but final integer-weight verification failed.
3. success: ternary enabled and bilinear algorithm verified.

The script keeps completely hands-off w.r.t hyper-parameters; it merely patches the
PyTorch random seed used inside 3x3.py before importing/executing the script.
That avoids modifying the original file.
"""
from __future__ import annotations
import argparse
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Tuple
from process_success_runs import success_log_text_2_txt_file
import datetime

sys.stdin.reconfigure(encoding='utf-8')
sys.stdout.reconfigure(encoding='utf-8')

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "3x3.py"


def _worker(seed: int) -> Tuple[int, str, str]:
    """Run one training instance with the given *seed*.

    Returns (seed, classification, full_output).
    """
    code = f"""
import runpy, torch, sys, pathlib
seed = {seed}
# Ensure the script's directory is on the import path so that `model.py`, etc. are discoverable
script_path = pathlib.Path(r"{SCRIPT}")
sys.path.insert(0, str(script_path.parent))

# Patch torch.manual_seed so calls inside 3x3.py use *this* seed
_orig_manual_seed = torch.manual_seed

def _patched(_=None):
    return _orig_manual_seed(seed)

torch.manual_seed = _patched
_orig_manual_seed(seed)  # set seed before script begins
runpy.run_path(str(script_path), run_name="__main__")
"""
    # Execute in a fresh interpreter so global namespace is clean
    proc = subprocess.run(
        [sys.executable, "-u", "-c", code],
        capture_output=True,
        text=True,
        cwd=SCRIPT.parent,  # run from the script's directory to further ease relative imports
    )
    output = proc.stdout + proc.stderr

    ternary_enabled = "Enabled ternary regularisation" in output
    verified = "Bilinear algorithm verified" in output

    if not ternary_enabled:
        status = "failed_to_converge"
    elif not verified or proc.returncode != 0:
        status = "failed_to_ternarize"
    else:
        status = "success"

    return seed, status, output


def main() -> None:
    parser = argparse.ArgumentParser(description="Run 3x3.py many times in parallel and classify outcomes.")
    parser.add_argument("-n", "--runs", type=int, default=8, help="Number of runs (different seeds) to launch. Default: 8")
    parser.add_argument("--start", type=int, default=0, help="Starting seed value (inclusive). Default: 0")
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=os.cpu_count() or 1,
        help="Maximum concurrent processes. Default: #CPU cores",
    )
    parser.add_argument(
        "--logdir",
        type=Path,
        default=Path("3x3_run_logs"),
        help="Directory to store individual run logs.",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default="mmul-3x3/",
        help="Prefix (folder) within the bucket for uploaded logs.",
    )
    args = parser.parse_args()

    args.logdir.mkdir(parents=True, exist_ok=True)

    def run_batch(start_seed: int) -> int:
        """Run one batch starting at *start_seed* (inclusive). Returns next start_seed."""
        
        # print(f"\nLaunching batch: seeds {start_seed} – {start_seed + args.runs - 1} with up to {args.jobs} parallel jobs…\n")

        futures = {}
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            for seed in range(start_seed, start_seed + args.runs):
                fut = pool.submit(_worker, seed)
                futures[fut] = seed

            results = {
                "success": [],
                "failed_to_ternarize": [],
                "failed_to_converge": [],
            }

            for fut in as_completed(futures):
                seed, status, out = fut.result()
                results[status].append(seed)

                # Save log
                if status == "success":
                    success_log_text_2_txt_file(seed, out)
                    #timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
                    #log_path = args.logdir / f"seed_{seed:04d}_{status}_{timestamp}.log"
                    #log_path.write_text(out)
                    #print(f"Seed {seed:4d}: {status}  → log saved to {log_path}")

        # Summary
        print("\nBatch summary:")
        for k, seeds in results.items():
            print(f"  {k:20s}: {len(seeds):3d}  seeds → {sorted(seeds)}")

        return start_seed + args.runs

    # Continuous batches until user aborts
    first_seed = next_seed = args.start
    first_time = datetime.datetime.now()
    try:
        while True:
            next_seed = run_batch(next_seed)
            # print stats
            total_time = datetime.datetime.now() - first_time
            num_seeds = next_seed - first_seed
            time_per_seed = total_time / num_seeds
            seconds_per_seed = time_per_seed.total_seconds()
            seeds_per_day = 24 * 60 * 60 / seconds_per_seed
            print(f"next seed = {next_seed}, {num_seeds:6} seeds processed, total time = {total_time}, time per seed = {time_per_seed}, seeds per day = {seeds_per_day:4.0f}")
    except KeyboardInterrupt:
        print("\nInterrupted by user. Exiting.")


if __name__ == "__main__":
    main()