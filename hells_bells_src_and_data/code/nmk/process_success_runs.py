# process_success_runs.py
"""Automate post-processing of successful 3x3 training runs.

This script performs four high-level stages that can be toggled on/off by
setting the RUN_* flags in the __main__ section:

1. scan_logs:   locate *_success_*.log files and extract the ternary weight
                block printed by 3x3.py (as produced by verification.snap()).
                The matrices are saved to algorithms/other/ in the required
                format, and their naïve additive cost is computed via
                additive_cost() from verification.py.

2. run_reducer: invoke the compiled C binary ``fmm_add_reduction`` twice on
                every matrix file - once with mode "v" (vanilla greedy) and
                once with mode "p" (potential greedy).  The stdout/stderr of
                each invocation are captured into reduction_logs/.

3. summarise:   summarise the additive cost of every matrix file that exists
                in algorithms/other/.

All heavy-lifting lives in small, composable helper functions so any single
stage can be re-used easily from a notebook / REPL.

"""
from __future__ import annotations
import re
import subprocess
from pathlib import Path
from typing import Tuple, List
import numpy as np
import sys
import io
sys.stdin.reconfigure(encoding='utf-8')
sys.stdout.reconfigure(encoding='utf-8')

# Resolve project root relative to this file so paths work from any CWD
ROOT_DIR = Path(__file__).resolve().parent

# Ensure project root is on the import path for absolute imports
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

# Default to export directory containing successful runs
print(ROOT_DIR)
LOG_DIR = ROOT_DIR / "3x3_run_logs"
MATRIX_DIR = ROOT_DIR / "successful_algorithms"
REDUCTION_LOG_DIR = ROOT_DIR / "reduced_successful_algorithms_2"
FMM_REDUCER = ROOT_DIR / "fmm_add_reduction.exe"  # compiled binary

# Regex shortcuts
# New logs are plain per-task exports; treat every *.log as a success.
LOG_GLOB = "*.log"
BLOCK_RE = re.compile(
    r"----- paste-ready weight block -----\n(.*?)\n----------- end block --------------",
    re.DOTALL,
)
# Seed line emitted by entrypoint.py: "=== Running seed <N> ==="
#SEED_LINE_RE = re.compile(r"===\s*Running seed\s+(\d+)\s*===")
SEED_LINE_RE = re.compile(r"\d+")

SUCCESS_RE = re.compile(r"success")

# New regex patterns for logs that only contain Python literal arrays (without the
# paste-ready block).  We capture the *inner* list-of-lists so that we can feed
# it straight into ``ast.literal_eval``.
WA_ARR_RE = re.compile(r"WA\s*=\s*np\.array\((\[[\s\S]*?\])\s*,\s*dtype=np\.int8", re.MULTILINE)
WB_ARR_RE = re.compile(r"WB\s*=\s*np\.array\((\[[\s\S]*?\])\s*,\s*dtype=np\.int8", re.MULTILINE)
WC_ARR_RE = re.compile(r"WC\s*=\s*np\.array\((\[[\s\S]*?\])\s*,\s*dtype=np\.int8", re.MULTILINE)

from verification import additive_cost  # type: ignore


def _parse_weight_block(block: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (WA, WB, WC) from the triple-section *block* string.

    The incoming block has the format emitted by verification.snap():
        WA^T rows (space-separated ints)
        "#"
        WB^T rows
        "#"
        WC rows
    """
    sections = [s.strip() for s in block.strip().split("#")]
    assert len(sections) == 3, "Expected exactly three sections in weight block"

    def _to_int_matrix(txt: str) -> np.ndarray:
        rows = [[int(tok) for tok in line.split()] for line in txt.splitlines() if line.strip()]
        return np.array(rows, dtype=np.int8)

    WA_T = _to_int_matrix(sections[0])
    WB_T = _to_int_matrix(sections[1])
    WC = _to_int_matrix(sections[2])

    WA = WA_T.T  # (r, n^2)
    WB = WB_T.T
    # WC is already (n^2, r)

    return WA, WB, WC


def _success_log_text_2_txt_file_output(text: str) -> str:
    m = BLOCK_RE.search(text)
    if m:
        # Classic paste-ready block present: use it as-is.
        block = m.group(1)
    else:
        # Fallback to Python literal arrays (newer log format).
        import ast  # local import to avoid polluting global namespace unnecessarily

        wa_m = WA_ARR_RE.search(text)
        wb_m = WB_ARR_RE.search(text)
        wc_m = WC_ARR_RE.search(text)

        if not (wa_m and wb_m and wc_m):
            raise ValueError(f"Could not find weight arrays in {log_path}")

        # Safely evaluate the captured list-of-lists strings.
        WA = np.array(ast.literal_eval(wa_m.group(1)), dtype=np.int8)
        WB = np.array(ast.literal_eval(wb_m.group(1)), dtype=np.int8)
        WC = np.array(ast.literal_eval(wc_m.group(1)), dtype=np.int8)

        # Reconstruct the optimiser-friendly text block expected downstream.
        WA_block = "\n".join(" ".join(map(str, row)) for row in WA.T)
        WB_block = "\n".join(" ".join(map(str, row)) for row in WB.T)
        WC_block = "\n".join(" ".join(map(str, row)) for row in WC)
        block = f"{WA_block}\n#\n{WB_block}\n#\n{WC_block}"

    return block



def success_log_text_2_txt_file(n:int, m:int, k:int, r:int, seed : int, text: str):
    block = _success_log_text_2_txt_file_output(text)
    WA, WB, WC = _parse_weight_block(block)
    naive_adds = int(additive_cost(WA, WB, WC))
    out_name = f"seed{seed:06d}-{n}{m}{k}-{r}-{naive_adds}.txt"
    out_path = MATRIX_DIR / out_name
    try:
        out_path.write_text(block + "\n", encoding="utf-8")
        #print(f"Extracted seed {seed:4d} -> {out_name}  (additions={naive_adds})")
    except Exception as e:
        print(f"[WARN] Failed to write output to file {out_name}: {e}")



def _extract_success_log(log_path: Path) -> Tuple[int, str]:
    """Return (seed, weight_block) from a successful run log."""
    success_log = SUCCESS_RE.search(log_path.name)
    if not success_log:
      return

    text = log_path.read_text(encoding="utf-8", errors="replace")
    block = _success_log_text_2_txt_file_output(text)

    seed_line = SEED_LINE_RE.search(log_path.name)
    seed = int(seed_line.group(0))
#    print(seed)

#    seed_line = SEED_LINE_RE.search(text)
    seed_line = SEED_LINE_RE.search(log_path.name)
#    print(int(seed_line.group(0)))
    if not seed_line:
        raise ValueError(f"Failed to find seed banner in {log_path.name}")
#    seed = int(seed_line.group(1))
    seed = int(seed_line.group(0))
    return seed, block


def scan_logs(log_dir: Path = LOG_DIR, matrix_dir: Path = MATRIX_DIR) -> List[Path]:
    """Extract matrices from all *_success_*.log files.

    Returns list of newly created *.txt files."""
    matrix_dir.mkdir(parents=True, exist_ok=True)
    created: List[Path] = []

    for log_path in sorted(log_dir.glob(LOG_GLOB)):
        try:
            seed, block = _extract_success_log(log_path)
            WA, WB, WC = _parse_weight_block(block)
            naive_adds = int(additive_cost(WA, WB, WC))
            out_name = f"seed{seed:06d}-333-23-{naive_adds}.txt"
            out_path = matrix_dir / out_name
            out_path.write_text(block + "\n", encoding="utf-8")
            print(f"Extracted seed {seed:4d} -> {out_path}  (additions={naive_adds})")
            created.append(out_path)
        except Exception as e:
            print(f"[WARN] Skipping {log_path.name}: {e}")

    return created


def _run_single_reduction(txt_file: Path, mode: str, log_dir: Path = REDUCTION_LOG_DIR) -> Path:
    """Run fmm_add_reduction on *txt_file* with given *mode* (v / p)."""
    assert mode in {"v", "p"}, "mode must be 'v' or 'p'"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{txt_file.stem}.{mode}.log"

    if mode == "p":
        mode = "p 0 .5 50"
    cmd = [str(FMM_REDUCER), str(txt_file), mode]
    print("$", " ".join(cmd))
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
        output, _ = proc.communicate()
        log_path.write_text(output, encoding="utf-8")
        if proc.returncode != 0:
            print(f"  -> reducer exited with code {proc.returncode}")
    return log_path


def run_reducer(matrix_dir: Path = MATRIX_DIR) -> None:
    """Invoke the C++ reducer on every .txt matrix found."""
    if not FMM_REDUCER.is_file():
        print(f"Reducer binary not found at {FMM_REDUCER}. Skipping reduction stage.")
        return
    txt_files = sorted(matrix_dir.glob("seed*-333-23-*.txt"))
    if not txt_files:
        print("No matrix files found – run scan_logs() first.")
        return

    for txt in txt_files:
        for mode in ("v", "p"):
            _run_single_reduction(txt, mode)


def summarise(matrix_dir: Path = MATRIX_DIR, log_dir: Path = REDUCTION_LOG_DIR) -> None:
    """Print python-vs-CPP additive-cost summary for each matrix file.

    For every `seed*-333-23-*.txt` file we look for matching reduction logs:
      reduction_logs/<stem>.v.log and .p.log

    We report:
        • naive additions computed in Python (from filename)
        • naive additions parsed from the CPP logs (one per mode)
        • reduced additions after reduction (per mode)
    """

    # Regex patterns inside reduction logs
    RE_CPP_NAIVE = re.compile(r"Algorithm uses\s+\d+\s+\+\s+\d+\s+\+\s+\d+\s+=\s+(\d+)\s+additions \[naive\]")
    RE_CPP_REDUCED = re.compile(r"Algorithm uses\s+\d+\s+\+\s+\d+\s+\+\s+\d+\s+=\s+(\d+)\s+additions after reduction")

    def _parse_cpp_numbers(log_path: Path) -> tuple[int | None, int | None]:
        if not log_path.is_file() or log_path.stat().st_size == 0:
            return None, None
        text = log_path.read_text(errors="replace")
        naive = RE_CPP_NAIVE.search(text)
        reduced = RE_CPP_REDUCED.search(text)
        return (int(naive.group(1)) if naive else None, int(reduced.group(1)) if reduced else None)

    header = (
        f"{'matrix':35s}  {'py':>4s}  {'v':>4s} {'v_red':>5s}  "
        f"{'p':>4s} {'p_red':>5s}  {'delta_v':>3s} {'delta_p':>3s}"
    )
    print("\n" + header)
    print("-" * len(header))

    for txt in sorted(matrix_dir.glob("seed*-333-23-*.txt")):
        # python naive from filename
        m_cost = re.search(r"-(\d+)\.txt$", txt.name)
        py_naive = int(m_cost.group(1)) if m_cost else None

        stem = txt.stem  # without .txt
        v_log = log_dir / f"{stem}.v.log"
        p_log = log_dir / f"{stem}.p.log"

        v_naive, v_red = _parse_cpp_numbers(v_log)
        p_naive, p_red = _parse_cpp_numbers(p_log)

        # Differences
        dv = (py_naive - v_red) if (py_naive is not None and v_red is not None) else None
        dp = (py_naive - p_red) if (py_naive is not None and p_red is not None) else None

        def _fmt(x: int | None) -> str:
            return f"{x:4d}" if x is not None else "----"

        print(
            f"{txt.name:35s}  "
            f"{_fmt(py_naive)}  {_fmt(v_naive)} {_fmt(v_red)}  "
            f"{_fmt(p_naive)} {_fmt(p_red)}  "
            f"{_fmt(dv) if dv is not None else ' --'} {_fmt(dp) if dp is not None else ' --'}"
        )

        # Track minimal reduced scores (allowing ties across modes)
        for mode_label, red in (("v", v_red), ("p", p_red)):
            if red is None:
                continue
            if 'min_red' not in locals() or red < min_red:
                # new record – start fresh list
                min_red = red
                winners = [(txt.name, mode_label, red)]
            elif red == min_red:
                # tie – append
                winners.append((txt.name, mode_label, red))

    if 'min_red' in locals():
        print("\nLowest reduced-addition run(s):")
        for name, mode_label, score in winners:
            print(f"  {name}  ({mode_label}) → {score} additions after reduction")
    else:
        print("\nNo reduction logs with parsed results found.")


if __name__ == "__main__":
    import sys
    logfile = open(str(ROOT_DIR / "process_success_runs.log"), "w")
    class Tee:
        def __init__(self, *streams):
            self.streams = streams
        def write(self, data):
            for s in self.streams:
                s.write(data)
        def flush(self):
            for s in self.streams:
                s.flush()
    sys.stdout = Tee(sys.stdout, logfile)
    sys.stderr = Tee(sys.stderr, logfile)

    RUN_SCAN_LOGS = True     # extract matrices from *_success_*.log
    RUN_REDUCER = False       # run fmm_add_reduction on each matrix (v & p)
    RUN_SUMMARY = True       # print additive-cost summary at the end

    if RUN_SCAN_LOGS:
        scan_logs()

    if RUN_REDUCER:
        run_reducer()

    if RUN_SUMMARY:
        summarise() 