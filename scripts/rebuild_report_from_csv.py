"""Regenerate comparison_<gateset>_report.{md,json} from an already-written
CSV, without re-running the benchmark. Used after hand-merging the new
fast-sweep exact_len/exact_cost rows with the original (unaffected)
qiskit/BQSKit baseline rows into results/comparison/comparison_ion_trap.csv.

Usage:
    PYTHONPATH=src python scripts/rebuild_report_from_csv.py ion_trap
"""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_comparison import GATE_TYPES, _build_report, _stats_for  # noqa: E402


def main() -> None:
    gateset = sys.argv[1] if len(sys.argv) > 1 else "ion_trap"
    outdir = Path("results/comparison")
    csv_path = outdir / f"comparison_{gateset}.csv"

    results = []
    with csv_path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            results.append({
                "seed": int(row["seed"]),
                "method": row["method"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "counts": {
                    "RX": int(row["RX"]), "RY": int(row["RY"]), "RZ": int(row["RZ"]),
                    "RXX": int(row["RXX"]), "CZ": int(row["CZ"]),
                },
                "twq": int(row["twq"]),
                "secs": float(row["runtime_s"]),
                "ok": row["ok"] == "True",
                "verifier": row["verifier"],
            })

    methods = sorted({r["method"] for r in results},
                      key=lambda m: (m not in ("exact_len", "exact_cost", "numeric_len", "numeric_cost"), m))
    stats = {m: _stats_for(results, m, GATE_TYPES[gateset]) for m in methods}
    budgets = {r["secs"] for r in results if r["method"] in ("exact_len", "exact_cost", "numeric_len", "numeric_cost")}
    meta = {
        "gateset": gateset,
        "num_circuits": len({r["seed"] for r in results if r["method"] in methods[:1]}) or 100,
        "num_qubits": 4,
        "length": 300,
        "budget": 30.0 if gateset == "ion_trap" else 60.0,
        "depths": {},
        "max_block_len": None,
        "rz_pass": gateset == "nisq",
        "weights": None,
        "restarts": 1,
        "backend": "ram",
        "rf_gate": False,
        "hybrid": False,
        "fast_sweep": True,
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "wall_sec": 0.0,
        "with_numeric": False,
        "method_order": methods,
    }
    report = _build_report(gateset, results, stats, meta)
    (outdir / f"comparison_{gateset}_report.md").write_text(report, encoding="utf-8")
    (outdir / f"comparison_{gateset}_report.json").write_text(
        json.dumps({"meta": meta, "methods": {m: s for m, s in stats.items()}}, indent=2),
        encoding="utf-8",
    )
    print(f"rebuilt report for {gateset}: {len(results)} rows, methods {methods}")


if __name__ == "__main__":
    main()
