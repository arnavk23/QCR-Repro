from __future__ import annotations

import argparse
import csv
import statistics
import time
from pathlib import Path

from qcr_repro.gates import circuit_unitary
from qcr_repro.qasm import parse_qasm_subset, write_qasm_subset
from qcr_repro.reducer import reduce_with_lookup
from qcr_repro.token_pool import TokenPool
from qcr_repro.unitary import equivalent_up_to_global_phase


def parse_int_list(raw: str) -> list[int]:
    return [int(item.strip()) for item in raw.split(",") if item.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark the Python MATLAB-demo reducer across multiple seeds/configurations.")
    parser.add_argument("--input", required=True, help="Input QASM file")
    parser.add_argument("--output-dir", default="results/demo_sweep/strict", help="Directory for CSV and optional reduced outputs")
    parser.add_argument("--iters", default="500,1000,1500", help="Comma-separated iteration counts")
    parser.add_argument("--depths", default="3,4", help="Comma-separated graph depths")
    parser.add_argument("--seeds", default="1,5,10", help="Comma-separated RNG seeds")
    parser.add_argument("--local-qubits", type=int, default=3)
    parser.add_argument("--max-block", type=int, default=7)
    parser.add_argument("--write-reduced", action="store_true", help="Write reduced QASM for each run")
    parser.add_argument("--atol", type=float, default=1e-5, help="Equivalence tolerance")
    args = parser.parse_args()

    iters_list = parse_int_list(args.iters)
    depth_list = parse_int_list(args.depths)
    seed_list = parse_int_list(args.seeds)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    num_qubits, input_gates = parse_qasm_subset(args.input)
    # reduce_with_lookup snaps QASM's 4-decimal-rounded angles to the pool's
    # exact values before reducing (TokenPool.snap); verify against that
    # snapped reference rather than the raw QASM unitary, so the ~1e-4 QASM
    # rounding artifact isn't mistaken for a reduction error (same
    # convention used to validate the exact engine's demo-circuit result).
    snap_pool = TokenPool(num_qubits=1, gate_set="ion_trap")
    snapped_gates = [snap_pool.snap(gate) for gate in input_gates]
    input_u = circuit_unitary(num_qubits, snapped_gates)

    rows: list[dict[str, object]] = []

    for depth in depth_list:
        for iters in iters_list:
            for seed in seed_list:
                start = time.time()
                reduced_gates, stats = reduce_with_lookup(
                    input_gates,
                    num_qubits=num_qubits,
                    local_qubits=args.local_qubits,
                    max_block_len=args.max_block,
                    graph_depth=depth,
                    iterations=iters,
                    seed=seed,
                )
                runtime = time.time() - start

                reduced_u = circuit_unitary(num_qubits, reduced_gates)
                equivalent = equivalent_up_to_global_phase(input_u, reduced_u, atol=args.atol)

                reduction_ratio = len(reduced_gates) / len(input_gates)

                row = {
                    "depth": depth,
                    "iterations": iters,
                    "seed": seed,
                    "start_len": len(input_gates),
                    "end_len": len(reduced_gates),
                    "replacements": stats.replacements,
                    "runtime_sec": round(runtime, 4),
                    "equivalent": equivalent,
                    "reduction_ratio": round(reduction_ratio, 4),
                }
                rows.append(row)

                if args.write_reduced:
                    out_file = output_dir / f"reduced_d{depth}_i{iters}_s{seed}.qasm"
                    write_qasm_subset(out_file, num_qubits, reduced_gates)

                print(
                    f"depth={depth} iters={iters} seed={seed} | "
                    f"{len(input_gates)}->{len(reduced_gates)} gates | "
                    f"runtime={runtime:.2f}s | eq={equivalent}"
                )

    csv_path = output_dir / "benchmark_reducer.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=[
                "depth",
                "iterations",
                "seed",
                "start_len",
                "end_len",
                "replacements",
                "runtime_sec",
                "equivalent",
                "reduction_ratio",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    end_lengths = [int(row["end_len"]) for row in rows]
    runtimes = [float(row["runtime_sec"]) for row in rows]
    eq_rate = sum(bool(row["equivalent"]) for row in rows) / len(rows)

    summary_path = output_dir / "benchmark_summary.txt"
    with summary_path.open("w", encoding="utf-8") as file:
        file.write("Benchmark Summary\n")
        file.write(f"Runs: {len(rows)}\n")
        file.write(f"Input gates: {len(input_gates)}\n")
        file.write(f"Best end_len: {min(end_lengths)}\n")
        file.write(f"Mean end_len: {statistics.mean(end_lengths):.2f}\n")
        file.write(f"Mean runtime_sec: {statistics.mean(runtimes):.2f}\n")
        file.write(f"Median runtime_sec: {statistics.median(runtimes):.2f}\n")
        file.write(f"Equivalence pass rate: {eq_rate:.3f}\n")

    print(f"Saved CSV: {csv_path}")
    print(f"Saved summary: {summary_path}")


if __name__ == "__main__":
    main()
