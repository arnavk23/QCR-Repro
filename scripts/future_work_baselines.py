from __future__ import annotations

import csv
import json
import time
from pathlib import Path

from qcr_repro.gates import circuit_unitary
from qcr_repro.qasm_io import parse_qasm_subset
from qcr_repro.reducer import reduce_with_lookup
from qcr_repro.unitary_utils import equivalent_up_to_global_phase

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "future_work"
LONG_QASM = ROOT / "paper_demo" / "QCOptimDemo" / "longcode10.txt"


def gate_count_ops(qc) -> int:
    return sum(int(v) for v in qc.count_ops().values())


def run_local_method(num_qubits, gates):
    start = time.time()
    reduced, stats = reduce_with_lookup(
        gates,
        num_qubits=num_qubits,
        local_qubits=3,
        max_block_len=7,
        graph_depth=4,
        iterations=1500,
        seed=10,
    )
    runtime = time.time() - start
    ok = equivalent_up_to_global_phase(circuit_unitary(num_qubits, gates), circuit_unitary(num_qubits, reduced), atol=1e-5)
    return {
        "method": "local_lookup_strict",
        "config": "depth=4,iters=1500,seed=10",
        "end_len": len(reduced),
        "runtime_sec": round(runtime, 4),
        "equivalent_1e5": ok,
        "notes": "Python implementation method",
    }


def run_qiskit_baselines(num_qubits, gates):
    rows = []
    try:
        from qiskit import QuantumCircuit, transpile
        from qiskit.circuit.library import RXGate, RYGate, RZGate, RXXGate

        qc = QuantumCircuit(num_qubits)
        for g in gates:
            if g.name == "RX":
                qc.append(RXGate(g.theta), [g.qubits[0]])
            elif g.name == "RY":
                qc.append(RYGate(g.theta), [g.qubits[0]])
            elif g.name == "RZ":
                qc.append(RZGate(g.theta), [g.qubits[0]])
            elif g.name == "RXX":
                qc.append(RXXGate(g.theta), [g.qubits[0], g.qubits[1]])
            else:
                raise ValueError(f"Unsupported gate: {g.name}")

        basis = ["rx", "ry", "rz", "rxx"]
        for level in [0, 1, 2, 3]:
            st = time.time()
            t = transpile(qc, basis_gates=basis, optimization_level=level)
            rt = time.time() - st
            rows.append(
                {
                    "method": f"qiskit_opt{level}",
                    "config": f"basis={','.join(basis)}",
                    "end_len": gate_count_ops(t),
                    "runtime_sec": round(rt, 4),
                    "equivalent_1e5": "N/A",
                    "notes": "Transpile baseline",
                }
            )
    except Exception as exc:
        rows.append(
            {
                "method": "qiskit_baseline",
                "config": "unavailable",
                "end_len": -1,
                "runtime_sec": -1,
                "equivalent_1e5": "N/A",
                "notes": f"Failed: {type(exc).__name__}: {exc}",
            }
        )
    return rows


def run_bqskit_baseline(num_qubits, gates):
    rows = []
    try:
        import bqskit
        rows.append(
            {
                "method": "bqskit_import_check",
                "config": f"version={getattr(bqskit, '__version__', 'unknown')}",
                "end_len": -1,
                "runtime_sec": -1,
                "equivalent_1e5": "N/A",
                "notes": "BQSKit installed; full pass parity script requires dedicated pass mapping for this gate set.",
            }
        )
    except Exception as exc:
        rows.append(
            {
                "method": "bqskit_baseline",
                "config": "unavailable",
                "end_len": -1,
                "runtime_sec": -1,
                "equivalent_1e5": "N/A",
                "notes": f"Failed: {type(exc).__name__}: {exc}",
            }
        )
    return rows


def hardware_metrics_probe():
    result = {
        "hardware_available": False,
        "provider": "none",
        "backend": "none",
        "reason": "No authenticated hardware provider configured in this environment.",
    }
    try:
        from qiskit_ibm_runtime import QiskitRuntimeService

        service = QiskitRuntimeService()
        backend = service.least_busy(simulator=False, operational=True)
        props = backend.properties()
        result = {
            "hardware_available": True,
            "provider": "qiskit_ibm_runtime",
            "backend": backend.name,
            "avg_t1_us": round(sum(q.t1 for q in props.qubits) / len(props.qubits) * 1e6, 3),
            "avg_t2_us": round(sum(q.t2 for q in props.qubits) / len(props.qubits) * 1e6, 3),
            "pending_jobs": backend.status().pending_jobs,
        }
    except Exception:
        pass
    return result


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    num_qubits, gates = parse_qasm_subset(LONG_QASM)

    rows = []
    rows.append(run_local_method(num_qubits, gates))
    rows.extend(run_qiskit_baselines(num_qubits, gates))
    rows.extend(run_bqskit_baseline(num_qubits, gates))

    csv_path = OUT_DIR / "baseline_parity_results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["method", "config", "end_len", "runtime_sec", "equivalent_1e5", "notes"],
        )
        writer.writeheader()
        writer.writerows(rows)

    hw = hardware_metrics_probe()
    hw_path = OUT_DIR / "hardware_metrics.json"
    hw_path.write_text(json.dumps(hw, indent=2), encoding="utf-8")

    md_path = OUT_DIR / "future_work_status.md"
    with md_path.open("w", encoding="utf-8") as f:
        f.write("# Future Work Integration Status\n\n")
        f.write("## Baseline parity\n")
        f.write(f"- Results CSV: {csv_path.name}\n")
        f.write("- Methods executed: local implementation + Qiskit baselines + BQSKit availability probe.\n\n")
        f.write("## Hardware-level metrics\n")
        if hw.get("hardware_available"):
            f.write(f"- Backend: {hw['backend']}\n")
            f.write(f"- Avg T1 (us): {hw['avg_t1_us']}\n")
            f.write(f"- Avg T2 (us): {hw['avg_t2_us']}\n")
            f.write(f"- Pending jobs: {hw['pending_jobs']}\n")
        else:
            f.write(f"- Not available: {hw.get('reason', 'unknown reason')}\n")
        f.write("\n## Note\n")
        f.write("This folder contains a concrete first implementation of the requested future-work direction.\n")

    print(f"Wrote {csv_path}")
    print(f"Wrote {hw_path}")
    print(f"Wrote {md_path}")


if __name__ == "__main__":
    main()
