"""Head-to-head comparison benchmark against the paper (NJP 27, 104509, 2025).

Runs the paper's statistical protocol (Tables 6/7): 100 random 4-qubit length-300 circuits per gate set under a fixed per-circuit time budget, our reducers vs qiskit/BQSKit baselines, with WIN/LOSE verdicts against the paper's reported numbers.

Usage:
    PYTHONPATH=src python scripts/benchmark_comparison.py --gateset ion_trap --num-circuits 100 --budget 30
    PYTHONPATH=src python scripts/benchmark_comparison.py --gateset nisq --budget 60 [--deep] [--hybrid] [--rf-gate] [--no-baselines]

Results (CSV + markdown + JSON) are written to --outdir."""

from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
import os
import statistics
import time
from pathlib import Path

from qcr_repro.circuits import count_gates, random_circuit
from qcr_repro.database import load_or_build_database
from qcr_repro.exact_database import load_or_build_exact
from qcr_repro.exact_reducer import reduce_circuit_exact, verify_exact
from qcr_repro.gates import circuit_unitary
from qcr_repro.hybrid import HybridDatabase
from qcr_repro.reducer import reduce_circuit
from qcr_repro.rf_gate import RfGate, RfGatedDatabase
from qcr_repro.unitary import equivalent_up_to_global_phase

# --------------------------------------------------------------------------- #
# configuration
# --------------------------------------------------------------------------- #

# Depth per number of wires a block may touch (wN dM = N wires, depth M).
# The ion-trap databases are already cached at these sizes; the NISQ default
# raises the 3-wire graph to depth 5 (the paper's main tool, ~140k nodes) and
# the 2-wire graph to depth 6 -- both are rebuilt once and cached on disk.
ION_DEPTHS = {1: 12, 2: 10, 3: 7, 4: 5}
NISQ_DEPTHS_DEFAULT = {1: 12, 2: 6, 3: 5, 4: 4}
NISQ_DEPTHS_CACHED = {1: 12, 2: 5, 3: 4, 4: 4}  # reuse the pre-built pickle only
NISQ_DEPTHS_DEEP = {1: 14, 2: 8, 3: 6, 4: 4}  # slower builds, more reductions

# Depths for the exact Clifford sub-pool graphs used by --hybrid.  The
# Clifford sub-pool {RX, RZ @ +/-pi/2, CZ} is far smaller than the full NISQ
# pool, so its graphs can go deeper at a fraction of the build cost.
HYBRID_EXACT_DEPTHS = {1: 12, 2: 8, 3: 6, 4: 5}

# Sampling weights matching the paper's input compositions.
ION_WEIGHTS = None  # uniform over the 30-token pool -> Table 6 input row
NISQ_WEIGHTS = {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}  # -> Table 7 input row

# Paper's "Ours" results (mean over 100 runs, 4 qubits, length 300) and the
# paper's reported baseline means (sums of the per-type means in Tables 6/7).
PAPER_TABLES = {
    "ion_trap": {
        "label": "Table 6",
        "ours": {"RX": (10, 3), "RY": (29, 6), "RZ": (29, 5), "RXX": (43, 8)},
        "baselines": {"qiskit_l1": 196, "qiskit_l2": 204, "qiskit_l3": 204,
                      "bqskit_l2": 154, "bqskit_l3": 129, "bqskit_l4": 126},
    },
    "nisq": {
        "label": "Table 7",
        "ours": {"RX": (45, 6), "RZ": (19, 4), "CZ": (43, 6)},
        "baselines": {"qiskit_l1": 196, "qiskit_l2": 149, "qiskit_l3": 149,
                      "bqskit_l2": 214, "bqskit_l3": 164, "bqskit_l4": 168},
    },
}

GATE_TYPES = {"ion_trap": ["RX", "RY", "RZ", "RXX"], "nisq": ["RX", "RZ", "CZ"]}
BASELINE_METHODS = ["qiskit_l1", "qiskit_l2", "qiskit_l3", "bqskit_l2", "bqskit_l3", "bqskit_l4"]


def _paper_total(gateset: str) -> float:
    return sum(mean for mean, _ in PAPER_TABLES[gateset]["ours"].values())


def _depths_for(gateset: str, depth_mode: str) -> dict[int, int]:
    if gateset == "ion_trap":
        return ION_DEPTHS
    if depth_mode == "cached":
        return NISQ_DEPTHS_CACHED
    if depth_mode == "deep":
        return NISQ_DEPTHS_DEEP
    return NISQ_DEPTHS_DEFAULT


def methods_for(gateset: str, with_numeric: bool) -> list[str]:
    if gateset == "ion_trap":
        methods = ["exact_len", "exact_cost"]
        if with_numeric:
            methods.append("numeric_len")
        return methods
    return ["numeric_len", "numeric_cost"]


# --------------------------------------------------------------------------- #
# baseline compilers (qiskit / BQSKit) with lazy imports
# --------------------------------------------------------------------------- #

_QISKIT = {"state": "unknown"}


def _qiskit_transpile(gates, num_qubits: int, gateset: str, level: int):
    """Transpile a pool circuit with qiskit to the target basis at level L.

    Returns (counts, total, ok) or None if qiskit is unavailable/failed.
    """
    if _QISKIT["state"] == "missing":
        return None
    try:
        if _QISKIT["state"] != "ok":
            from qiskit import QuantumCircuit, transpile
            from qiskit.circuit.library import RXGate, RYGate, RZGate, CZGate, RXXGate
            from qiskit.quantum_info import Operator

            _QISKIT.update({"state": "ok", "QuantumCircuit": QuantumCircuit,
                            "transpile": transpile, "RXGate": RXGate, "RYGate": RYGate,
                            "RZGate": RZGate, "CZGate": CZGate, "RXXGate": RXXGate,
                            "Operator": Operator})
        Q = _QISKIT
        qc = Q["QuantumCircuit"](num_qubits)
        # qiskit numbers qubits little-endian (qubit 0 = LSB); our engine is
        # big-endian, so reverse the wire labels for the unitary comparison.
        rev = lambda qubits: [num_qubits - 1 - q for q in qubits]
        for g in gates:
            if g.name == "RX":
                qc.append(Q["RXGate"](g.theta), rev(g.qubits))
            elif g.name == "RY":
                qc.append(Q["RYGate"](g.theta), rev(g.qubits))
            elif g.name == "RZ":
                qc.append(Q["RZGate"](g.theta), rev(g.qubits))
            elif g.name == "RXX":
                qc.append(Q["RXXGate"](g.theta), rev(g.qubits))
            elif g.name == "CZ":
                qc.append(Q["CZGate"](), rev(g.qubits))
        basis = ["rx", "ry", "rz", "rxx"] if gateset == "ion_trap" else ["rx", "rz", "cz"]
        t = Q["transpile"](qc, basis_gates=basis, optimization_level=level)
        ops = t.count_ops()
        total = sum(n for name, n in ops.items() if name not in ("global_phase", "id", "delay"))
        counts = {
            "RX": ops.get("rx", 0), "RY": ops.get("ry", 0), "RZ": ops.get("rz", 0),
            "RXX": ops.get("rxx", 0), "CZ": ops.get("cz", 0),
        }
        u1 = Q["Operator"](t).data
        ok = equivalent_up_to_global_phase(circuit_unitary(num_qubits, gates), u1, atol=1e-5)
        return counts, total, ok
    except Exception:
        _QISKIT["state"] = "missing"
        return None


_BQSKIT = {"state": "unknown"}


def _bqskit_compile(gates, num_qubits: int, gateset: str, level: int):
    """Compile a pool circuit with BQSKit (optimization levels 2..4).

    Returns (counts, total, ok) or None if BQSKit is unavailable/failed.
    """
    if _BQSKIT["state"] == "missing":
        return None
    try:
        if _BQSKIT["state"] != "ok":
            import numpy as np
            from bqskit import Circuit, MachineModel, compile
            from bqskit.compiler.gateset import GateSet
            from bqskit.ir.gates import RXGate, RYGate, RZGate, CZGate, RXXGate

            _BQSKIT.update({"state": "ok", "np": np, "Circuit": Circuit,
                            "MachineModel": MachineModel, "compile": compile,
                            "GateSet": GateSet, "RXGate": RXGate, "RYGate": RYGate,
                            "RZGate": RZGate, "CZGate": CZGate, "RXXGate": RXXGate})
        B = _BQSKIT
        c = B["Circuit"](num_qubits)
        for g in gates:
            if g.name == "RX":
                c.append_gate(B["RXGate"]().with_all_frozen_params([g.theta]), [g.qubits[0]])
            elif g.name == "RY":
                c.append_gate(B["RYGate"]().with_all_frozen_params([g.theta]), [g.qubits[0]])
            elif g.name == "RZ":
                c.append_gate(B["RZGate"]().with_all_frozen_params([g.theta]), [g.qubits[0]])
            elif g.name == "RXX":
                c.append_gate(B["RXXGate"]().with_all_frozen_params([g.theta]), [g.qubits[0], g.qubits[1]])
            elif g.name == "CZ":
                c.append_gate(B["CZGate"](), [g.qubits[0], g.qubits[1]])
        basis = [B["RXGate"](), B["RYGate"](), B["RZGate"](), B["RXXGate"]()] if gateset == "ion_trap" \
            else [B["RXGate"](), B["RZGate"](), B["CZGate"]()]
        model = B["MachineModel"](num_qubits, gate_set=B["GateSet"](list(basis)))
        out = B["compile"](c, model=model, optimization_level=level)
        counts = {}
        for gate, n in dict(out.gate_counts).items():
            name = type(gate).__name__.replace("Gate", "")
            counts[name] = counts.get(name, 0) + int(n)
        total = out.num_operations
        u1 = B["np"].asarray(out.get_unitary())
        ok = equivalent_up_to_global_phase(circuit_unitary(num_qubits, gates), u1, atol=1e-5)
        return counts, total, ok
    except Exception:
        _BQSKIT["state"] = "missing"
        return None


# --------------------------------------------------------------------------- #
# per-process database cache (avoids pickling multi-100 MB DBs to workers)
# --------------------------------------------------------------------------- #

_DB_CACHE: dict = {}


def _load_db(kind: str, gateset: str, depths: dict[int, int], backend: str = "ram"):
    key = (kind, gateset, tuple(sorted(depths.items())), backend)
    cached = _DB_CACHE.get(key)
    if cached is not None:
        return cached
    if kind == "exact":
        db = load_or_build_exact(gateset, depths)
    else:
        db = load_or_build_database(gateset, depths, backend=backend)
    _DB_CACHE[key] = db
    return db


def _count_twq(gates) -> int:
    return sum(1 for g in gates if len(g.qubits) == 2)


def _worker(args):
    (gateset, method, num_qubits, length, budget_s, seed, weights, rz_pass, depths,
     max_block_len, restarts, backend, rf_gate, hybrid) = args
    gates, _ = random_circuit(num_qubits, length, gateset, seed=seed, weights=weights)
    if method.startswith("qiskit_"):
        res = _qiskit_transpile(gates, num_qubits, gateset, int(method[-1]))
        if res is None:
            return {"seed": seed, "method": method, "start": len(gates), "end": -1,
                    "counts": {}, "twq": -1, "secs": 0.0, "ok": False, "verifier": "unavailable"}
        counts, total, ok = res
        return {"seed": seed, "method": method, "start": len(gates), "end": total, "counts": counts,
                "twq": counts.get("RXX", 0) + counts.get("CZ", 0), "secs": 0.0, "ok": ok,
                "verifier": "qiskit-1e-5"}
    if method.startswith("bqskit_"):
        res = _bqskit_compile(gates, num_qubits, gateset, int(method[-1]))
        if res is None:
            return {"seed": seed, "method": method, "start": len(gates), "end": -1,
                    "counts": {}, "twq": -1, "secs": 0.0, "ok": False, "verifier": "unavailable"}
        counts, total, ok = res
        return {"seed": seed, "method": method, "start": len(gates), "end": total, "counts": counts,
                "twq": counts.get("RXX", 0) + counts.get("CZ", 0), "secs": 0.0, "ok": ok,
                "verifier": "bqskit-1e-5"}
    cost_aware = method in ("exact_cost", "numeric_cost")
    restarts = max(1, restarts)
    if method.startswith("exact_"):
        db = _load_db("exact", gateset, depths)
        t0 = time.time()
        best_r = None
        best_key = None
        for i in range(restarts):
            r, stats = reduce_circuit_exact(
                gates, num_qubits, db, budget_s, seed + i * 10000,
                cost_aware=cost_aware, max_block_len=max_block_len,
            )
            key = (_count_twq(r), len(r)) if cost_aware else (len(r),)
            if best_r is None or key < best_key:
                best_r, best_key, best_stats = r, key, stats
        r, stats = best_r, best_stats
        ok = verify_exact(gates, r, num_qubits)
        verifier = "exact"
        secs = stats.runtime_sec  # excludes DB load, like the numeric path below
    else:
        # SQLite-backed databases are re-opened in each worker (forked
        # processes must not share a parent's SQLite connection).
        if backend == "sqlite":
            db = load_or_build_database(gateset, depths, backend="sqlite")
        else:
            db = _load_db("numeric", gateset, depths, "ram")
        if hybrid:
            if gateset != "nisq":
                raise ValueError("--hybrid is implemented for the NISQ pool only")
            exact_clifford = _load_db("exact", "nisq_clifford", HYBRID_EXACT_DEPTHS)
            numeric_gate = RfGate() if rf_gate else None
            db = HybridDatabase(db, exact_clifford, numeric_gate=numeric_gate)
        elif rf_gate:
            db = RfGatedDatabase(db, RfGate())
        t0 = time.time()
        best_r = None
        best_key = None
        for i in range(restarts):
            r, _passes, _reds = reduce_circuit(
                gates,
                num_qubits,
                db,
                budget_s,
                seed + i * 10000,
                rz_pass=rz_pass,
                cost_aware=cost_aware,
                max_block_len=max_block_len,
            )
            key = (_count_twq(r), len(r)) if cost_aware else (len(r),)
            if best_r is None or key < best_key:
                best_r, best_key = r, key
        r = best_r
        ok = equivalent_up_to_global_phase(
            circuit_unitary(num_qubits, gates), circuit_unitary(num_qubits, r), atol=1e-5
        )
        verifier = "numeric-1e-5"
        secs = time.time() - t0
    return {
        "seed": seed,
        "method": method,
        "start": len(gates),
        "end": len(r),
        "counts": count_gates(r),
        "twq": sum(1 for g in r if len(g.qubits) == 2),
        "secs": secs,
        "ok": ok,
        "verifier": verifier,
    }


def _probe_compiler_availability() -> tuple[bool, bool]:
    """Check qiskit/BQSKit importability in a fresh interpreter."""
    import subprocess
    import sys
    import textwrap

    code = textwrap.dedent(
        """
        out = []
        for mod in ("qiskit", "bqskit"):
            try:
                __import__(mod)
                out.append("1")
            except Exception:
                out.append("0")
        print("".join(out))
        """
    )
    try:
        r = subprocess.run([sys.executable, "-c", code], capture_output=True,
                           text=True, timeout=120)
        flags = (r.stdout or "").strip()
        return flags[:1] == "1", flags[1:2] == "1"
    except Exception:
        return False, False


def _run_tasks(tasks: list, workers: int) -> list[dict]:
    """Run tasks, falling back to less exotic multiprocessing, then sequential."""
    if workers <= 1 or len(tasks) <= 1:
        return [_worker(t) for t in tasks]
    for ctx_name in ("fork", "spawn"):
        try:
            ctx = mp.get_context(ctx_name)
            pool = ctx.Pool(processes=workers)
        except Exception:
            continue
        try:
            results = pool.map(_worker, tasks)
            pool.close()
            pool.join()
            return results
        except Exception as exc:  # child crash -> fall through to next strategy
            print(f"  [warn] multiprocessing ({ctx_name}) failed: {exc}; retrying...")
            try:
                pool.terminate()
            except Exception:
                pass
    print("  [warn] falling back to sequential execution")
    return [_worker(t) for t in tasks]


# --------------------------------------------------------------------------- #
# reporting
# --------------------------------------------------------------------------- #

def _mean_std(values) -> tuple[float, float]:
    if not values:
        return float("nan"), float("nan")
    return statistics.mean(values), statistics.pstdev(values)


def _stats_for(results: list[dict], method: str, gate_types: list[str]) -> dict:
    rows = [r for r in results if r["method"] == method and r["end"] >= 0]
    if not rows:
        return {}
    ends = [r["end"] for r in rows]
    twqs = [r["twq"] for r in rows]
    secs = [r["secs"] for r in rows]
    counts = {name: [r["counts"].get(name, 0) for r in rows] for name in gate_types}
    return {
        "n": len(rows),
        "end_mean": statistics.mean(ends),
        "end_std": statistics.pstdev(ends),
        "end_best": min(ends),
        "twq_mean": statistics.mean(twqs),
        "twq_std": statistics.pstdev(twqs),
        "counts": {name: _mean_std(vals) for name, vals in counts.items()},
        "secs_mean": statistics.mean(secs),
        "ok_rate": sum(r["ok"] for r in rows) / len(rows),
    }


def _build_report(gateset: str, results: list[dict], stats: dict, meta: dict) -> str:
    paper = PAPER_TABLES[gateset]
    types = GATE_TYPES[gateset]
    lines = [
        "# Comparison benchmark report",
        "",
        f"- Gate set: `{gateset}` (paper {paper['label']}, 4 qubits, length 300)",
        f"- Circuits per method: {meta['num_circuits']}, per-circuit budget: {meta['budget']}s",
        f"- Generated: {meta['date']}",
        f"- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)",
        "",
        "Paper reference ('Ours', mean +/- std over 100 runs):",
        "",
        "| method | " + " | ".join(types) + " | total |",
        "|---" * (len(types) + 2) + "|",
        "| **paper** | "
        + " | ".join(f"{m:.0f} (+- {s:.0f})" for m, s in paper["ours"].values())
        + f" | {_paper_total(gateset):.0f} |",
        "",
        "All methods on identical circuits (mean +/- std):",
        "",
        "| method | " + " | ".join(types) + " | total | vs paper Ours | two-qubit | time (s) |",
        "|---" * (len(types) + 5) + "|",
    ]
    for method in meta["method_order"]:
        st = stats.get(method)
        if not st:
            continue
        cells = " | ".join(f"{st['counts'][n][0]:.1f} (+- {st['counts'][n][1]:.1f})" for n in types)
        if method in paper["baselines"]:
            note = f"base (paper {paper['baselines'][method]:.0f})"
        elif method in ("exact_len", "exact_cost", "numeric_len", "numeric_cost"):
            margin = st["end_mean"] - _paper_total(gateset)
            note = "WIN" if margin < 0 else ("TIE" if margin == 0 else "LOSE")
            note += f" ({margin:+.1f})"
        else:
            note = "--"
        lines.append(
            f"| **{method}** | {cells} | {st['end_mean']:.1f} (+- {st['end_std']:.1f}) | "
            f"{note} | {st['twq_mean']:.1f} | {st['secs_mean']:.1f} |"
        )

    lines += ["", "Baseline fidelity check (our means vs paper's reported baseline means):", "",
              "| baseline | our total | paper total | delta |", "|---|---:|---:|---:|"]
    for method in paper["baselines"]:
        st = stats.get(method)
        if not st:
            lines.append(f"| {method} | (skipped) | {paper['baselines'][method]:.0f} | -- |")
            continue
        delta = st["end_mean"] - paper["baselines"][method]
        lines.append(f"| {method} | {st['end_mean']:.1f} | {paper['baselines'][method]:.0f} | {delta:+.1f} |")

    lines += ["", "Verdict notes:", ""]
    for method in ("exact_len", "exact_cost", "numeric_len", "numeric_cost"):
        st = stats.get(method)
        if not st:
            continue
        paper_twq = paper["ours"][types[-1]][0]
        margin = st["end_mean"] - _paper_total(gateset)
        verdict = "WIN" if margin < 0 else ("TIE" if margin == 0 else "LOSE")
        twq_verdict = "WIN" if st["twq_mean"] < paper_twq else "TIE" if st["twq_mean"] == paper_twq else "LOSE"
        lines.append(
            f"- `{method}`: total {verdict} vs paper ({st['end_mean']:.1f} vs "
            f"{_paper_total(gateset):.0f}, {st['end_mean'] - _paper_total(gateset):+.1f}); "
            f"two-qubit {twq_verdict} ({st['twq_mean']:.1f} vs paper {paper_twq}); "
            f"equivalence pass rate {st['ok_rate']:.3f}; best {st['end_best']}."
        )
    if gateset == "nisq":
        lines += [
            "",
            "Note: NISQ inputs are CZ-weighted (weights RX:1, RZ:1, CZ:2) to match the paper's",
            "Table 7 input composition (RX~108, RZ~109, CZ~82). The RZ-across-CZ pass is enabled",
            "(iterated to a fixpoint).",
        ]
        if meta.get("hybrid"):
            lines += [
                "Exact/numeric hybrid enabled: Clifford-only windows (RX/RZ at +/-pi/2, CZ) are",
                "reduced by the exact symplectic engine at deep graph depths; only non-Clifford",
                "(pi/4) windows hit the numeric database.",
            ]
        if meta.get("rf_gate"):
            lines += [
                "V3-style RF-gated lookup enabled: an online classifier skips lookups on blocks",
                "predicted irreducible, freeing budget for blocks that actually reduce.",
            ]
    lines += [
        "",
        "Timing caveat: the \"time (s)\" column is the per-circuit budget cap -- each reducer",
        "loops until its budget is exhausted. It is a cutoff, not a convergence time, and is",
        "not directly comparable to the paper's Table 2 (a different task: reducing 100-gate",
        "circuits to ~50, ~38 s for their best variant).",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--gateset", type=str, required=True, choices=["ion_trap", "nisq"])
    parser.add_argument("--num-circuits", type=int, default=None, help="circuits per method (default 100, or 8 with --quick)")
    parser.add_argument("--num-qubits", type=int, default=4)
    parser.add_argument("--length", type=int, default=300)
    parser.add_argument("--budget", type=float, default=None,
                        help="per-circuit time budget (s, default 30 ion_trap / 60 nisq, or 10 with --quick)")
    parser.add_argument("--seed-base", type=int, default=0)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--rz-pass", action="store_true", help="RZ-across-CZ pass (always on for nisq)")
    parser.add_argument("--numeric", action="store_true", help="also run numeric reducer for ion_trap")
    parser.add_argument("--depths", type=str, default="", help='override DB depths, e.g. "1:12,2:8,3:5,4:4"')
    parser.add_argument("--deep", action="store_true", help="use deeper NISQ graphs (slower build, more reductions)")
    parser.add_argument("--backend", type=str, default=None, choices=["ram", "sqlite"],
                        help="lookup DB storage backend (default: ram; --deep implies sqlite, which builds the\n                        deep graphs within a laptop's memory)")
    parser.add_argument("--rf-gate", action="store_true",
                        help="V3-style RandomForest-gated lookup: learn which blocks reduce and skip\n                        useless DB lookups (needs optional extra [ml])")
    parser.add_argument("--hybrid", action="store_true",
                        help="NISQ: solve Clifford-only windows (RX/RZ at +/-pi/2, CZ) with the exact\n                        symplectic engine; only non-Clifford (pi/4) windows hit the numeric DB")
    parser.add_argument("--block-len", type=int, default=None, help="max sweep window length (default 8, or 10 for nisq/--deep)")
    parser.add_argument("--no-baselines", action="store_true", help="skip qiskit/BQSKit baseline compilers")
    parser.add_argument("--baselines-only", action="store_true", help="run only the qiskit/BQSKit baselines (no reducers, no DB build)")
    parser.add_argument("--bqskit-circuits", type=int, default=None,
                        help="max circuits for BQSKit baselines (slow; default = num-circuits)")
    parser.add_argument("--bqskit-levels", type=str, default="2,3,4",
                        help="BQSKit optimization levels to run, comma-separated (default 2,3,4)")
    parser.add_argument("--restarts", type=int, default=1,
                        help="best-of-K independent reductions per circuit (more compute, better results)")
    parser.add_argument("--outdir", type=str, default="results/comparison")
    parser.add_argument("--quick", action="store_true", help="8 circuits, 10 s budget (smoke test)")
    args = parser.parse_args()

    if args.quick:
        args.num_circuits = args.num_circuits or 8
        args.budget = args.budget or 10.0
    else:
        args.num_circuits = args.num_circuits or 100
        # NISQ's search space is much larger than ion trap's, so it gets a
        # larger honestly-reported default budget (the paper admits their
        # approach "takes much longer (up to several minutes)" per circuit).
        args.budget = args.budget or (60.0 if args.gateset == "nisq" else 30.0)

    weights = ION_WEIGHTS if args.gateset == "ion_trap" else NISQ_WEIGHTS
    rz_pass = args.rz_pass or args.gateset == "nisq"

    backend = args.backend or ("sqlite" if args.deep else "ram")
    if args.hybrid and args.gateset != "nisq":
        parser.error("--hybrid is implemented for the NISQ pool only")

    depths = _depths_for(args.gateset, "deep" if args.deep else "default")
    if args.depths:
        depths = {}
        for item in args.depths.split(","):
            w, d = item.split(":")
            depths[int(w)] = int(d)

    max_block_len = args.block_len
    if max_block_len is None:
        max_block_len = 10 if (args.deep or args.gateset == "nisq") else 8

    if args.gateset == "nisq" and not args.depths:
        mode = "deep" if args.deep else "default"
        print(
            f"[note] NISQ {mode} depths {depths}: the 1/2/3-wire graphs will be rebuilt once and "
            f"cached (.cache/, backend {backend}). Use --depths 1:12,2:5,3:4,4:4 to reuse the "
            "pre-built database. (--deep implies the sqlite backend so the deep graphs build "
            "within a laptop's memory; they build more slowly but reduce more.)",
            flush=True,
        )

    bqskit_levels = {int(x) for x in args.bqskit_levels.split(",") if x.strip()}
    if args.baselines_only:
        ours_methods = []
        baselines = [] if args.no_baselines else list(BASELINE_METHODS)
    else:
        ours_methods = methods_for(args.gateset, args.numeric)
        baselines = [] if args.no_baselines else list(BASELINE_METHODS)
    baselines = [m for m in baselines if not m.startswith("bqskit_") or int(m[-1]) in bqskit_levels]
    methods = ours_methods + baselines
    seeds = [args.seed_base + s for s in range(args.num_circuits)]

    # --- pre-build / load databases in the parent (single build, shared COW) ---
    for kind in ("exact", "numeric"):
        if any(m.startswith(kind) for m in ours_methods):
            t0 = time.time()
            _load_db(kind, args.gateset, depths, backend)
            print(f"[{args.gateset}] {kind} database ready ({time.time() - t0:.1f}s, "
                  f"depths {depths}, backend {backend})", flush=True)
    if args.hybrid:
        t0 = time.time()
        _load_db("exact", "nisq_clifford", HYBRID_EXACT_DEPTHS)
        print(f"[nisq] hybrid exact (Clifford sub-pool) database ready ({time.time() - t0:.1f}s)",
              flush=True)

    # --- availability probes for baselines (fresh interpreter: importing
    # qiskit/BQSKit in the parent would poison the fork used by the worker
    # pool -- the workers import them lazily instead) ---
    if baselines:
        qiskit_ok, bqskit_ok = _probe_compiler_availability()
        for method in baselines:
            avail = qiskit_ok if method.startswith("qiskit_") else bqskit_ok
            print(f"[{args.gateset}] {method}: {'ok' if avail else 'SKIPPED (compiler unavailable)'}",
                  flush=True)

    bqskit_max = args.bqskit_circuits if args.bqskit_circuits is not None else args.num_circuits
    tasks = []
    for m in methods:
        for s in seeds:
            if m.startswith("bqskit_") and s - args.seed_base >= bqskit_max:
                continue
            tasks.append((args.gateset, m, args.num_qubits, args.length, args.budget, s,
                          weights, rz_pass, depths, max_block_len, max(1, args.restarts),
                          backend, args.rf_gate, args.hybrid))

    print(
        f"[{args.gateset}] {args.num_circuits} circuits x {args.length} gates (q{args.num_qubits}), "
        f"budget {args.budget}s, methods {methods}, depths {depths}",
        flush=True,
    )
    t0 = time.time()
    workers = args.workers or min(len(tasks), os.cpu_count() or 4)
    results = _run_tasks(tasks, workers)
    wall = time.time() - t0

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / f"comparison_{args.gateset}.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["method", "seed", "start", "end", "runtime_s", "ok", "verifier", "twq",
                         "RX", "RY", "RZ", "RXX", "CZ"])
        for r in sorted(results, key=lambda r: (r["method"], r["seed"])):
            c = r["counts"]
            writer.writerow([r["method"], r["seed"], r["start"], r["end"], round(r["secs"], 3),
                             r["ok"], r["verifier"], r["twq"],
                             c.get("RX", 0), c.get("RY", 0), c.get("RZ", 0), c.get("RXX", 0), c.get("CZ", 0)])

    stats = {m: _stats_for(results, m, GATE_TYPES[args.gateset]) for m in methods}
    meta = {
        "gateset": args.gateset,
        "num_circuits": args.num_circuits,
        "num_qubits": args.num_qubits,
        "length": args.length,
        "budget": args.budget,
        "depths": depths,
        "max_block_len": max_block_len,
        "rz_pass": rz_pass,
        "weights": weights,
        "restarts": max(1, args.restarts),
        "backend": backend,
        "rf_gate": args.rf_gate,
        "hybrid": args.hybrid,
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "wall_sec": round(wall, 1),
        "with_numeric": args.numeric,
        "method_order": methods,
    }
    report = _build_report(args.gateset, results, stats, meta)
    (outdir / f"comparison_{args.gateset}_report.md").write_text(report, encoding="utf-8")
    (outdir / f"comparison_{args.gateset}_report.json").write_text(
        json.dumps({"meta": meta, "methods": {m: s for m, s in stats.items()}}, indent=2),
        encoding="utf-8",
    )

    # console summary ------------------------------------------------------- #
    print(f"[{args.gateset}] wall {wall:.1f}s across {workers} workers", flush=True)
    for method in methods:
        st = stats[method]
        if not st:
            continue
        types = GATE_TYPES[args.gateset]
        per_type = " ".join(f"{st['counts'][n][0]:.1f}" for n in types)
        if method in PAPER_TABLES[args.gateset]["baselines"]:
            note = f"paper base {PAPER_TABLES[args.gateset]['baselines'][method]:.0f}"
        else:
            margin = st["end_mean"] - _paper_total(args.gateset)
            verdict = "WIN" if margin < 0 else ("TIE" if margin == 0 else "LOSE")
            note = f"vs paper {verdict} ({margin:+.1f})"
        print(
            f"  {method:<12} end {st['end_mean']:6.1f} (+- {st['end_std']:4.1f})  best {st['end_best']:3d}  "
            f"twq {st['twq_mean']:5.1f}  [{per_type}]  {note}  "
            f"ok {st['ok_rate']:.3f}  time {st['secs_mean']:5.1f}s",
            flush=True,
        )
    print(f"  saved CSV: {csv_path}")
    print(f"  saved report: {outdir / f'comparison_{args.gateset}_report.md'}")


if __name__ == "__main__":
    main()
