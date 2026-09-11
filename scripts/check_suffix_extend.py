"""SuffixExtendedDatabase is a unitary-preserving, gate-set-agnostic reach
extension: undoing up to max_extra trailing gates from a block's unitary and
re-querying the same graph, to reach past its own max_depth without
rebuilding it deeper. Checks (1) correctness on direct probed windows for
both gate sets, and (2) an A/B of reduce_circuit with vs without it.

Usage:
    PYTHONPATH=src python scripts/check_suffix_extend.py
"""

from __future__ import annotations

import sys
import time

from qcr_repro.circuits import random_circuit
from qcr_repro.database import load_or_build_database
from qcr_repro.gates import circuit_unitary
from qcr_repro.reducer import reduce_circuit
from qcr_repro.suffix_extend import SuffixExtendedDatabase
from qcr_repro.unitary import equivalent_up_to_global_phase

ION_DEPTHS = {1: 12, 2: 10, 3: 7, 4: 5}
NISQ_DEPTHS = {1: 12, 2: 6, 3: 5, 4: 4}

# rz_pass is NISQ-only (RZ-across-CZ gathering assumes RX is the only
# non-diagonal gate, false for ion_trap's RY/RXX) -- matches
# benchmark_comparison.py's own rz_pass = args.rz_pass or gateset == "nisq".
CASES = [
    ("ion_trap", None, ION_DEPTHS, False),
    ("nisq", {"RX": 1.0, "RZ": 1.0, "CZ": 2.0}, NISQ_DEPTHS, True),
]


def check_try_reduce_unitary_preserving() -> int:
    failures = 0
    for gate_set, weights, depths, _rz_pass in CASES:
        db = SuffixExtendedDatabase(load_or_build_database(gate_set, depths), max_extra=1)
        for seed in range(1, 4):
            gates, _ = random_circuit(4, 100, gate_set, seed=seed, weights=weights)
            n = len(gates)
            checked = 0
            for start in range(0, n, 9):
                for length in range(2, min(8, n - start) + 1):
                    block = gates[start : start + length]
                    u0 = circuit_unitary(4, block)
                    cand = db.try_reduce(block)
                    checked += 1
                    if cand is None:
                        continue
                    u1 = circuit_unitary(4, cand)
                    if not equivalent_up_to_global_phase(u0, u1, atol=1e-8):
                        print(f"  FAIL {gate_set} seed{seed} start{start} len{length}: "
                              f"try_reduce changed the unitary (cand len {len(cand)})")
                        failures += 1
            print(f"  ok   {gate_set} seed{seed}: {checked} windows probed, "
                  f"direct_hits={db.direct_hits} extended_hits={db.extended_hits}")
    return failures


def check_reduce_circuit_ab(num_seeds: int = 8, budget_s: float = 6.0) -> int:
    failures = 0
    for gate_set, weights, depths, rz_pass in CASES:
        base_db = load_or_build_database(gate_set, depths)
        ext_db = SuffixExtendedDatabase(load_or_build_database(gate_set, depths), max_extra=1)
        total_base = 0
        total_ext = 0
        t_base = 0.0
        t_ext = 0.0
        print(f"-- {gate_set} --")
        for seed in range(1, num_seeds + 1):
            gates, _ = random_circuit(4, 300, gate_set, seed=seed, weights=weights)
            u0 = circuit_unitary(4, gates)

            t0 = time.time()
            base, pb, _ = reduce_circuit(
                list(gates), 4, base_db, budget_s=budget_s, seed=seed, rz_pass=rz_pass
            )
            t_base += time.time() - t0

            t0 = time.time()
            ext, pe, _ = reduce_circuit(
                list(gates), 4, ext_db, budget_s=budget_s, seed=seed, rz_pass=rz_pass
            )
            t_ext += time.time() - t0

            ok_base = equivalent_up_to_global_phase(u0, circuit_unitary(4, base), atol=1e-5)
            ok_ext = equivalent_up_to_global_phase(u0, circuit_unitary(4, ext), atol=1e-5)
            status = "ok  " if (ok_base and ok_ext) else "FAIL"
            total_base += len(base)
            total_ext += len(ext)
            print(f"  {status} seed{seed}: baseline {len(base)} (passes {pb}, ok={ok_base}) vs "
                  f"suffix_ext {len(ext)} (passes {pe}, ok={ok_ext})")
            if not (ok_base and ok_ext):
                failures += 1
        print(f"  mean baseline {total_base / num_seeds:.1f}  mean suffix_ext {total_ext / num_seeds:.1f}  "
              f"delta {total_ext / num_seeds - total_base / num_seeds:+.1f}  "
              f"direct_hits={ext_db.direct_hits} extended_hits={ext_db.extended_hits}  "
              f"wall-clock base {t_base:.1f}s ext {t_ext:.1f}s")
    return failures


def main() -> None:
    print("== 1) try_reduce is unitary-preserving on probed windows (ion_trap + nisq) ==")
    failures = check_try_reduce_unitary_preserving()
    print("== 2) reduce_circuit A/B (baseline vs SuffixExtendedDatabase, same seed/budget) ==")
    failures += check_reduce_circuit_ab()
    print(f"\n{'ALL CHECKS PASSED' if failures == 0 else f'{failures} CHECK(S) FAILED'}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
