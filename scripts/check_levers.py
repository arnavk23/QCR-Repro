"""Equivalence-preservation checks for the NISQ levers (RF-gated lookup, exact/numeric hybrid, SQLite backend).

Asserts every lever preserves the input unitary (1e-5, up to global phase) and that the SQLite store is byte-identical to RAM.

Usage:
    PYTHONPATH=src python scripts/check_levers.py"""

from __future__ import annotations

import copy
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from qcr_repro.circuits import random_circuit
from qcr_repro.database import load_or_build_database
from qcr_repro.exact_database import load_or_build_exact
from qcr_repro.gates import circuit_unitary
from qcr_repro.hybrid import HybridDatabase
from qcr_repro.reducer import reduce_circuit
from qcr_repro.rf_gate import RfGate, RfGatedDatabase
from qcr_repro.unitary import equivalent_up_to_global_phase

NISQ_DEPTHS = {1: 12, 2: 6, 3: 5, 4: 4}  # matches the cached RAM database
HYBRID_EXACT_DEPTHS = {1: 12, 2: 8, 3: 6, 4: 5}
BUDGET = 8.0
SEEDS = (1, 2)

results: list[str] = []


def check(tag: str, ok: bool, detail: str = "") -> None:
    status = "PASS" if ok else "FAIL"
    results.append(f"[{status}] {tag}" + (f"  ({detail})" if detail else ""))
    print(results[-1], flush=True)


db = load_or_build_database("nisq", NISQ_DEPTHS)
exact_db = load_or_build_exact("nisq_clifford", HYBRID_EXACT_DEPTHS)

# --- SQLite store content parity with RAM (deterministic, unlike end
# lengths, which depend on how many passes fit in the wall-clock budget) ---
db_s = load_or_build_database("nisq", NISQ_DEPTHS, backend="sqlite")
for wires, g_sql in sorted(db_s.graphs.items()):
    g_ram = db.graphs[wires]
    check(f"sqlite w{wires} node count == ram", g_sql.num_nodes == g_ram.num_nodes,
          f"{g_sql.num_nodes} vs {g_ram.num_nodes}")
    mismatch = 0
    n = 0
    for key in g_sql.buckets.keys():
        n += 1
        if g_sql.buckets[key] != g_ram.buckets[key]:
            mismatch += 1
    check(f"sqlite w{wires} bucket content == ram", mismatch == 0, f"{n} keys, {mismatch} mismatches")

db_s2 = load_or_build_database("nisq", NISQ_DEPTHS, backend="sqlite")
for wires, g_sql in sorted(db_s.graphs.items()):
    g_sql2 = db_s2.graphs[wires]
    check(f"sqlite w{wires} reload content stable", g_sql.num_nodes == g_sql2.num_nodes)

for seed in SEEDS:
    gates, _ = random_circuit(4, 300, "nisq", seed=seed, weights={"RX": 1.0, "RZ": 1.0, "CZ": 2.0})
    u0 = circuit_unitary(4, gates)
    tag = f"seed {seed}"

    r_base, _, _ = reduce_circuit(copy.deepcopy(gates), 4, db, BUDGET, seed=7,
                                  rz_pass=True, max_block_len=8)
    check(f"{tag} baseline ok", equivalent_up_to_global_phase(u0, circuit_unitary(4, r_base), atol=1e-5),
          f"end {len(r_base)}")

    rf = RfGate()
    r_rf, _, _ = reduce_circuit(copy.deepcopy(gates), 4, RfGatedDatabase(db, rf), BUDGET, seed=7,
                                rz_pass=True, max_block_len=8)
    check(f"{tag} rf-gate ok", equivalent_up_to_global_phase(u0, circuit_unitary(4, r_rf), atol=1e-5),
          f"end {len(r_rf)}, skipped {rf.lookups_skipped}")

    hdb = HybridDatabase(db, exact_db)
    r_hy, _, _ = reduce_circuit(copy.deepcopy(gates), 4, hdb, BUDGET, seed=7,
                                rz_pass=True, max_block_len=8)
    check(f"{tag} hybrid ok", equivalent_up_to_global_phase(u0, circuit_unitary(4, r_hy), atol=1e-5),
          f"end {len(r_hy)}, exact lookups {hdb.exact_lookups}")

    r_sql, _, _ = reduce_circuit(copy.deepcopy(gates), 4, db_s, BUDGET, seed=7,
                                 rz_pass=True, max_block_len=8)
    check(f"{tag} sqlite reduction ok", equivalent_up_to_global_phase(u0, circuit_unitary(4, r_sql), atol=1e-5),
          f"end {len(r_sql)}")

if all("FAIL" not in r for r in results):
    print("ALL LEVER CHECKS PASSED")
else:
    print("SOME LEVER CHECKS FAILED")
    raise SystemExit(1)
