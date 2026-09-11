from __future__ import annotations

from qcr_repro.circuits import random_circuit
from qcr_repro.database import load_or_build_database
from qcr_repro.gates import circuit_unitary
from qcr_repro.reducer import _sweep_reduce, cluster_single_qubit, reduce_single_wire_runs
from qcr_repro.symmetry import SymmetricDatabase
from qcr_repro.unitary import equivalent_up_to_global_phase

ION_DEPTHS = {1: 12, 2: 10, 3: 7, 4: 5}

db = load_or_build_database("ion_trap", ION_DEPTHS)
sym = SymmetricDatabase(load_or_build_database("ion_trap", ION_DEPTHS))

gates, _ = random_circuit(4, 30, "ion_trap", seed=1)
u0 = circuit_unitary(4, gates)

base = list(gates)
one_wire_base = db.graphs.get(1)
cluster_single_qubit(base, 4)
print("after cluster:", len(base), equivalent_up_to_global_phase(u0, circuit_unitary(4, base), atol=1e-8))
reduce_single_wire_runs(base, one_wire_base)
print("after single-wire collapse:", len(base), equivalent_up_to_global_phase(u0, circuit_unitary(4, base), atol=1e-8))
n_replaced = _sweep_reduce(base, 4, db, 8)
print("after sweep (baseline db):", len(base), "replaced", n_replaced,
      equivalent_up_to_global_phase(u0, circuit_unitary(4, base), atol=1e-8))

sym_gates = list(gates)
one_wire_sym = sym.graphs.get(1)
cluster_single_qubit(sym_gates, 4)
print("after cluster (sym):", len(sym_gates), equivalent_up_to_global_phase(u0, circuit_unitary(4, sym_gates), atol=1e-8))
reduce_single_wire_runs(sym_gates, one_wire_sym)
print("after single-wire collapse (sym):", len(sym_gates), equivalent_up_to_global_phase(u0, circuit_unitary(4, sym_gates), atol=1e-8))

# Step the sweep one window at a time to localize the first bad replacement.
working = list(sym_gates)
u_before_sweep = circuit_unitary(4, working)
pos = 0
max_block_len = 8
step = 0
while pos < len(working):
    hi = min(max_block_len, len(working) - pos)
    replaced = False
    for length in range(hi, 1, -1):
        block = working[pos:pos + length]
        u_block = circuit_unitary(4, block)
        candidate = sym.try_reduce(block)
        if candidate is not None and len(candidate) < length:
            u_cand = circuit_unitary(4, candidate)
            ok = equivalent_up_to_global_phase(u_block, u_cand, atol=1e-6)
            step += 1
            if not ok:
                print(f"STEP {step}: pos={pos} length={length} -> {len(candidate)}  BLOCK MISMATCH  ok={ok}")
                print("  block:", [(g.name, g.qubits, g.theta) for g in block])
                print("  candidate:", [(g.name, g.qubits, g.theta) for g in candidate])
                raise SystemExit(1)
            working[pos:pos + length] = candidate
            replaced = True
            break
    if replaced:
        pos = max(0, pos - 1)
    else:
        pos += 1

print("manual step-through sweep ok, final len", len(working),
      equivalent_up_to_global_phase(u_before_sweep, circuit_unitary(4, working), atol=1e-6))
