# DAG block-compaction comparison report

- Gate set: `ion_trap` (paper Table 6, 4 qubits, length 300)
- Circuits per method: 20 (seeds 1-20), per-circuit budget: 30.0s
- Verifier: numeric 1e-5 (input vs output unitary up to global phase)
- Paper 'Ours' (mean over 100 runs): total 111, two-qubit 43

| method | input | after prepass | final total (mean +/- std) | two-qubit | vs paper | t-stat | p | runtime (s) | exact |
|---|---|---|---|---|---|---|---|---|---|
| baseline | 300 | 300 | 68.5 (+/- 15.4) | 28.8 | WIN (-38.2%) | -12.3 | 0 | 30.2 | OK |
| **prepass+dag_compact** | 300 | 284 | 66.1 (+/- 14.9) | 28.2 | WIN (-40.5%) | -13.5 | 0 | 30.1 | OK |

`prepass+dag_compact` vs `baseline` at equal 30.0s budget: 66.1 vs 68.5 gates (-3.6%).

## What dag_compact does

`src/dag.py` deterministically reorders the circuit so every <=3-wire block becomes contiguous, using the true per-wire dependency DAG rather than physical adjacency in the gate list -- generalizing Qiskit-style 2-qubit block collection to k<=3-wire blocks. Any two gates that swap position under this reordering act on disjoint wires, so it is a valid topological order and preserves the circuit's unitary exactly (scripts/check_dag_compact.py). It replaces reliance on the stochastic transport_shuffle/shuffle_commuting_pairs passes for exposing reducible windows: the previous pipeline could only find a window if random shuffling happened to bring its gates physically adjacent; dag_compact finds every such window in one deterministic pass.