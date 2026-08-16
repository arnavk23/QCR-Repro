# Comparison benchmark report

- Gate set: `nisq` (paper Table 7, 4 qubits, length 300)
- Circuits per method: 8, per-circuit budget: 20.0s
- Generated: 2026-08-17 02:25:54
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RZ | CZ | total |
|---|---|---|---|---|
| **paper** | 45 (+- 6) | 19 (+- 4) | 43 (+- 6) | 107 |

All methods on identical circuits (mean +/- std):

| method | RX | RZ | CZ | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|
| **numeric_len** | 68.1 (+- 6.6) | 43.6 (+- 5.6) | 51.6 (+- 4.4) | 163.4 (+- 12.8) | LOSE (+56.4) | 51.6 | 83.6 |
| **numeric_cost** | 71.6 (+- 6.5) | 48.2 (+- 6.6) | 54.6 (+- 3.8) | 174.5 (+- 14.2) | LOSE (+67.5) | 54.6 | 43.7 |

Baseline fidelity check (our means vs paper's reported baseline means):

| baseline | our total | paper total | delta |
|---|---:|---:|---:|
| qiskit_l1 | (skipped) | 196 | -- |
| qiskit_l2 | (skipped) | 149 | -- |
| qiskit_l3 | (skipped) | 149 | -- |
| bqskit_l2 | (skipped) | 214 | -- |
| bqskit_l3 | (skipped) | 164 | -- |
| bqskit_l4 | (skipped) | 168 | -- |

Verdict notes:

- `numeric_len`: total LOSE vs paper (163.4 vs 107, +56.4); two-qubit LOSE (51.6 vs paper 43); equivalence pass rate 1.000; best 142.
- `numeric_cost`: total LOSE vs paper (174.5 vs 107, +67.5); two-qubit LOSE (54.6 vs paper 43); equivalence pass rate 1.000; best 149.

Note: NISQ inputs are CZ-weighted (weights RX:1, RZ:1, CZ:2) to match the paper's
Table 7 input composition (RX~108, RZ~109, CZ~82). The RZ-across-CZ pass is enabled
(iterated to a fixpoint).

Timing caveat: the "time (s)" column is the per-circuit budget cap -- each reducer
loops until its budget is exhausted. It is a cutoff, not a convergence time, and is
not directly comparable to the paper's Table 2 (a different task: reducing 100-gate
circuits to ~50, ~38 s for their best variant).
