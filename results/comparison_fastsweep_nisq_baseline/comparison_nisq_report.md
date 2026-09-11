# Comparison benchmark report

- Gate set: `nisq` (paper Table 7, 4 qubits, length 300)
- Circuits per method: 100, per-circuit budget: 60.0s
- Generated: 2026-09-08 15:11:55
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RZ | CZ | total |
|---|---|---|---|---|
| **paper** | 45 (+- 6) | 19 (+- 4) | 43 (+- 6) | 107 |

All methods on identical circuits (mean +/- std):

| method | RX | RZ | CZ | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|
| **numeric_len** | 68.2 (+- 5.9) | 42.2 (+- 5.0) | 50.3 (+- 5.7) | 160.7 (+- 12.5) | LOSE (+53.7) | 50.3 | 61.0 |
| **numeric_cost** | 68.8 (+- 5.7) | 42.4 (+- 5.3) | 49.9 (+- 5.6) | 161.1 (+- 12.5) | LOSE (+54.1) | 49.9 | 61.0 |

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

- `numeric_len`: total LOSE vs paper (160.7 vs 107, +53.7); two-qubit LOSE (50.3 vs paper 43); equivalence pass rate 1.000; best 130.
- `numeric_cost`: total LOSE vs paper (161.1 vs 107, +54.1); two-qubit LOSE (49.9 vs paper 43); equivalence pass rate 1.000; best 131.

Note: NISQ inputs are CZ-weighted (weights RX:1, RZ:1, CZ:2) to match the paper's
Table 7 input composition (RX~108, RZ~109, CZ~82). The RZ-across-CZ pass is enabled
(iterated to a fixpoint).

Timing caveat: the "time (s)" column is the per-circuit budget cap -- each reducer
loops until its budget is exhausted. It is a cutoff, not a convergence time, and is
not directly comparable to the paper's Table 2 (a different task: reducing 100-gate
circuits to ~50, ~38 s for their best variant).
