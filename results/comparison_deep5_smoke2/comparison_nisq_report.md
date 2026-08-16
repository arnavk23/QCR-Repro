# Comparison benchmark report

- Gate set: `nisq` (paper Table 7, 4 qubits, length 300)
- Circuits per method: 8, per-circuit budget: 20.0s
- Generated: 2026-08-17 02:38:49
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RZ | CZ | total |
|---|---|---|---|---|
| **paper** | 45 (+- 6) | 19 (+- 4) | 43 (+- 6) | 107 |

All methods on identical circuits (mean +/- std):

| method | RX | RZ | CZ | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|
| **numeric_len** | 67.8 (+- 6.8) | 43.4 (+- 5.7) | 51.4 (+- 4.5) | 162.5 (+- 13.2) | LOSE (+55.5) | 51.4 | 22.3 |
| **numeric_cost** | 68.4 (+- 6.7) | 43.2 (+- 6.3) | 51.1 (+- 3.9) | 162.8 (+- 13.5) | LOSE (+55.8) | 51.1 | 21.8 |

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

- `numeric_len`: total LOSE vs paper (162.5 vs 107, +55.5); two-qubit LOSE (51.4 vs paper 43); equivalence pass rate 1.000; best 139.
- `numeric_cost`: total LOSE vs paper (162.8 vs 107, +55.8); two-qubit LOSE (51.1 vs paper 43); equivalence pass rate 1.000; best 141.

Note: NISQ inputs are CZ-weighted (weights RX:1, RZ:1, CZ:2) to match the paper's
Table 7 input composition (RX~108, RZ~109, CZ~82). The RZ-across-CZ pass is enabled
(iterated to a fixpoint).

Timing caveat: the "time (s)" column is the per-circuit budget cap -- each reducer
loops until its budget is exhausted. It is a cutoff, not a convergence time, and is
not directly comparable to the paper's Table 2 (a different task: reducing 100-gate
circuits to ~50, ~38 s for their best variant).
