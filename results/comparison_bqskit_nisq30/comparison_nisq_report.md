# Comparison benchmark report

- Gate set: `nisq` (paper Table 7, 4 qubits, length 300)
- Circuits per method: 30, per-circuit budget: 60.0s
- Generated: 2026-09-11 15:47:45
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RZ | CZ | total |
|---|---|---|---|---|
| **paper** | 45 (+- 6) | 19 (+- 4) | 43 (+- 6) | 107 |

All methods on identical circuits (mean +/- std):

| method | RX | RZ | CZ | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|
| **numeric_len** | 66.9 (+- 6.5) | 42.4 (+- 5.9) | 51.2 (+- 5.6) | 160.5 (+- 14.1) | LOSE (+53.5) | 51.2 | 61.4 |
| **numeric_cost** | 67.9 (+- 6.8) | 42.6 (+- 5.6) | 51.1 (+- 5.8) | 161.7 (+- 13.9) | LOSE (+54.7) | 51.1 | 61.3 |
| **qiskit_l1** | 61.6 (+- 4.9) | 68.8 (+- 4.5) | 69.4 (+- 6.7) | 199.8 (+- 12.4) | base (paper 196) | 69.4 | 0.0 |
| **qiskit_l2** | 63.2 (+- 5.0) | 43.2 (+- 4.9) | 51.3 (+- 6.1) | 157.7 (+- 12.3) | base (paper 149) | 51.3 | 0.0 |
| **qiskit_l3** | 63.1 (+- 4.9) | 43.2 (+- 4.9) | 51.3 (+- 6.1) | 157.6 (+- 12.2) | base (paper 149) | 51.3 | 0.0 |
| **bqskit_l2** | 72.7 (+- 7.1) | 91.0 (+- 11.7) | 65.2 (+- 8.1) | 228.9 (+- 19.3) | base (paper 214) | 65.2 | 0.0 |
| **bqskit_l3** | 70.9 (+- 7.4) | 91.5 (+- 11.4) | 50.8 (+- 7.2) | 213.2 (+- 23.5) | base (paper 164) | 50.8 | 0.0 |

Baseline fidelity check (our means vs paper's reported baseline means):

| baseline | our total | paper total | delta |
|---|---:|---:|---:|
| qiskit_l1 | 199.8 | 196 | +3.8 |
| qiskit_l2 | 157.7 | 149 | +8.7 |
| qiskit_l3 | 157.6 | 149 | +8.6 |
| bqskit_l2 | 228.9 | 214 | +14.9 |
| bqskit_l3 | 213.2 | 164 | +49.2 |
| bqskit_l4 | (skipped) | 168 | -- |

Verdict notes:

- `numeric_len`: total LOSE vs paper (160.5 vs 107, +53.5); two-qubit LOSE (51.2 vs paper 43); equivalence pass rate 1.000; best 132.
- `numeric_cost`: total LOSE vs paper (161.7 vs 107, +54.7); two-qubit LOSE (51.1 vs paper 43); equivalence pass rate 1.000; best 132.

Note: NISQ inputs are CZ-weighted (weights RX:1, RZ:1, CZ:2) to match the paper's
Table 7 input composition (RX~108, RZ~109, CZ~82). The RZ-across-CZ pass is enabled
(iterated to a fixpoint).

Timing caveat: the "time (s)" column is the per-circuit budget cap -- each reducer
loops until its budget is exhausted. It is a cutoff, not a convergence time, and is
not directly comparable to the paper's Table 2 (a different task: reducing 100-gate
circuits to ~50, ~38 s for their best variant).
