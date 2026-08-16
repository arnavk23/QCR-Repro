# Comparison benchmark report

- Gate set: `nisq` (paper Table 7, 4 qubits, length 300)
- Circuits per method: 100, per-circuit budget: 60.0s
- Generated: 2026-08-17 03:24:44
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RZ | CZ | total |
|---|---|---|---|---|
| **paper** | 45 (+- 6) | 19 (+- 4) | 43 (+- 6) | 107 |

All methods on identical circuits (mean +/- std):

| method | RX | RZ | CZ | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|
| **numeric_len** | 68.0 (+- 6.1) | 42.5 (+- 5.1) | 50.2 (+- 5.4) | 160.7 (+- 13.0) | LOSE (+53.7) | 50.2 | 63.8 |
| **numeric_cost** | 68.7 (+- 6.2) | 42.9 (+- 5.3) | 49.5 (+- 5.5) | 161.0 (+- 13.1) | LOSE (+54.0) | 49.5 | 64.1 |

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

- `numeric_len`: total LOSE vs paper (160.7 vs 107, +53.7); two-qubit LOSE (50.2 vs paper 43); equivalence pass rate 1.000; best 132.
- `numeric_cost`: total LOSE vs paper (161.0 vs 107, +54.0); two-qubit LOSE (49.5 vs paper 43); equivalence pass rate 1.000; best 130.

Note: NISQ inputs are CZ-weighted (weights RX:1, RZ:1, CZ:2) to match the paper's
Table 7 input composition (RX~108, RZ~109, CZ~82). The RZ-across-CZ pass is enabled
(iterated to a fixpoint).

Timing caveat: the "time (s)" column is the per-circuit budget cap -- each reducer
loops until its budget is exhausted. It is a cutoff, not a convergence time, and is
not directly comparable to the paper's Table 2 (a different task: reducing 100-gate
circuits to ~50, ~38 s for their best variant).
