# Comparison benchmark report

- Gate set: `nisq` (paper Table 7, 4 qubits, length 300)
- Circuits per method: 10, per-circuit budget: 60.0s
- Generated: 2026-09-11 08:11:16
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RZ | CZ | total |
|---|---|---|---|---|
| **paper** | 45 (+- 6) | 19 (+- 4) | 43 (+- 6) | 107 |

All methods on identical circuits (mean +/- std):

| method | RX | RZ | CZ | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|
| **numeric_len** | 65.2 (+- 7.1) | 40.6 (+- 5.6) | 50.1 (+- 5.2) | 155.9 (+- 14.5) | LOSE (+48.9) | 50.1 | 60.4 |
| **numeric_cost** | 66.5 (+- 7.3) | 41.9 (+- 5.8) | 50.7 (+- 5.4) | 159.1 (+- 15.3) | LOSE (+52.1) | 50.7 | 61.0 |
| **qiskit_l1** | 61.0 (+- 4.1) | 68.9 (+- 4.8) | 72.5 (+- 5.3) | 202.4 (+- 10.0) | base (paper 196) | 72.5 | 0.0 |
| **qiskit_l2** | 61.0 (+- 3.7) | 41.6 (+- 4.8) | 49.8 (+- 5.8) | 152.4 (+- 11.3) | base (paper 149) | 49.8 | 0.0 |
| **qiskit_l3** | 61.0 (+- 3.7) | 41.6 (+- 4.8) | 49.8 (+- 5.8) | 152.4 (+- 11.3) | base (paper 149) | 49.8 | 0.0 |
| **bqskit_l2** | 71.3 (+- 7.7) | 92.7 (+- 8.6) | 64.8 (+- 6.3) | 228.8 (+- 18.6) | base (paper 214) | 64.8 | 0.0 |
| **bqskit_l3** | 71.1 (+- 7.5) | 95.4 (+- 11.6) | 49.2 (+- 7.1) | 215.7 (+- 24.9) | base (paper 164) | 49.2 | 0.0 |

Baseline fidelity check (our means vs paper's reported baseline means):

| baseline | our total | paper total | delta |
|---|---:|---:|---:|
| qiskit_l1 | 202.4 | 196 | +6.4 |
| qiskit_l2 | 152.4 | 149 | +3.4 |
| qiskit_l3 | 152.4 | 149 | +3.4 |
| bqskit_l2 | 228.8 | 214 | +14.8 |
| bqskit_l3 | 215.7 | 164 | +51.7 |
| bqskit_l4 | (skipped) | 168 | -- |

Verdict notes:

- `numeric_len`: total LOSE vs paper (155.9 vs 107, +48.9); two-qubit LOSE (50.1 vs paper 43); equivalence pass rate 1.000; best 132.
- `numeric_cost`: total LOSE vs paper (159.1 vs 107, +52.1); two-qubit LOSE (50.7 vs paper 43); equivalence pass rate 1.000; best 135.

Note: NISQ inputs are CZ-weighted (weights RX:1, RZ:1, CZ:2) to match the paper's
Table 7 input composition (RX~108, RZ~109, CZ~82). The RZ-across-CZ pass is enabled
(iterated to a fixpoint).

Timing caveat: the "time (s)" column is the per-circuit budget cap -- each reducer
loops until its budget is exhausted. It is a cutoff, not a convergence time, and is
not directly comparable to the paper's Table 2 (a different task: reducing 100-gate
circuits to ~50, ~38 s for their best variant).
