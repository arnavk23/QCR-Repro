# Comparison benchmark report

- Gate set: `nisq` (paper Table 7, 4 qubits, length 300)
- Circuits per method: 100, per-circuit budget: 60.0s
- Generated: 2026-09-16 22:25:32
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RZ | CZ | total |
|---|---|---|---|---|
| **paper** | 45 (+- 6) | 19 (+- 4) | 43 (+- 6) | 107 |

All methods on identical circuits (mean +/- std):

| method | RX | RZ | CZ | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|
| **numeric_len** | 67.5 (+- 5.7) | 38.7 (+- 4.7) | 49.1 (+- 5.4) | 155.3 (+- 12.5) | LOSE (+48.3) | 49.1 | 60.8 |
| **numeric_cost** | 67.6 (+- 5.7) | 38.6 (+- 4.6) | 48.9 (+- 5.4) | 155.1 (+- 12.4) | LOSE (+48.1) | 48.9 | 61.0 |
| **qiskit_l1** | 62.1 (+- 4.8) | 70.1 (+- 5.4) | 69.7 (+- 7.1) | 201.8 (+- 12.0) | base (paper 196) | 69.7 | 0.0 |
| **qiskit_l2** | 63.1 (+- 5.0) | 44.0 (+- 4.4) | 50.8 (+- 5.6) | 157.8 (+- 11.6) | base (paper 149) | 50.8 | 0.0 |
| **qiskit_l3** | 63.0 (+- 5.0) | 44.0 (+- 4.4) | 50.8 (+- 5.6) | 157.8 (+- 11.6) | base (paper 149) | 50.8 | 0.0 |

Baseline fidelity check (our means vs paper's reported baseline means):

| baseline | our total | paper total | delta |
|---|---:|---:|---:|
| qiskit_l1 | 201.8 | 196 | +5.8 |
| qiskit_l2 | 157.8 | 149 | +8.8 |
| qiskit_l3 | 157.8 | 149 | +8.8 |
| bqskit_l2 | (skipped) | 214 | -- |
| bqskit_l3 | (skipped) | 164 | -- |
| bqskit_l4 | (skipped) | 168 | -- |

Verdict notes:

- `numeric_len`: total LOSE vs paper (155.3 vs 107, +48.3); two-qubit LOSE (49.1 vs paper 43); equivalence pass rate 1.000; best 126.
- `numeric_cost`: total LOSE vs paper (155.1 vs 107, +48.1); two-qubit LOSE (48.9 vs paper 43); equivalence pass rate 1.000; best 126.

Note: NISQ inputs are CZ-weighted (weights RX:1, RZ:1, CZ:2) to match the paper's
Table 7 input composition (RX~108, RZ~109, CZ~82). The RZ-across-CZ pass is enabled
(iterated to a fixpoint).
Exact/numeric hybrid enabled: Clifford-only windows (RX/RZ at +/-pi/2, CZ) are
reduced by the exact symplectic engine at deep graph depths; only non-Clifford
(pi/4) windows hit the numeric database.

Timing caveat: the "time (s)" column is the per-circuit budget cap -- each reducer
loops until its budget is exhausted. It is a cutoff, not a convergence time, and is
not directly comparable to the paper's Table 2 (a different task: reducing 100-gate
circuits to ~50, ~38 s for their best variant).
