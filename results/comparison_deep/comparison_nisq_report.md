# Comparison benchmark report

- Gate set: `nisq` (paper Table 7, 4 qubits, length 300)
- Circuits per method: 100, per-circuit budget: 60.0s
- Generated: 2026-08-13 05:33:10
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RZ | CZ | total |
|---|---|---|---|---|
| **paper** | 45 (+- 6) | 19 (+- 4) | 43 (+- 6) | 107 |

All methods on identical circuits (mean +/- std):

| method | RX | RZ | CZ | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|
| **numeric_len** | 67.2 (+- 6.0) | 40.7 (+- 5.1) | 49.4 (+- 5.4) | 157.3 (+- 12.7) | LOSE (+50.3) | 49.4 | 61.1 |
| **numeric_cost** | 67.8 (+- 5.7) | 41.1 (+- 5.1) | 48.4 (+- 5.4) | 157.2 (+- 12.6) | LOSE (+50.2) | 48.4 | 61.2 |
| **qiskit_l1** | 62.1 (+- 4.8) | 70.1 (+- 5.4) | 69.7 (+- 7.1) | 201.8 (+- 12.0) | base (paper 196) | 69.7 | 0.0 |
| **qiskit_l2** | 63.1 (+- 5.0) | 44.1 (+- 4.4) | 50.8 (+- 5.6) | 158.0 (+- 11.5) | base (paper 149) | 50.8 | 0.0 |
| **qiskit_l3** | 63.0 (+- 5.0) | 44.1 (+- 4.5) | 50.8 (+- 5.6) | 157.9 (+- 11.6) | base (paper 149) | 50.8 | 0.0 |

Baseline fidelity check (our means vs paper's reported baseline means):

| baseline | our total | paper total | delta |
|---|---:|---:|---:|
| qiskit_l1 | 201.8 | 196 | +5.8 |
| qiskit_l2 | 158.0 | 149 | +9.0 |
| qiskit_l3 | 157.9 | 149 | +8.9 |
| bqskit_l2 | (skipped) | 214 | -- |
| bqskit_l3 | (skipped) | 164 | -- |
| bqskit_l4 | (skipped) | 168 | -- |

Verdict notes:

- `numeric_len`: total LOSE vs paper (157.3 vs 107, +50.3); two-qubit LOSE (49.4 vs paper 43); equivalence pass rate 1.000; best 125.
- `numeric_cost`: total LOSE vs paper (157.2 vs 107, +50.2); two-qubit LOSE (48.4 vs paper 43); equivalence pass rate 1.000; best 121.

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
