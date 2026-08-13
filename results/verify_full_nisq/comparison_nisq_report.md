# Comparison benchmark report

- Gate set: `nisq` (paper Table 7, 4 qubits, length 300)
- Circuits per method: 100, per-circuit budget: 15.0s
- Generated: 2026-08-13 01:07:31
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RZ | CZ | total |
|---|---|---|---|---|
| **paper** | 45 (+- 6) | 19 (+- 4) | 43 (+- 6) | 107 |

All methods on identical circuits (mean +/- std):

| method | RX | RZ | CZ | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|
| **numeric_len** | 69.3 (+- 6.0) | 43.8 (+- 5.4) | 50.9 (+- 5.6) | 164.0 (+- 12.7) | LOSE (+57.0) | 50.9 | 15.8 |
| **numeric_cost** | 69.7 (+- 6.0) | 44.2 (+- 5.6) | 50.7 (+- 5.6) | 164.6 (+- 13.3) | LOSE (+57.6) | 50.7 | 16.1 |

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

- `numeric_len`: total LOSE vs paper (164.0 vs 107, +57.0); two-qubit LOSE (50.9 vs paper 43); equivalence pass rate 1.000; best 135.
- `numeric_cost`: total LOSE vs paper (164.6 vs 107, +57.6); two-qubit LOSE (50.7 vs paper 43); equivalence pass rate 1.000; best 133.

Note: NISQ inputs are CZ-weighted (weights RX:1, RZ:1, CZ:2) to match the paper's
Table 7 input composition (RX~108, RZ~109, CZ~82). The RZ-across-CZ pass is enabled.
