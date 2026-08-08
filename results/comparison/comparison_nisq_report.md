# Comparison benchmark report

- Gate set: `nisq` (paper Table 7, 4 qubits, length 300)
- Circuits per method: 100, per-circuit budget: 15.0s
- Generated: 2026-08-09 04:37:01
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RZ | CZ | total |
|---|---|---|---|---|
| **paper** | 45 (+- 6) | 19 (+- 4) | 43 (+- 6) | 107 |

All methods on identical circuits (mean +/- std):

| method | RX | RZ | CZ | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|
| **numeric_len** | 69.3 (+- 6.0) | 43.7 (+- 5.2) | 50.9 (+- 5.5) | 163.8 (+- 12.5) | LOSE (+56.8) | 50.9 | 15.7 |
| **numeric_cost** | 70.0 (+- 7.5) | 44.7 (+- 8.2) | 50.8 (+- 6.3) | 165.4 (+- 18.7) | LOSE (+58.4) | 50.8 | 16.0 |

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

- `numeric_len`: total LOSE vs paper (163.8 vs 107, +56.8); two-qubit LOSE (50.9 vs paper 43); equivalence pass rate 1.000; best 135.
- `numeric_cost`: total LOSE vs paper (165.4 vs 107, +58.4); two-qubit LOSE (50.8 vs paper 43); equivalence pass rate 1.000; best 131.

Note: NISQ inputs are CZ-weighted (weights RX:1, RZ:1, CZ:2) to match the paper's
Table 7 input composition (RX~108, RZ~109, CZ~82). The RZ-across-CZ pass is enabled.
