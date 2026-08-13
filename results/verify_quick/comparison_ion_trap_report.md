# Comparison benchmark report

- Gate set: `ion_trap` (paper Table 6, 4 qubits, length 300)
- Circuits per method: 8, per-circuit budget: 10.0s
- Generated: 2026-08-13 00:48:02
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RY | RZ | RXX | total |
|---|---|---|---|---|---|
| **paper** | 10 (+- 3) | 29 (+- 6) | 29 (+- 5) | 43 (+- 8) | 111 |

All methods on identical circuits (mean +/- std):

| method | RX | RY | RZ | RXX | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|---|
| **exact_len** | 6.5 (+- 1.1) | 25.9 (+- 4.9) | 25.9 (+- 5.3) | 39.9 (+- 6.7) | 98.1 (+- 15.1) | WIN (-12.9) | 39.9 | 10.2 |
| **exact_cost** | 9.8 (+- 1.7) | 24.2 (+- 5.1) | 27.5 (+- 5.9) | 33.5 (+- 4.2) | 95.0 (+- 12.4) | WIN (-16.0) | 33.5 | 10.3 |

Baseline fidelity check (our means vs paper's reported baseline means):

| baseline | our total | paper total | delta |
|---|---:|---:|---:|
| qiskit_l1 | (skipped) | 196 | -- |
| qiskit_l2 | (skipped) | 204 | -- |
| qiskit_l3 | (skipped) | 204 | -- |
| bqskit_l2 | (skipped) | 154 | -- |
| bqskit_l3 | (skipped) | 129 | -- |
| bqskit_l4 | (skipped) | 126 | -- |

Verdict notes:

- `exact_len`: total WIN vs paper (98.1 vs 111, -12.9); two-qubit WIN (39.9 vs paper 43); equivalence pass rate 1.000; best 62.
- `exact_cost`: total WIN vs paper (95.0 vs 111, -16.0); two-qubit WIN (33.5 vs paper 43); equivalence pass rate 1.000; best 69.
