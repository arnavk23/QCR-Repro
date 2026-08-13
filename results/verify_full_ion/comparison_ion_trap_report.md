# Comparison benchmark report

- Gate set: `ion_trap` (paper Table 6, 4 qubits, length 300)
- Circuits per method: 100, per-circuit budget: 15.0s
- Generated: 2026-08-13 00:57:49
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RY | RZ | RXX | total |
|---|---|---|---|---|---|
| **paper** | 10 (+- 3) | 29 (+- 6) | 29 (+- 5) | 43 (+- 8) | 111 |

All methods on identical circuits (mean +/- std):

| method | RX | RY | RZ | RXX | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|---|
| **exact_len** | 7.5 (+- 3.3) | 23.2 (+- 5.3) | 23.3 (+- 5.6) | 36.5 (+- 6.8) | 90.5 (+- 17.4) | WIN (-20.5) | 36.5 | 15.5 |
| **exact_cost** | 9.2 (+- 3.7) | 25.2 (+- 5.5) | 24.2 (+- 6.0) | 32.5 (+- 5.9) | 91.0 (+- 16.7) | WIN (-20.0) | 32.5 | 15.8 |

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

- `exact_len`: total WIN vs paper (90.5 vs 111, -20.5); two-qubit WIN (36.5 vs paper 43); equivalence pass rate 1.000; best 52.
- `exact_cost`: total WIN vs paper (91.0 vs 111, -20.0); two-qubit WIN (32.5 vs paper 43); equivalence pass rate 1.000; best 59.
