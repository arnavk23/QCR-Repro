# Comparison benchmark report

- Gate set: `ion_trap` (paper Table 6, 4 qubits, length 300)
- Circuits per method: 100, per-circuit budget: 15.0s
- Generated: 2026-08-09 04:19:24
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RY | RZ | RXX | total |
|---|---|---|---|---|---|
| **paper** | 10 (+- 3) | 29 (+- 6) | 29 (+- 5) | 43 (+- 8) | 111 |

All methods on identical circuits (mean +/- std):

| method | RX | RY | RZ | RXX | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|---|
| **exact_len** | 4.6 (+- 1.9) | 20.8 (+- 5.1) | 20.9 (+- 5.4) | 32.7 (+- 6.4) | 79.0 (+- 15.3) | WIN (-32.0) | 32.7 | 15.1 |
| **exact_cost** | 5.7 (+- 2.4) | 22.1 (+- 5.4) | 21.8 (+- 4.9) | 29.0 (+- 4.9) | 78.5 (+- 13.5) | WIN (-32.5) | 29.0 | 15.2 |

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

- `exact_len`: total WIN vs paper (79.0 vs 111, -32.0); two-qubit WIN (32.7 vs paper 43); equivalence pass rate 1.000; best 42.
- `exact_cost`: total WIN vs paper (78.5 vs 111, -32.5); two-qubit WIN (29.0 vs paper 43); equivalence pass rate 1.000; best 51.
