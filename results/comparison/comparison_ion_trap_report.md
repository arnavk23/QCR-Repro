# Comparison benchmark report

- Gate set: `ion_trap` (paper Table 6, 4 qubits, length 300)
- Circuits per method: 24, per-circuit budget: 15.0s
- Generated: 2026-08-09 02:36:46
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RY | RZ | RXX | total |
|---|---|---|---|---|---|
| **paper** | 10 (+- 3) | 29 (+- 6) | 29 (+- 5) | 43 (+- 8) | 111 |

All methods on identical circuits (mean +/- std):

| method | RX | RY | RZ | RXX | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|---|
| **exact_len** | 4.5 (+- 1.8) | 20.3 (+- 5.2) | 22.0 (+- 5.2) | 32.7 (+- 6.7) | 79.5 (+- 15.5) | WIN (-31.5) | 32.7 | 15.2 |
| **exact_cost** | 4.5 (+- 2.1) | 21.5 (+- 4.8) | 19.6 (+- 5.9) | 27.4 (+- 5.0) | 73.0 (+- 14.6) | WIN (-38.0) | 27.4 | 15.2 |
| **qiskit_l1** | 25.8 (+- 4.4) | 38.7 (+- 5.7) | 41.0 (+- 4.4) | 56.6 (+- 8.0) | 162.2 (+- 14.9) | base (paper 196) | 56.6 | 0.0 |
| **qiskit_l2** | 37.7 (+- 5.9) | 35.8 (+- 5.7) | 41.2 (+- 5.2) | 45.7 (+- 6.2) | 160.4 (+- 13.0) | base (paper 204) | 45.7 | 0.0 |
| **qiskit_l3** | 37.6 (+- 5.9) | 35.6 (+- 5.8) | 40.9 (+- 4.9) | 45.5 (+- 6.3) | 159.5 (+- 13.4) | base (paper 204) | 45.5 | 0.0 |

Baseline fidelity check (our means vs paper's reported baseline means):

| baseline | our total | paper total | delta |
|---|---:|---:|---:|
| qiskit_l1 | 162.2 | 196 | -33.8 |
| qiskit_l2 | 160.4 | 204 | -43.6 |
| qiskit_l3 | 159.5 | 204 | -44.5 |
| bqskit_l2 | (skipped) | 154 | -- |
| bqskit_l3 | (skipped) | 129 | -- |
| bqskit_l4 | (skipped) | 126 | -- |

Verdict notes:

- `exact_len`: total WIN vs paper (79.5 vs 111, -31.5); two-qubit WIN (32.7 vs paper 43); equivalence pass rate 1.000; best 53.
- `exact_cost`: total WIN vs paper (73.0 vs 111, -38.0); two-qubit WIN (27.4 vs paper 43); equivalence pass rate 1.000; best 49.
