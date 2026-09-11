# Comparison benchmark report

- Gate set: `ion_trap` (paper Table 6, 4 qubits, length 300)
- Circuits per method: 30, per-circuit budget: 30.0s
- Generated: 2026-09-11 00:04:11
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RY | RZ | RXX | total |
|---|---|---|---|---|---|
| **paper** | 10 (+- 3) | 29 (+- 6) | 29 (+- 5) | 43 (+- 8) | 111 |

All methods on identical circuits (mean +/- std):

| method | RX | RY | RZ | RXX | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|---|
| **exact_len** | 2.5 (+- 1.2) | 18.7 (+- 4.6) | 17.2 (+- 4.2) | 28.0 (+- 5.7) | 66.4 (+- 12.0) | WIN (-44.6) | 28.0 | 30.1 |
| **exact_cost** | 3.4 (+- 2.0) | 19.4 (+- 4.1) | 19.6 (+- 5.5) | 25.4 (+- 4.9) | 67.8 (+- 12.3) | WIN (-43.2) | 25.4 | 30.1 |
| **qiskit_l1** | 25.5 (+- 4.5) | 38.7 (+- 5.3) | 41.4 (+- 4.6) | 57.1 (+- 7.4) | 162.7 (+- 14.1) | base (paper 196) | 57.1 | 0.0 |
| **qiskit_l2** | 35.8 (+- 6.5) | 34.8 (+- 6.2) | 41.0 (+- 4.8) | 45.3 (+- 5.9) | 157.0 (+- 12.6) | base (paper 204) | 45.3 | 0.0 |
| **qiskit_l3** | 35.4 (+- 6.7) | 34.5 (+- 6.2) | 40.9 (+- 4.6) | 44.8 (+- 6.2) | 155.6 (+- 13.6) | base (paper 204) | 44.8 | 0.0 |
| **bqskit_l2** | 49.7 (+- 8.7) | 0.0 (+- 0.0) | 68.8 (+- 10.4) | 36.6 (+- 5.7) | 155.1 (+- 23.7) | base (paper 154) | 36.6 | 0.0 |
| **bqskit_l3** | 43.1 (+- 7.7) | 0.0 (+- 0.0) | 59.6 (+- 9.2) | 30.3 (+- 5.3) | 133.0 (+- 21.1) | base (paper 129) | 30.3 | 0.0 |
| **bqskit_l4** | 39.3 (+- 7.2) | 0.0 (+- 0.0) | 55.9 (+- 8.8) | 26.4 (+- 4.2) | 121.6 (+- 18.9) | base (paper 126) | 26.4 | 0.0 |

Baseline fidelity check (our means vs paper's reported baseline means):

| baseline | our total | paper total | delta |
|---|---:|---:|---:|
| qiskit_l1 | 162.7 | 196 | -33.3 |
| qiskit_l2 | 157.0 | 204 | -47.0 |
| qiskit_l3 | 155.6 | 204 | -48.4 |
| bqskit_l2 | 155.1 | 154 | +1.1 |
| bqskit_l3 | 133.0 | 129 | +4.0 |
| bqskit_l4 | 121.6 | 126 | -4.4 |

Verdict notes:

- `exact_len`: total WIN vs paper (66.4 vs 111, -44.6); two-qubit WIN (28.0 vs paper 43); equivalence pass rate 1.000; best 46.
- `exact_cost`: total WIN vs paper (67.8 vs 111, -43.2); two-qubit WIN (25.4 vs paper 43); equivalence pass rate 1.000; best 38.

Timing caveat: the "time (s)" column is the per-circuit budget cap -- each reducer
loops until its budget is exhausted. It is a cutoff, not a convergence time, and is
not directly comparable to the paper's Table 2 (a different task: reducing 100-gate
circuits to ~50, ~38 s for their best variant).
