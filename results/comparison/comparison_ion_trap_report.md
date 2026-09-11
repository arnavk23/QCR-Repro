# Comparison benchmark report

- Gate set: `ion_trap` (paper Table 6, 4 qubits, length 300)
- Circuits per method: 8, per-circuit budget: 10.0s
- Generated: 2026-09-10 22:30:21
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RY | RZ | RXX | total |
|---|---|---|---|---|---|
| **paper** | 10 (+- 3) | 29 (+- 6) | 29 (+- 5) | 43 (+- 8) | 111 |

All methods on identical circuits (mean +/- std):

| method | RX | RY | RZ | RXX | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|---|
| **exact_len** | 5.0 (+- 1.2) | 22.5 (+- 5.8) | 20.5 (+- 4.7) | 34.5 (+- 5.4) | 82.5 (+- 13.0) | WIN (-28.5) | 34.5 | 10.1 |
| **exact_cost** | 6.4 (+- 2.8) | 23.4 (+- 2.6) | 24.4 (+- 5.9) | 30.4 (+- 3.8) | 84.5 (+- 11.1) | WIN (-26.5) | 30.4 | 10.3 |
| **qiskit_l1** | 25.5 (+- 5.3) | 39.4 (+- 4.9) | 43.0 (+- 3.6) | 59.1 (+- 8.5) | 167.0 (+- 15.7) | base (paper 196) | 59.1 | 0.0 |
| **qiskit_l2** | 37.2 (+- 7.7) | 36.8 (+- 5.4) | 41.0 (+- 5.6) | 48.5 (+- 5.2) | 163.5 (+- 13.9) | base (paper 204) | 48.5 | 0.0 |
| **qiskit_l3** | 37.1 (+- 7.8) | 36.6 (+- 5.5) | 40.9 (+- 5.5) | 48.4 (+- 5.2) | 163.0 (+- 14.2) | base (paper 204) | 48.4 | 0.0 |
| **bqskit_l2** | 55.3 (+- 3.7) | 0.0 (+- 0.0) | 79.0 (+- 4.5) | 42.3 (+- 2.6) | 176.7 (+- 8.7) | base (paper 154) | 42.3 | 0.0 |

Baseline fidelity check (our means vs paper's reported baseline means):

| baseline | our total | paper total | delta |
|---|---:|---:|---:|
| qiskit_l1 | 167.0 | 196 | -29.0 |
| qiskit_l2 | 163.5 | 204 | -40.5 |
| qiskit_l3 | 163.0 | 204 | -41.0 |
| bqskit_l2 | 176.7 | 154 | +22.7 |
| bqskit_l3 | (skipped) | 129 | -- |
| bqskit_l4 | (skipped) | 126 | -- |

Verdict notes:

- `exact_len`: total WIN vs paper (82.5 vs 111, -28.5); two-qubit WIN (34.5 vs paper 43); equivalence pass rate 1.000; best 61.
- `exact_cost`: total WIN vs paper (84.5 vs 111, -26.5); two-qubit WIN (30.4 vs paper 43); equivalence pass rate 1.000; best 60.

Timing caveat: the "time (s)" column is the per-circuit budget cap -- each reducer
loops until its budget is exhausted. It is a cutoff, not a convergence time, and is
not directly comparable to the paper's Table 2 (a different task: reducing 100-gate
circuits to ~50, ~38 s for their best variant).
