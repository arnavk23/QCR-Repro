# Comparison benchmark report

- Gate set: `ion_trap` (paper Table 6, 4 qubits, length 300)
- Circuits per method: 100, per-circuit budget: 30.0s
- Generated: 2026-09-08 21:47:54
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RY | RZ | RXX | total |
|---|---|---|---|---|---|
| **paper** | 10 (+- 3) | 29 (+- 6) | 29 (+- 5) | 43 (+- 8) | 111 |

All methods on identical circuits (mean +/- std):

| method | RX | RY | RZ | RXX | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|---|
| **exact_len** | 7.8 (+- 6.0) | 22.4 (+- 6.1) | 22.6 (+- 6.5) | 35.1 (+- 7.8) | 87.9 (+- 22.2) | WIN (-23.1) | 35.1 | 1213.8 |
| **exact_cost** | 7.1 (+- 4.5) | 23.7 (+- 5.5) | 22.7 (+- 5.6) | 30.5 (+- 5.7) | 84.0 (+- 17.0) | WIN (-27.0) | 30.5 | 490.3 |
| **numeric_len** | 3.8 (+- 1.7) | 20.2 (+- 4.4) | 19.3 (+- 5.5) | 31.4 (+- 6.3) | 74.6 (+- 13.9) | WIN (-36.4) | 31.4 | 30.6 |

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

- `exact_len`: total WIN vs paper (87.9 vs 111, -23.1); two-qubit WIN (35.1 vs paper 43); equivalence pass rate 1.000; best 42.
- `exact_cost`: total WIN vs paper (84.0 vs 111, -27.0); two-qubit WIN (30.5 vs paper 43); equivalence pass rate 1.000; best 52.
- `numeric_len`: total WIN vs paper (74.6 vs 111, -36.4); two-qubit WIN (31.4 vs paper 43); equivalence pass rate 1.000; best 41.

Timing caveat: the "time (s)" column is the per-circuit budget cap -- each reducer
loops until its budget is exhausted. It is a cutoff, not a convergence time, and is
not directly comparable to the paper's Table 2 (a different task: reducing 100-gate
circuits to ~50, ~38 s for their best variant).
