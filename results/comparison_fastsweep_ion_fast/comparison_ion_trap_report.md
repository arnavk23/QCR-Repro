# Comparison benchmark report

- Gate set: `ion_trap` (paper Table 6, 4 qubits, length 300)
- Circuits per method: 100, per-circuit budget: 30.0s
- Generated: 2026-09-08 22:09:49
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RY | RZ | RXX | total |
|---|---|---|---|---|---|
| **paper** | 10 (+- 3) | 29 (+- 6) | 29 (+- 5) | 43 (+- 8) | 111 |

All methods on identical circuits (mean +/- std):

| method | RX | RY | RZ | RXX | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|---|
| **exact_len** | 4.9 (+- 2.2) | 21.1 (+- 5.0) | 21.1 (+- 5.4) | 33.1 (+- 6.4) | 80.2 (+- 15.3) | WIN (-30.8) | 33.1 | 30.4 |
| **exact_cost** | 6.4 (+- 2.7) | 23.1 (+- 5.2) | 22.4 (+- 5.1) | 30.0 (+- 5.1) | 82.0 (+- 14.1) | WIN (-29.0) | 30.0 | 30.7 |
| **numeric_len** | 4.0 (+- 1.8) | 19.6 (+- 5.3) | 19.5 (+- 4.3) | 31.0 (+- 6.1) | 74.0 (+- 13.8) | WIN (-37.0) | 31.0 | 30.6 |

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

- `exact_len`: total WIN vs paper (80.2 vs 111, -30.8); two-qubit WIN (33.1 vs paper 43); equivalence pass rate 1.000; best 42.
- `exact_cost`: total WIN vs paper (82.0 vs 111, -29.0); two-qubit WIN (30.0 vs paper 43); equivalence pass rate 1.000; best 52.
- `numeric_len`: total WIN vs paper (74.0 vs 111, -37.0); two-qubit WIN (31.0 vs paper 43); equivalence pass rate 1.000; best 44.

Timing caveat: the "time (s)" column is the per-circuit budget cap -- each reducer
loops until its budget is exhausted. It is a cutoff, not a convergence time, and is
not directly comparable to the paper's Table 2 (a different task: reducing 100-gate
circuits to ~50, ~38 s for their best variant).
