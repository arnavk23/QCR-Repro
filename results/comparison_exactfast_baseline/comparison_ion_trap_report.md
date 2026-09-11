# Comparison benchmark report

- Gate set: `ion_trap` (paper Table 6, 4 qubits, length 300)
- Circuits per method: 100, per-circuit budget: 30.0s
- Generated: 2026-09-08 22:38:45
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RY | RZ | RXX | total |
|---|---|---|---|---|---|
| **paper** | 10 (+- 3) | 29 (+- 6) | 29 (+- 5) | 43 (+- 8) | 111 |

All methods on identical circuits (mean +/- std):

| method | RX | RY | RZ | RXX | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|---|
| **exact_len** | 5.2 (+- 2.1) | 21.4 (+- 4.9) | 21.2 (+- 5.4) | 33.4 (+- 6.3) | 81.1 (+- 15.2) | WIN (-29.9) | 33.4 | 30.4 |
| **exact_cost** | 4.9 (+- 2.3) | 21.3 (+- 4.9) | 21.2 (+- 5.1) | 28.3 (+- 4.7) | 75.7 (+- 13.2) | WIN (-35.3) | 28.3 | 30.3 |

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

- `exact_len`: total WIN vs paper (81.1 vs 111, -29.9); two-qubit WIN (33.4 vs paper 43); equivalence pass rate 1.000; best 42.
- `exact_cost`: total WIN vs paper (75.7 vs 111, -35.3); two-qubit WIN (28.3 vs paper 43); equivalence pass rate 1.000; best 50.

Timing caveat: the "time (s)" column is the per-circuit budget cap -- each reducer
loops until its budget is exhausted. It is a cutoff, not a convergence time, and is
not directly comparable to the paper's Table 2 (a different task: reducing 100-gate
circuits to ~50, ~38 s for their best variant).
