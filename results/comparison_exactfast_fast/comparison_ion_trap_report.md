# Comparison benchmark report

- Gate set: `ion_trap` (paper Table 6, 4 qubits, length 300)
- Circuits per method: 100, per-circuit budget: 30.0s
- Generated: 2026-09-08 22:53:07
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RY | RZ | RXX | total |
|---|---|---|---|---|---|
| **paper** | 10 (+- 3) | 29 (+- 6) | 29 (+- 5) | 43 (+- 8) | 111 |

All methods on identical circuits (mean +/- std):

| method | RX | RY | RZ | RXX | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|---|
| **exact_len** | 2.4 (+- 1.3) | 18.2 (+- 5.1) | 17.2 (+- 4.6) | 28.2 (+- 6.7) | 66.1 (+- 14.4) | WIN (-44.9) | 28.2 | 30.0 |
| **exact_cost** | 2.9 (+- 1.6) | 18.8 (+- 5.1) | 19.2 (+- 5.0) | 25.9 (+- 5.1) | 66.7 (+- 12.9) | WIN (-44.3) | 25.9 | 30.1 |

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

- `exact_len`: total WIN vs paper (66.1 vs 111, -44.9); two-qubit WIN (28.2 vs paper 43); equivalence pass rate 1.000; best 32.
- `exact_cost`: total WIN vs paper (66.7 vs 111, -44.3); two-qubit WIN (25.9 vs paper 43); equivalence pass rate 1.000; best 37.

Timing caveat: the "time (s)" column is the per-circuit budget cap -- each reducer
loops until its budget is exhausted. It is a cutoff, not a convergence time, and is
not directly comparable to the paper's Table 2 (a different task: reducing 100-gate
circuits to ~50, ~38 s for their best variant).
