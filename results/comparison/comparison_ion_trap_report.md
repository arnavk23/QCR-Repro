# Comparison benchmark report

- Gate set: `ion_trap` (paper Table 6, 4 qubits, length 300)
- Circuits per method: 100, per-circuit budget: 30.0s
- Generated: 2026-08-13 03:55:49
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RY | RZ | RXX | total |
|---|---|---|---|---|---|
| **paper** | 10 (+- 3) | 29 (+- 6) | 29 (+- 5) | 43 (+- 8) | 111 |

All methods on identical circuits (mean +/- std):

| method | RX | RY | RZ | RXX | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|---|
| **exact_len** | 3.8 (+- 1.7) | 19.2 (+- 4.9) | 19.5 (+- 4.9) | 30.9 (+- 6.4) | 73.5 (+- 14.4) | WIN (-37.5) | 30.9 | 30.1 |
| **exact_cost** | 4.0 (+- 1.7) | 20.2 (+- 4.7) | 20.1 (+- 4.7) | 27.2 (+- 4.5) | 71.6 (+- 11.8) | WIN (-39.4) | 27.2 | 30.2 |
| **qiskit_l1** | 26.2 (+- 4.6) | 39.2 (+- 5.2) | 43.1 (+- 4.9) | 58.5 (+- 7.1) | 167.1 (+- 14.7) | base (paper 196) | 58.5 | 0.0 |
| **qiskit_l2** | 37.1 (+- 6.7) | 36.5 (+- 6.2) | 43.2 (+- 5.2) | 45.2 (+- 6.2) | 162.0 (+- 14.7) | base (paper 204) | 45.2 | 0.0 |
| **qiskit_l3** | 36.8 (+- 6.8) | 36.1 (+- 6.4) | 43.1 (+- 5.1) | 44.8 (+- 6.3) | 160.8 (+- 14.8) | base (paper 204) | 44.8 | 0.0 |

Baseline fidelity check (our means vs paper's reported baseline means):

| baseline | our total | paper total | delta |
|---|---:|---:|---:|
| qiskit_l1 | 167.1 | 196 | -28.9 |
| qiskit_l2 | 162.0 | 204 | -42.0 |
| qiskit_l3 | 160.8 | 204 | -43.2 |
| bqskit_l2 | (skipped) | 154 | -- |
| bqskit_l3 | (skipped) | 129 | -- |
| bqskit_l4 | (skipped) | 126 | -- |

Verdict notes:

- `exact_len`: total WIN vs paper (73.5 vs 111, -37.5); two-qubit WIN (30.9 vs paper 43); equivalence pass rate 1.000; best 41.
- `exact_cost`: total WIN vs paper (71.6 vs 111, -39.4); two-qubit WIN (27.2 vs paper 43); equivalence pass rate 1.000; best 45.

Timing caveat: the "time (s)" column is the per-circuit budget cap -- each reducer
loops until its budget is exhausted. It is a cutoff, not a convergence time, and is
not directly comparable to the paper's Table 2 (a different task: reducing 100-gate
circuits to ~50, ~38 s for their best variant).
