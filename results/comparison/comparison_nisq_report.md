# Comparison benchmark report

- Gate set: `nisq` (paper Table 7, 4 qubits, length 300)
- Circuits per method: 12, per-circuit budget: 15.0s
- Generated: 2026-08-09 02:56:50
- Verifier: exact symplectic (ion trap) / numeric 1e-5 (NISQ)

Paper reference ('Ours', mean +/- std over 100 runs):

| method | RX | RZ | CZ | total |
|---|---|---|---|---|
| **paper** | 45 (+- 6) | 19 (+- 4) | 43 (+- 6) | 107 |

All methods on identical circuits (mean +/- std):

| method | RX | RZ | CZ | total | vs paper Ours | two-qubit | time (s) |
|---|---|---|---|---|---|---|---|
| **numeric_len** | 65.8 (+- 6.0) | 42.2 (+- 4.8) | 50.3 (+- 5.6) | 158.4 (+- 13.4) | LOSE (+51.4) | 50.3 | 47.0 |
| **numeric_cost** | 66.9 (+- 6.1) | 42.6 (+- 4.6) | 50.2 (+- 5.3) | 159.8 (+- 12.9) | LOSE (+52.8) | 50.2 | 47.0 |
| **qiskit_l1** | 60.5 (+- 4.0) | 68.8 (+- 4.8) | 71.4 (+- 5.9) | 200.7 (+- 11.3) | base (paper 196) | 71.4 | 0.0 |
| **qiskit_l2** | 61.2 (+- 3.8) | 41.8 (+- 4.3) | 50.4 (+- 5.6) | 153.3 (+- 10.5) | base (paper 149) | 50.4 | 0.0 |
| **qiskit_l3** | 61.2 (+- 3.8) | 41.8 (+- 4.3) | 50.4 (+- 5.6) | 153.3 (+- 10.5) | base (paper 149) | 50.4 | 0.0 |

Baseline fidelity check (our means vs paper's reported baseline means):

| baseline | our total | paper total | delta |
|---|---:|---:|---:|
| qiskit_l1 | 200.7 | 196 | +4.7 |
| qiskit_l2 | 153.3 | 149 | +4.3 |
| qiskit_l3 | 153.3 | 149 | +4.3 |
| bqskit_l2 | (skipped) | 214 | -- |
| bqskit_l3 | (skipped) | 164 | -- |
| bqskit_l4 | (skipped) | 168 | -- |

Verdict notes:

- `numeric_len`: total LOSE vs paper (158.4 vs 107, +51.4); two-qubit LOSE (50.3 vs paper 43); equivalence pass rate 1.000; best 137.
- `numeric_cost`: total LOSE vs paper (159.8 vs 107, +52.8); two-qubit LOSE (50.2 vs paper 43); equivalence pass rate 1.000; best 139.

Note: NISQ inputs are CZ-weighted (weights RX:1, RZ:1, CZ:2) to match the paper's
Table 7 input composition (RX~108, RZ~109, CZ~82). The RZ-across-CZ pass is enabled.
