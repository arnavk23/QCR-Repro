# Pre-pass + batched sweep comparison report

- Gate set: `ion_trap` (paper Table 6, 4 qubits, length 300)
- Circuits per method: 8, per-circuit budget: 10.0s
- Exacter: numeric 1e-5 (input vs output unitary up to global phase)
- Paper 'Ours' (mean over 100 runs): total 111, two-qubit 43

| method | input | after prepass | final total | two-qubit | vs paper Ours | prepass removed | runtime (s) | exact |
|---|---|---|---|---|---|---|---|---|
| **baseline** | 300 | 300 | 90.8 (+/- 14.7) | 37.4 | WIN (-18.2%) | 0.0 | 10.1 | OK |
| **prepass** | 300 | 284 | 92.1 (+/- 15.8) | 36.8 | WIN (-17.0%) | 15.5 | 10.1 | OK |
| **prepass+batched** | 300 | 284 | 89.5 (+/- 14.3) | 36.9 | WIN (-19.4%) | 15.5 | 10.0 | OK |

## Effect of the pre-passes

- `algebraic_merge`: fuses adjacent same-axis rotations whose sum snaps to a pool angle;drops exact cancellations (e.g. RZ(a)RZ(-a) = I).
- `zx_cancellations`: drops adjacent same-pair CZ CZ = I; for the diagonal-CZ pools, gathers RZ gates across CZ so runs become adjacent and fuse.
- `apply_prepass` is applied to fixpoint *before* the database loop; the output stays in the discrete pool, so the DB loop is unaffected except that the input is shorter.

## Wall-clock vs the paper (contextual)

- Paper Table 2 (reported, 100 length-100 circuits, 3-4 qubit): V1-random 199.0s, V2-db 55.0s, V3-rf 38.0s.
- Our pipelines above ran on identical random circuits of length 300 on this machine; the same budget was given to every method, so the comparison of interest is the final length at fixed budget.

## Batched sweep vs scalar sweep (sweep-only microbenchmark)

- The vectorized batched sweep is bit-identical to the scalar sweep (verified in scripts/check_batched_matches_scalar.py) and measures ~1.5-1.7x faster on the length-300 ion-trap fixpoint in our microbenchmark; the remaining per-window cost is the SHA-256 digest.