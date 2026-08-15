# Optimization-Driven Quantum Circuit Reduction

Python implementation of
**Bodo Rosenhahn, Tobias J Osborne, Christoph Hirche, "Optimization driven quantum circuit reduction," New J. Phys. 27 (2025) 104509**

covering the paper's local term-replacement scheme (variants V1–V3), plus a
bit-exact symplectic (Clifford) engine and a dependency-graph based sweep
extension that are not part of the paper.

Headline comparison (`results/comparison/`, 4 qubits, 300-gate random
circuits, identical circuits per method, 100 circuits per method):

| Gate set | Paper "Ours" | Our reducers | Difference |
|---|---|---|---|
| Ion trap (RX/RY/RZ/RXX) | 111 gates (RXX 43) | 73.5 exact_len / 71.6 exact_cost (RXX 30.9 / 27.2) | −34% / −36% |
| NISQ (RX/RZ/CZ) | 107 gates (CZ 43) | 160.9 numeric_len / 160.5 numeric_cost (CZ 50.3 / 49.6) | gap remains |

The NISQ qiskit baselines match the paper's (158 vs 149 at L2/L3), so the
comparison is fair. The gap is statistically significant; input composition,
database depth, time budget, and the paper's own random-sampling loop are
ruled out as explanations (report §2.2). It is attributed to database
factorization coverage on 3-wire blocks — see "NISQ levers" below.

## Repository layout

```
matlab_demo/     reference MATLAB demo from the paper (QCOptimDemo)
src/             Python package (`qcr_repro`): gate/token models, compute-graph
                 DB, exact symplectic engine, reducers, QASM I/O
scripts/         benchmarks, figure generation, report builders
results/         benchmark outputs, organized by protocol
  demo_sweep/      paper-style sweep on the demo circuit (strict/, loose/)
  comparison/      head-to-head comparison benchmark (CSV + reports)
  comparison_deep/ deep-database NISQ comparison runs
  dag_compact/     DAG block-compaction comparison benchmark
  rf_sampling/     paper V2/V3 loop reimplementation benchmark
paper/           reference material (citation.bib)
figures/         generated figures (1-6 protocol runs, 7-9 comparison)
report/          results report (LaTeX + Markdown), data
```

## Quickstart

```bash
python -m pip install -e .
# optional baseline compilers (only for scripts that compare against them):
python -m pip install -e ".[baselines]"

python scripts/benchmark_comparison.py --gateset ion_trap --num-circuits 2 --budget 5 --no-baselines   # quick smoke run
```

## Head-to-head comparison with the paper

The flagship benchmark is `scripts/benchmark_comparison.py`. It replicates the
paper's Tables 6/7 protocol (4 qubits, length 300, matching input composition
per gate set), runs our reducers and the paper's baseline compilers on the
same circuits, and writes a comparison report with per-method verdicts:

```bash
python scripts/benchmark_comparison.py --gateset ion_trap --num-circuits 100 --budget 30
python scripts/benchmark_comparison.py --gateset nisq    --num-circuits 100 --budget 30
# deeper NISQ databases (slower one-time build, better NISQ results):
python scripts/benchmark_comparison.py --gateset nisq --deep
```

The first run per gate set builds the lookup databases (cached in the
gitignored `.cache/`): ion trap is quick, NISQ takes about a minute at
default depths, and `--deep` (3-wire depth 6) needs more memory than a
typical laptop provides by default (disk-backed, see below).

Outputs land in `results/comparison/`: per-run CSVs, a Markdown report with
per-type means, WIN/LOSE verdicts against the paper's numbers, and a baseline
fidelity check against the paper's reported qiskit/BQSKit means.

**Timing note.** The "time (s)" column in the comparison reports is the
per-circuit budget cap (each reducer loops until it is exhausted), not a
convergence time, and is not comparable to the paper's Table 2 (a different
task: 100-gate circuits reduced to ~50, ~38 s for their best variant).

Reducers:

- **Exact symplectic engine** (`src/symplectic.py`): Clifford-pool lookups are
  keyed by a bit-exact signed tableau. Replacements are equivalent up to
  global phase by construction, and final verification is exact.
- **Cost-aware objective**: every lookup minimizes (two-qubit count, length),
  the decoherence objective the paper defers to future work, including
  equal-length rewrites that cut RXX/CZ count.
- **Structural passes**: exhaustive window sweeps, single-qubit clustering and
  1-wire collapse, transport shuffling, an RZ-across-CZ pass (NISQ), and
  equivalence-class escape moves with restart-from-best.

## Pre-passes and the batched sweep

Two additions shrink and accelerate the database loop without changing its
results:

- `src/prepass.py` — deterministic pre-passes applied before the DB loop:
  adjacent same-axis rotation fusion (RZ(a)RZ(b) = RZ(a+b), same for
  RX/RY/RXX, kept only when the sum snaps to a pool angle or cancels to the
  identity), plus the cheapest ZX-calculus rules (adjacent CZ·CZ = I and, for
  the diagonal-CZ pools, RZ-gathering across CZ so runs fuse). Every rule uses
  exact angle arithmetic; the output stays pool-representable and the input
  unitary is preserved (`scripts/check_batched_vs_scalar.py`).
- `src/batched.py` — a vectorized version of the exhaustive sweep: all
  window unitaries of a pass are computed with batched matmuls and
  vectorized phase normalization instead of per-window Python matrix
  products. It is bit-identical to the scalar sweep and measures 1.5–1.7×
  faster on the length-300 ion-trap fixpoint; the remaining per-window cost
  is the SHA-256 digest.

Both are opt-in flags of `reduce_circuit` (`algebraic`, `zx`, `use_batched`),
so existing pipelines are unchanged by default.

```bash
python scripts/benchmark_prepass.py --num-circuits 12 --budget 10 --length 300
python scripts/benchmark_prepass.py --gateset nisq --num-circuits 8 --budget 10
```

Outputs land in `results/prepass/`.

## DAG block-compaction (`src/dag.py`)

The exhaustive sweep only finds a reduction if the gates it needs are
physically adjacent in the flat gate list; otherwise it depends on
`transport_shuffle`/`shuffle_commuting_pairs` shuffling them together by
chance. `collect_wire_blocks`/`compact_by_blocks` build the circuit's true
per-wire dependency DAG (a gate depends on the most recent gate touching each
of its wires) and reorder the gate list so every block of gates touching
≤3 wires becomes contiguous, generalizing the 2-qubit block collection used
by production compilers (e.g. Qiskit's `Collect2qBlocks`) to k≤3-wire
blocks. Any two gates that swap position under this reordering act on
disjoint wires, so it is a valid topological order of the dependency DAG and
preserves the circuit's unitary exactly, verified on random ion-trap and
NISQ circuits plus adversarial edge cases (bridging blocks that must split,
single-wire runs, empty circuits) in `scripts/check_dag_compact.py`.

Recompacting on every iteration of the existing stochastic search loop, the
first integration attempt, gave *worse* results at equal budget than the
baseline (ion trap: mean 87.4 vs 83.3 gates, 7 seeds, length 300, 10 s
budget). `compact_by_blocks` is a pure function of the gate list, so
re-running it on an unchanged, already-compacted list is a no-op;
interleaving it with `transport_shuffle` displaced half of the loop's
randomized-diversity passes — the actual mechanism for escaping local
optima — without contributing anything on the passes where nothing had
changed. The fix is to run `compact_by_blocks` only as a deterministic
pre-pass before the stochastic loop starts and leave the loop itself
unchanged (`dag_compact` flag on `reduce_circuit`, `src/reducer.py`).

With that fix, `dag_compact` improves both gate sets at equal budget. The
effect is largest at short budgets and shrinks as the budget grows toward
the paper's own protocol length:

| gate set | budget | delta vs baseline | n |
|---|---|---:|---|
| ion trap (numeric path) | 8 s | −4.5% to −10.0% (4 settings) | 5 seeds × 4 `max_block_len` settings |
| ion trap (numeric path) | 30 s (paper protocol) | −3.6% (68.5 → 66.1) | 20 seeds |
| NISQ | 8 s | −0.2% to −1.7% (4 settings) | 5 seeds × 4 `max_block_len` settings |
| NISQ | 30 s (paper protocol) | −0.1% (151.8 → 151.6, within ±13–14 std) | 20 seeds |

`dag_compact` reaches a given length faster rather than reaching a better
one: given enough time, the baseline's own stochastic search partly catches
up to what the deterministic pre-pass exposes on the first pass. On NISQ
that catch-up is close to complete by 30 s.

Scope: on ion trap, `dag_compact` improves the numeric lookup pipeline
(`numeric_len`); it does not change the stronger exact symplectic-engine
result (`exact_len`/`exact_cost`, 73.5/71.6 gates) that is the headline
ion-trap comparison against the paper's 111. On NISQ, `numeric_len`/
`numeric_cost` is the headline comparison, so the gain there is real but
does not close the gap. The size difference between the two gate sets is
consistent with the existing NISQ diagnosis: `dag_compact` exposes
substantially more candidate windows (mean block length 6–14 gates versus
the sweep's default 8-gate cap, up to 30–54 at length 300), but on NISQ the
lookup database cannot reduce most of them, so exposing more windows barely
moves the final length.

Reproduce: `scripts/check_dag_compact.py` (correctness),
`scripts/diag_dag_maxlen.py` (the `max_block_len` sweep above),
`scripts/benchmark_dag_compact.py` (full baseline/prepass/dag/prepass+dag
comparison with WIN/LOSE verdicts and t-tests, same protocol as
`benchmark_comparison.py`).

## Paper protocol benchmark (numeric, tolerance-based)

- `scripts/benchmark_demo_sweep.py` — depth/iteration/seed sweep on the demo circuit
  (`results/demo_sweep/strict/`, `results/demo_sweep/loose/`)
- `scripts/benchmark_protocol.py` — 100-run Table 6/7-style stats
- `scripts/benchmark_exact.py` — head-to-head numeric vs exact reducers
- `scripts/reduce_demo_circuits.py` — reduce the MATLAB demo QASM files
- `scripts/build_protocol_report.py` — combined strict/loose report tables

## Figures and report

- `scripts/generate_figures.py` regenerates all figures into `figures/`
  (figures 7–9 are the head-to-head comparisons).
- `report/results_report.tex` (and a GitHub-renderable
  `report/results_report.md`) documents the implementation and the
  comparison results.

## NISQ levers (closing the Table 7 gap)

The NISQ gap is real and statistically significant. Input composition,
database depth (`--deep` narrows the mean by only ~2%), time budget, and the
paper's own random-sampling loop (`scripts/benchmark_rf_sampling.py`) are
ruled out as explanations; the residual gap is attributed to database
factorization coverage. Levers implemented and opt-in from
`scripts/benchmark_comparison.py`:

```bash
# 1. larger per-gate-set time budgets (NISQ defaults to 60 s vs 30 s ion trap)
python scripts/benchmark_comparison.py --gateset nisq --budget 60

# 2. V3-style RF-gated lookup: an online RandomForest learns which blocks
#    actually reduce, skipping useless DB lookups and freeing budget for
#    useful ones (needs the optional extra: pip install -e ".[ml]")
python scripts/benchmark_comparison.py --gateset nisq --rf-gate

# 3. exact/numeric hybrid: Clifford-only windows (RX/RZ at +/-pi/2, CZ) are
#    reduced by the bit-exact symplectic engine at deep graph depths; only
#    genuinely non-Clifford (pi/4) windows hit the numeric database
python scripts/benchmark_comparison.py --gateset nisq --hybrid

# 4. disk-backed (SQLite) lookup databases so --deep runs on a laptop
python scripts/benchmark_comparison.py --gateset nisq --deep            # = --backend sqlite
python scripts/benchmark_comparison.py --gateset nisq --backend sqlite --depths 1:12,2:8,3:6,4:4
```

- `--rf-gate` implements the paper's V3 idea (RandomForest-gated lookup): an
  exact memo cache skips already-tried irreducible blocks at zero cost, and a
  lazily-trained classifier (sklearn, optional `[ml]` extra) scores novel
  ones. Skipping never changes equivalence, only the budget allocation. On
  this exhaustive-sweep reducer the classifier is neutral-to-negative for end
  length (`scripts/check_levers.py`), because the sweep already terminates
  when no window reduces, so gating causes premature convergence. The
  paper's V3 gates a random-sampling loop rather than an exhaustive sweep;
  our reimplementation of that loop (`scripts/benchmark_rf_sampling.py`)
  shows the same effect (179.5 vs 168.2), for the same reason.
- `--hybrid` routes every sweep window to the exact symplectic engine when it
  is fully Clifford (all angles ±π/2), falling back to the numeric database
  only for windows containing non-Clifford ±π/4 angles. The Clifford
  sub-pool graphs are small, so they can be built much deeper than the full
  NISQ pool.
- `--backend sqlite` (implied by `--deep`) stores the compute-graph bucket
  tables on disk (SQLite, stdlib) with a bounded in-memory LRU instead of
  all-in-RAM dicts, so deep builds trade memory for disk. Verified
  bit-identical to the RAM backend (`scripts/check_levers.py`).
- The cost-aware (two-qubit-first) objective already applies to CZ via
  `ReductionDatabase.try_reduce_cost` (the `numeric_cost` method), and the
  RZ-across-CZ pass iterates to a fixpoint per invocation.
- `scripts/check_levers.py` asserts every lever preserves the input unitary
  (1e-5) and that the SQLite backend matches RAM results exactly.

**Cross-check against the paper text.** Compute-graph construction, NISQ
angle discretization, CZ token canonicalization, and the 1e-5
replacement-acceptance tolerance all match the paper's description exactly,
ruling those out as sources of the gap. One implementation-specific
hypothesis — that exact hash-based node deduplication (`digest_decimals=10`
in `ComputeGraph`) fragments the compute graph relative to the paper's
tolerance-based (1e-5) node merging — was tested directly
(`scripts/diag_digest_decimals.py`) and ruled out: node counts are identical
from `digest_decimals=10` down to `4` at every depth tested, i.e. no
measurable floating-point fragmentation.

A further lever, not yet benchmarked: every NISQ configuration, including
`--deep`, caps the 4-wire graph at depth 4, even though Table 7's circuits
are 4 qubits wide, so a full-register-width block can only be reduced by the
4-wire graph, at any depth. A depth-5 4-wire graph (38-token pool) has been
built (disk-backed, `1,912,349` nodes, 1342 s): depth 4 reaches 167,449
nodes in about a minute, and depth 5 adds 1,744,900 more. The database
exists but is not yet wired into `NISQ_DEPTHS_DEEP` or re-benchmarked
against the NISQ comparison; that is the next step.
