# Optimization-Driven Quantum Circuit Reduction — Implementation & Comparison

This is a Python implementation of

**Bodo Rosenhahn, Tobias J Osborne, Christoph Hirche, "Optimization driven quantum circuit reduction," New J. Phys. 27 (2025) 104509**

implementing the paper's local term-replacement scheme (variants V1–V3) and
adding a bit-exact symplectic (Clifford) engine that reduces the paper's
ion-trap test circuits further than the paper reports (no tolerance anywhere).

Headline comparison (`results/comparison/`, 4 qubits, 300-gate random
circuits, identical circuits for every method, 100 circuits per method):

| Gate set | Paper "Ours" | Our reducers | Difference |
|---|---|---|---|
| Ion trap (RX/RY/RZ/RXX) | 111 gates (RXX 43) | **73.5 exact_len / 71.6 exact_cost (RXX 30.9 / 27.2)** | −34% / −36% |
| NISQ (RX/RZ/CZ) | 107 gates (CZ 43) | 160.9 numeric_len / 160.5 numeric_cost (CZ 50.3 / 49.6) | gap remains |

On NISQ our qiskit baselines match the paper's (158 vs 149 at L2/L3), so the
comparison is fair; the remaining gap is statistically significant, and input
composition, database depth, time budget and the paper's own random-sampling
loop are all ruled out as explanations (report §2.2). The gap is isolated to
database factorization coverage on 3-wire blocks; the "NISQ levers" section
below documents the remaining tools.

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
  rf_sampling/     paper V2/V3 loop reimplementation benchmark
paper/           reference material (the paper PDF, citation.bib)
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
per gate set), runs our reducers **and** the paper's baseline compilers on the
same circuits, and writes a comparison report with per-method verdicts:

```bash
python scripts/benchmark_comparison.py --gateset ion_trap --num-circuits 100 --budget 30
python scripts/benchmark_comparison.py --gateset nisq    --num-circuits 100 --budget 30
# deeper NISQ databases (slower one-time build, better NISQ results):
python scripts/benchmark_comparison.py --gateset nisq --deep
```

Note: the first run per gate set builds the lookup databases (cached in the
gitignored `.cache/`): ion trap is quick, NISQ takes ~1 minute at default
depths, and `--deep` (3-wire depth 6) needs more memory than a typical laptop.

Outputs land in `results/comparison/`: per-run CSVs, a Markdown report with
per-type means, WIN/LOSE verdicts vs the paper's numbers, and a baseline
fidelity check against the paper's reported qiskit/BQSKit means.

**Timing caveat.** The "time (s)" column in the comparison reports is the
per-circuit **budget cap** (each reducer loops until it is exhausted), not a
convergence time, and is not comparable to the paper's Table 2 (a different
task: 100-gate circuits reduced to ~50, ~38 s for their best variant).

Reducers:

- **Exact symplectic engine** (`src/symplectic.py`): Clifford-pool lookups are
  keyed by a bit-exact signed tableau — replacements are equivalent up to
  global phase by construction, and final verification is exact.
- **Cost-aware objective**: every lookup minimizes (two-qubit count, length) —
  the decoherence objective the paper defers to future work — including
  equal-length rewrites that cut RXX/CZ count.
- **Structural passes**: exhaustive window sweeps, single-qubit clustering and
  1-wire collapse, transport shuffling, an RZ-across-CZ pass (NISQ), and
  equivalence-class escape moves with restart-from-best.

## Pre-passes and the batched sweep

Two cheap optimizations shrink and accelerate the database loop without
changing its results:

- `src/prepass.py` — deterministic pre-passes applied *before* the
  DB loop: adjacent same-axis rotation fusion (RZ(a)RZ(b) = RZ(a+b), same for
  RX/RY/RXX, kept only when the sum snaps to a pool angle or cancels to the
  identity) plus the cheapest ZX-calculus rules (adjacent CZ·CZ = I and, for
  the diagonal-CZ pools, RZ-gathering across CZ so runs fuse). Every rule is
  exact angle arithmetic, the output stays pool-representable, and the input
  unitary is preserved (verified in
  `scripts/check_batched_vs_scalar.py`).
- `src/batched.py` — a vectorized version of the exhaustive sweep:
  all window unitaries of a pass are computed with batched matmuls and
  vectorized phase normalization instead of per-window Python matrix products.
  It is *bit-identical* to the scalar sweep (same protocol, same results) and
  measures ~1.5-1.7x faster on the length-300 ion-trap fixpoint; the remaining
  per-window cost is the SHA-256 digest.

Both are opt-in flags of `reduce_circuit` (`algebraic`, `zx`, `use_batched`),
so existing pipelines are unchanged by default. Run the head-to-head
comparison (baseline vs prepass vs prepass+batched, plus the paper's numbers):

    python scripts/benchmark_prepass.py --num-circuits 12 --budget 10 --length 300
    python scripts/benchmark_prepass.py --gateset nisq --num-circuits 8 --budget 10

Outputs land in `results/prepass/` (`comparison_prepass_report.md/csv/json`).

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

The NISQ gap is real and statistically significant, but experiments rule out
input composition, database depth (a `--deep` run narrows the mean by only
~2%), time budget, and the paper's own random-sampling loop
(`scripts/benchmark_rf_sampling.py`) as explanations; the residual gap is
attributed to database factorization coverage. Four levers remain implemented
and opt-in from `scripts/benchmark_comparison.py`:

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
  lazily-trained classifier (sklearn, optional `[ml]` extra) scores novel ones.
  Skipping never changes equivalence, only the budget allocation. **Measured
  caveat:** on this exhaustive-sweep reducer the classifier is
  neutral-to-negative for end length (`scripts/check_levers.py`, tuning
  probes) — the sweep already terminates when no window reduces, so gating
  causes premature convergence. The paper's V3 gated a *random-sampling* loop;
  we benchmarked our reimplementation of that loop too
  (`scripts/benchmark_rf_sampling.py`) and the classifier is worse there as
  well (179.5 vs 168.2), for the same premature-convergence reason.
- `--hybrid` routes every sweep window to the exact symplectic engine when it
is fully Clifford (all angles ±π/2), falling back to the numeric database
only for windows containing non-Clifford ±π/4 angles. The Clifford sub-pool
graphs are small, so they can be built much deeper than the full NISQ pool.
- `--backend sqlite` (implied by `--deep`) stores the compute-graph bucket
tables on disk (SQLite, stdlib) with a bounded in-memory LRU instead of
all-in-RAM dicts, so deep builds trade memory for disk. Verified
bit-identical to the RAM backend
(`scripts/check_levers.py`).
- The cost-aware (two-qubit-first) objective already applies to CZ via
`ReductionDatabase.try_reduce_cost` (the `numeric_cost` method), and the
RZ-across-CZ pass now iterates to a fixpoint per invocation.
- `scripts/check_levers.py` asserts every lever
preserves the input unitary (1e-5) and that the SQLite backend matches
RAM results exactly.
