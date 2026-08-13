# Optimization-Driven Quantum Circuit Reduction — Implementation & Comparison

This repository is a working implementation of

**Bodo Rosenhahn, Tobias J Osborne, Christoph Hirche, "Optimization driven quantum circuit reduction," New J. Phys. 27 (2025) 104509**

built around the paper's local term-replacement scheme (variants V1–V3), and it
adds an **original exact engine** that reduces the paper's own test circuits
further than the paper's reported results on the ion-trap gate set, while
staying bit-exact (no 1e-5 tolerance anywhere).

Headline comparison (`results/comparison/`, 4 qubits, 300-gate random
circuits, identical circuits for every method, 100 circuits per method):

| Gate set | Paper "Ours" | Our reducers | Difference |
|---|---|---|---|
| Ion trap (RX/RY/RZ/RXX) | 111 gates (RXX 43) | **73.5 exact_len / 71.6 exact_cost (RXX 30.9 / 27.2)** | −34% / −36% |
| NISQ (RX/RZ/CZ) | 107 gates (CZ 43) | 160.9 numeric_len / 160.5 numeric_cost (CZ 50.3 / 49.6) | gap remains |

On NISQ our qiskit baseline replication matches the paper's reported qiskit
numbers (158 vs 149 on L2/L3), confirming the protocol is faithful and the
comparison is fair; the remaining NISQ gap is database-depth limited on this
laptop. The NISQ levers below (RF-gated lookup, exact/numeric hybrid,
disk-backed databases, larger honest budgets) are the path to closing it.

## Repository layout

```
matlab_demo/         reference MATLAB demo from the paper (QCOptimDemo)
src/                Python package (`qcr_repro`): gate/token models
                    compute-graph DB, exact symplectic engine, reducers, QASM I/O
scripts/            benchmarks, figure generation, report builders
results/            benchmark outputs, organized by protocol
  demo_sweep/   paper-style sweep on the demo circuit
    strict/             strict tolerance (1e-5)
    loose/              loose tolerance (1e-3)
  comparison/            head-to-head comparison benchmark (CSV + reports)
  paper_tables/            paper tables from benchmark data (scripts/generate_paper_tables.py)
paper/              reference material (the paper PDF)
figures/            generated figures (1-6 protocol runs, 7-9 comparison)
report/             results report (LaTeX + Markdown), data
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
per-circuit **budget cap** (each reducer loops until its budget is exhausted)
— it is a cutoff, not a convergence time, so it must not be read as an
apple-to-apples timing comparison against the paper's Table 2 (which measures
a different task: reducing 100-gate circuits to ~50, ~38 s for their best
variant).

What makes our reducers strong:

- **Exact symplectic engine** (`src/symplectic.py`): Clifford-pool
  lookups are keyed by a bit-exact signed tableau — replacements are equivalent
  up to global phase *by construction*, and final verification is exact.
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

The NISQ gap is a search-space problem, not an algorithmic one: the exact
engine is Clifford-only, and the numeric fallback is limited by how deep a
lookup database a laptop can build in memory. Four levers are implemented and
opt-in from `scripts/benchmark_comparison.py`:

```bash
# 1. larger, honestly-reported budgets (NISQ defaults to 60 s vs 30 s ion trap)
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
exact memo cache remembers already-tried irreducible blocks (a guaranteed
win — repeated windows after transport shuffles are skipped at zero cost),
and a lazily-trained classifier (sklearn `RandomForestClassifier`, optional
`[ml]` extra) scores novel blocks. Skipping never changes equivalence — it
only reallocates budget. **Measured caveat:** on this exhaustive-sweep
reducer the classifier component is neutral-to-negative for end length
(tested at 10-60 s budgets, `scripts/check_levers.py`
and the tuning probes) — the sweep already terminates when no window
reduces, so gating causes premature convergence. The paper's V3 gated their
*random-sampling* loop, where DB lookups are the bottleneck; that is the
regime where the classifier helps.
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
