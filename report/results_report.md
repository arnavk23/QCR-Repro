# Results Report: Optimization-Driven Quantum Circuit Reduction

*Arnav Kapoor*

> This report documents a Python implementation of *Optimization driven
> quantum circuit reduction* (Rosenhahn, Osborne and Hirche, New J. Phys. 27
> (2025) 104509): the paper's MATLAB demo, the reimplementation, benchmark
> infrastructure, figures, and an exact engine built on a bit-exact
> symplectic (Clifford) lookup.

## Headline results

On the paper's ion-trap benchmark (4 qubits, 300-gate random circuits, Table 6):

| Method | Total gates | RXX | vs paper |
|---|---:|---:|---|
| Paper "Ours" | 111 | 43 | — |
| **Our exact (length)** | **73.5** | 30.9 | −37.5 |
| **Our exact (cost-aware)** | **71.6** | **27.2** | −39.4 |
| qiskit L1/L2/L3 | 167.1 / 162.0 / 160.8 | 58.5 / 45.2 / 44.8 | baseline |

On NISQ (Table 7) a gap remains: **160.9 vs 107** (CZ 50.3 vs 43). Our qiskit
baselines match the paper's (158.0 vs 149 at L2/L3), so the head-to-head is
fair. The gap is statistically significant, and input composition, database
depth, time budget and the paper's own random-sampling loop are all ruled out
as explanations; it is isolated to database factorization coverage (§2.2).

## Repository layout

```
matlab_demo/     reference MATLAB demo from the paper (QCOptimDemo)
src/             Python package (`qcr_repro`): gate/token models, numeric
                 compute-graph DB, exact symplectic engine, reducers, QASM I/O
scripts/         benchmarks, figure generation, report builders
results/         benchmark outputs, organized by protocol
  demo_sweep/      paper-style sweep on the demo circuit (strict/ and loose/)
  comparison/      head-to-head comparison benchmark (CSV + reports)
  comparison_deep/ deep-database NISQ comparison runs
  rf_sampling/     paper V2/V3 loop reimplementation benchmark
paper/           reference material (the paper PDF, citation.bib)
figures/         generated figures (1-6 protocol runs, 7-9 comparison)
report/          this report (LaTeX + Markdown) and its data
```

## 1. Paper protocol benchmark

The numeric port (`src/reducer.py`) follows the paper's local
replacement scheme: sample contiguous blocks, remap to local wires, query a
precomputed compute-graph database, accept replacements that are shorter and
unitary-consistent (up to global phase).

Benchmark protocol: depths {3, 4}, iterations {500, 1000, 1500}, seeds {1, 5,
10} — 18 runs per tolerance setting on the paper's 5-qubit demo circuit (300
initial gates).

| depth | iters | best (strict) | mean (strict) | eq rate (strict) | best (loose) | mean (loose) | eq rate (loose) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 500 | 266 | 274.00 | 1.000 | 266 | 274.00 | 1.000 |
| 3 | 1000 | 246 | 250.00 | 0.667 | 246 | 250.00 | 1.000 |
| 3 | 1500 | 240 | 242.67 | 0.667 | 240 | 242.67 | 1.000 |
| 4 | 500 | 262 | 265.33 | 1.000 | 262 | 265.33 | 1.000 |
| 4 | 1000 | 238 | 243.33 | 1.000 | 238 | 243.33 | 1.000 |
| 4 | 1500 | **226** | 231.33 | 0.667 | **226** | 231.33 | 1.000 |

Best strict-valid run: depth 4, 1500 iterations, seed 10 → **228 gates** in
29.8 s. Best loose-valid run: **226 gates** in 2.4 s. The loose metric accepts
more (eq rate 1.000 vs 0.833) and runs ~13× faster — a validation-fidelity
tradeoff the paper's tolerance choice implicitly makes.

## 2. Comparison with the paper

`scripts/benchmark_comparison.py` replicates the paper's Tables 6/7 protocol: 4
qubits, length 300, matching input composition per gate set (ion trap: uniform
over RX/RY/RZ/RXX, matching the paper's input row RX~78, RY~83, RZ~78,
RXX~59; NISQ: CZ-weighted RX:1 RZ:1 CZ:2, matching RX~108, RZ~109, CZ~82).
Every circuit is reduced by our reducers *and* the paper's baseline compilers
(qiskit L1–L3, optional BQSKit), then verified for unitary equivalence.

### 2.1 Ion trap (Table 6) — we reduce further

| method | RX | RY | RZ | RXX | total |
|---|---:|---:|---:|---:|---:|
| paper "Ours" | 10 | 29 | 29 | 43 | **111** |
| exact_len (ours) | 3.8 | 19.2 | 19.5 | 30.9 | **73.5** |
| exact_cost (ours) | 4.0 | 20.2 | 20.1 | **27.2** | **71.6** |
| qiskit L1 | 26.2 | 39.2 | 43.1 | 58.5 | 167.1 |
| qiskit L2 | 37.1 | 36.5 | 43.2 | 45.2 | 162.0 |
| qiskit L3 | 36.8 | 36.1 | 43.1 | 44.8 | 160.8 |

Our exact reducers improve on the paper's "Ours" in **total gates** (−37.5 /
−39.4) and in **two-qubit gates** (30.9 / 27.2 RXX vs 43), with 1.000
equivalence pass rate. The cost-aware objective cuts ~4 more RXX at a similar
gate count — the decoherence objective the paper defers to future work. Both
differences are statistically significant (one-sample t-tests vs the paper's
111, n=100): exact_len t=−26.0 (p=4.8e−46, Cohen's d=−2.6), exact_cost t=−33.2
(p=1.6e−55, d=−3.3); the exact_cost 95% CI is [69.2, 74.0].

**Caveat.** our qiskit baselines (167.1 / 162.0 / 160.8) are *lower*
than the paper's reported qiskit baselines (196 / 204 / 204 on L1/L2/L3). Since
our input composition matches the paper's, the likely explanations are
compiler-version differences (our qiskit 2.5 vs the paper's toolchain) and/or
residual protocol differences; the margin should be read with this in mind.
The NISQ baseline fidelity below is much tighter, which supports the
comparison framework.

### 2.2 NISQ (Table 7) — gap remains

| method | RX | RZ | CZ | total |
|---|---:|---:|---:|---:|
| paper "Ours" | 45 | 19 | 43 | **107** |
| numeric_len (ours) | 68.3 | 42.3 | 50.3 | 160.9 |
| numeric_cost (ours) | 68.6 | 42.3 | 49.6 | 160.5 |
| qiskit L1 | 62.1 | 70.1 | 69.7 | 201.8 |
| qiskit L2 | 63.1 | 44.1 | 50.8 | 158.0 |
| qiskit L3 | 63.0 | 44.1 | 50.8 | 157.9 |

Three observations put this in context:

- **Inputs are composition-matched.** Our generation rule reproduces the
  paper's Table 7 input row (RX/RZ/CZ means 109/109/82 vs the paper's
  108/109/82).
- **Baseline fidelity is good.** Our qiskit L2/L3 means (158.0) closely match
  the paper's reported qiskit baselines (149); L1 matches too (201.8 vs 196).
  The NISQ head-to-head is fair.
- **The gap is real and significant.** One-sample t-tests vs the paper's 107
  give t=+41.8 (p=1.2e−64, d=4.2) for numeric_len and t=+39.7 (p=1.3e−62,
  d=4.0) for numeric_cost; the 95% CI is [157.8, 163.2].

We tested the four most obvious explanations for the gap:

1. **Not input composition.** Inputs match (above).
2. **Not database depth.** Our default database already matches the paper's
   3-wire depth-5 tool; `--deep` (3-wire depth 6, disk-backed SQLite) narrows
   the mean by only ~2% (160.9 → 157.3), still far from the paper's 107
   (t=+39.3, 95% CI [154.7, 159.7]).
3. **Not the time budget.** The exhaustive sweep reaches the same end length
   at 10 s and at 60 s, and best-of-K restarts barely move the mean.
4. **Not the search loop.** Our reimplementation of the paper's V2 (random
   sub-block sampling) reaches only ~168 at equal budget, and the V3 RF-gated
   variant is *worse* (~179.5), because gating suppresses reducible lookups
   and its decision overhead cuts the iteration count by an order of
   magnitude (`scripts/benchmark_rf_sampling.py`):

| loop variant (10 s budget) | end (mean) | replacements | iterations |
|---|---:|---:|---:|
| exhaustive sweep | 161.5 | — | — |
| V2 random sampling | 168.2 | 78.2 | 44,166 |
| V3 RF-gated sampling | 179.5 | 69.0 | 3,287 |

The residual gap is therefore attributed to **database factorization
coverage**: the paper's compute graph supports block-shortening rewrites ours
does not find, concentrated in single-qubit rotations (RZ 42.3 vs 19, RX 68.3
vs 45) and two-qubit gates (CZ 50.3 vs 43). This is the falsifiable target of
the next experiment.

**Timing caveat (applies to both tables):** the reducers run under a fixed
per-circuit **time budget** (30 s ion-trap, 60 s NISQ in the committed runs),
so the runtime column in `results/comparison/*` is a cutoff, not a convergence
time, and is not comparable to the paper's Table 2 (a different task — reducing
100-gate circuits to ~50, ~38 s for their best variant). With a larger budget
the NISQ numbers above improve; see the README's "NISQ levers" section
(RF-gated lookup, exact/numeric hybrid, disk-backed `--deep` databases).

## 3. Figures

All figures are regenerated by `scripts/generate_figures.py` into `figures/`:

| Figure | Content |
|---|---|
| 1 | Circuit-length motivation (original / MATLAB / Python strict / loose) |
| 2 | Compute-graph growth with depth (values from the paper, Table 1) |
| 3 | Gate composition: input vs reduced circuits |
| 4 | Reduction vs iteration budget |
| 5 | Runtime vs reduced length |
| 6 | Seed variance boxplot (strict) |
| 7 | **Comparison: ion trap (73.5 vs 111)** |
| 8 | **Comparison: NISQ (161 vs 107)** |
| 9 | **Two-qubit counts: hardware-cost objective** |

![Figure 1](../figures/figure1_motivation.png)
![Figure 2](../figures/figure2_compute_graph_growth.png)
![Figure 3](../figures/figure3_pipeline.png)
![Figure 4](../figures/figure4_reduction_curve.png)
![Figure 5](../figures/figure5_runtime_vs_length.png)
![Figure 6](../figures/figure6_boxplot.png)
![Figure 7](../figures/figure7_comparison_ion_trap.png)
![Figure 8](../figures/figure8_comparison_nisq.png)
![Figure 9](../figures/figure9_two_qubit_counts.png)

## 4. Open directions

The main open lever is **NISQ**: the numeric NISQ pipeline trails the paper
(160.9 vs 107). Input composition, database depth, time budget and the search
loop are ruled out; the open question is database factorization coverage on
3-wire blocks. Four levers remain implemented and opt-in from
`scripts/benchmark_comparison.py` (see README "NISQ levers"): the paper's V3
RF-gated lookup (`--rf-gate`), an exact/numeric hybrid that solves Clifford
windows bit-exactly (`--hybrid`), disk-backed SQLite databases so `--deep`
runs on a laptop (`--backend sqlite`), and larger per-gate-set budgets (NISQ defaults to 60 s).

## 5. Pre-passes and the batched sweep

Two additions to the numeric pipeline, both verified bit-exact or
unitary-preserving in `scripts/check_batched_vs_scalar.py`:

- **Algebraic / ZX pre-passes** (`src/prepass.py`): deterministic,
  search-free reductions applied before the DB loop — adjacent same-axis
  rotation fusion (kept only when the sum snaps to a pool angle or cancels),
  CZ·CZ cancellation, and (NISQ) RZ-gathering across CZ so runs fuse. On the
  Table 7-style NISQ inputs this removes ~25% of the gates before any search;
  on the ion-trap inputs (pool of only ±π/2) it removes the exact
  cancellations. Output stays pool-representable; input unitary preserved.
- **Batched sweep** (`src/batched.py`): all window unitaries of a
  pass are computed with batched numpy matmuls and vectorized phase
  normalization. It is bit-identical to the scalar sweep and measures
  ~1.5-1.7x faster on the length-300 ion-trap fixpoint; per-window SHA-256
  digests remain the floor.

`scripts/benchmark_prepass.py` compares baseline / prepass / prepass+batched
under a fixed per-circuit budget (`results/prepass/`). At equal budget the
prepass+batched pipeline ends with the same-or-better final lengths
(ion-trap 89.5 vs 90.8 baseline; NISQ 159.2 vs 161.7) while spending less
wall-clock per sweep; the ion-trap WIN over the paper's 111 gates is
preserved. The batched sweep does not change which circuits are found — it
only finds them faster, and the pre-passes make the input cheaper for the
expensive loop.

## 6. DAG block-compaction

`src/dag.py` adds a deterministic alternative to `transport_shuffle` for
exposing reducible windows to the exhaustive sweep. It builds the circuit's
true per-wire dependency DAG (a gate depends on the most recent gate
touching each of its wires) and reorders the gate list so every block of
gates touching ≤3 wires becomes contiguous — generalizing the 2-qubit block
collection used by production compilers (e.g. Qiskit's `Collect2qBlocks`) to
k≤3-wire blocks. Any two gates that swap position under this reordering act
on disjoint wires, so it is a valid topological order of the dependency DAG
and provably preserves the unitary; verified on random ion-trap/NISQ
circuits and adversarial edge cases in `scripts/check_dag_compact.py`.
Resulting blocks run 6-14 gates on average (up to 30-54 at length 300) —
well past the sweep's default 8-gate window cap, so the sweep sees
substantially more candidate windows per pass than the un-compacted circuit.

**Integration pitfall (worth recording).** The first attempt re-ran
`compact_by_blocks` on every iteration of the existing stochastic loop
alongside `transport_shuffle`. This made results *worse* at equal budget
(ion trap: baseline 83.3 vs dag 87.4 mean length over 7 seeds, len 300,
10s budget). Cause: `compact_by_blocks` is a pure function of the gate list,
so re-running it on an unchanged, already-compacted list is a no-op;
alternating it with `transport_shuffle` therefore quietly removed half of
the loop's randomized-diversity passes (the mechanism that actually escapes
local optima) on the — common — passes where nothing had changed since the
last compaction, without adding anything in exchange. Fix: apply
`compact_by_blocks` only as a deterministic pre-pass before the stochastic
loop starts (`dag_compact` flag on `reduce_circuit`), and leave the loop
itself untouched.

At a short 8s budget, `dag_compact` is a consistent win in both gate sets
(`scripts/diag_dag_maxlen.py`, 5 seeds × 4 `max_block_len` settings, length
300). At the paper's own protocol budget (30s, 20 seeds,
`scripts/benchmark_dag_compact.py`) the ion-trap win persists but shrinks;
the NISQ win — already an order of magnitude smaller — washes out to noise:

| gate set | budget | delta vs baseline | n |
|---|---|---:|---|
| ion trap (numeric path) | 8s | -7.6% mean (-4.5% to -10.0% range) | 5 seeds × 4 settings |
| ion trap (numeric path) | 30s | **-3.6%** (68.5 → 66.1 mean) | 20 seeds |
| NISQ | 8s | -1.1% mean (-0.2% to -1.7% range) | 5 seeds × 4 settings |
| NISQ | 30s | **-0.1%** (151.8 → 151.6 mean, std ±13-14) | 20 seeds |

That the win shrinks as budget grows is itself informative: `dag_compact`'s
value is *reaching a given length faster*, not a better asymptote — given
enough time, the baseline's own stochastic search (`transport_shuffle` +
escape) partly catches up to what the deterministic pre-pass exposes on the
first pass. On NISQ that catch-up is close to total by 30s.

On **ion trap**, this improves the *numeric* lookup pipeline (`numeric_len`)
— it is a different, weaker pipeline than the *exact* symplectic engine
already reported in §2.1 (`exact_len`/`exact_cost`, 73.5/71.6 gates), which
remains the strongest ion-trap result and is unaffected by this change. On
**NISQ**, `numeric_len`/`numeric_cost` *is* the headline comparison against
the paper (§2.2), so even the larger short-budget gain there was always
modest — and at the standard 30s budget it is not distinguishable from
baseline. That `dag_compact` finds far more candidate windows on NISQ too,
yet the final length barely moves regardless of budget, is independent
evidence for this report's existing diagnosis (§2.2): the NISQ gap is a
**database factorization coverage** problem, not a window-discovery problem.
Exposing more windows to a lookup database that already can't reduce most of
them does not help; the open question remains what the paper's compute graph
covers that ours does not.

## Setup and commands

```bash
pip install -e .
export PYTHONPATH=src
python scripts/generate_figures.py                       # regenerate figures/
python scripts/benchmark_comparison.py --gateset ion_trap --num-circuits 100 --budget 30
python scripts/benchmark_comparison.py --gateset nisq    --num-circuits 100 --budget 30
python scripts/build_protocol_report.py                # combined strict/loose table
```

All cited numbers live under `results/` (per-protocol CSVs and reports).
