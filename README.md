# Optimization-Driven Quantum Circuit Reduction — Implementation & Comparison

This repository is a working, reproducible implementation of

**Bodo Rosenhahn, Tobias J Osborne, Christoph Hirche, "Optimization driven quantum circuit reduction," New J. Phys. 27 (2025) 104509**

built around the paper's local term-replacement scheme (variants V1–V3), and it
adds an **original exact engine** that reduces the paper's own test circuits
further than the paper's reported results on the ion-trap gate set, while
staying bit-exact (no 1e-5 tolerance anywhere).

Headline comparison (`results/comparison/`, 4 qubits, 300-gate random
circuits, identical circuits for every method):

| Gate set | Paper "Ours" | Our exact reducer | Difference |
|---|---|---|---|
| Ion trap (RX/RY/RZ/RXX) | 111 gates (RXX 43) | **73.0 gates (RXX 27.4)** | −34% |
| NISQ (RX/RZ/CZ) | 107 gates (CZ 43) | 159 gates (CZ 50) | gap remains |

On NISQ our qiskit baseline replication matches the paper's reported qiskit
numbers (153 vs 149 on L2/L3), confirming the protocol is faithful and the
comparison is fair; the remaining NISQ gap is database-memory limited on this
laptop (see below).

## Repository layout

```
paper_demo/         reference MATLAB demo from the paper (QCOptimDemo)
src/qcr_repro/      Python package: gate/token models, numeric compute-graph DB,
                    exact symplectic engine, reducers, QASM I/O
scripts/            benchmarks, figure generation, report builders
results/            benchmark outputs, organized by protocol
  paper_protocol_1e-5/   paper-style sweep, strict tolerance (1e-5)
  paper_protocol_1e-3/   paper-style sweep, loose tolerance (1e-3)
  baselines/             numeric-reducer protocol runs (paper-style replication)
  costaware_quick/       cost-aware exact engine quick runs
  comparison/            head-to-head comparison benchmark (CSV + reports)
figures/            generated figures (1-6 reproduction, 7-9 comparison)
report/             reproduction report (LaTeX + Markdown), data
future_work/        baseline parity + hardware-metrics probes
```

## Quickstart

```bash
python -m pip install -r requirements.txt   # numpy scipy qiskit bqskit ...
export PYTHONPATH=src
python scripts/benchmark_basic.py            # quick smoke run
```

## Head-to-head comparison with the paper

The flagship benchmark is `scripts/benchmark_comparison.py`. It replicates the
paper's Tables 6/7 protocol (4 qubits, length 300, matching input composition
per gate set), runs our reducers **and** the paper's baseline compilers on the
same circuits, and writes a comparison report with per-method verdicts:

```bash
PYTHONPATH=src python scripts/benchmark_comparison.py --gateset ion_trap --num-circuits 100 --budget 30
PYTHONPATH=src python scripts/benchmark_comparison.py --gateset nisq    --num-circuits 100 --budget 30
# deeper NISQ databases (slower one-time build, better NISQ results):
PYTHONPATH=src python scripts/benchmark_comparison.py --gateset nisq --deep
```

Note: the first run per gate set builds the lookup databases (cached in the
gitignored `.cache/`): ion trap is quick, NISQ takes ~1 minute at default
depths, and `--deep` (3-wire depth 6) needs more memory than a typical laptop.

Outputs land in `results/comparison/`: per-run CSVs, a Markdown report with
per-type means, WIN/LOSE verdicts vs the paper's numbers, and a baseline
fidelity check against the paper's reported qiskit/BQSKit means.

What makes our reducers strong:

- **Exact symplectic engine** (`src/qcr_repro/symplectic.py`): Clifford-pool
  lookups are keyed by a bit-exact signed tableau — replacements are equivalent
  up to global phase *by construction*, and final verification is exact.
- **Cost-aware objective**: every lookup minimizes (two-qubit count, length) —
  the decoherence objective the paper defers to future work — including
  equal-length rewrites that cut RXX/CZ count.
- **Structural passes**: exhaustive window sweeps, single-qubit clustering and
  1-wire collapse, transport shuffling, an RZ-across-CZ pass (NISQ), and
  equivalence-class escape moves with restart-from-best.

## Reproducing the paper protocol (numeric, tolerance-based)

- `scripts/benchmark_sweep.py` — depth/iteration/seed sweep on the demo circuit
  (`results/paper_protocol_1e-5/`, `results/paper_protocol_1e-3/`)
- `scripts/benchmark_paper_protocol.py` — 100-run Table 6/7-style stats
- `scripts/benchmark_exact.py` — head-to-head numeric vs exact reducers
- `scripts/benchmark_basic.py`, `scripts/benchmark_best_of.py` — quick runners
- `scripts/run_matlab_demo_port.py` — reduce the MATLAB demo QASM files
- `scripts/compare_qasm.py` — unitary equivalence check between two QASM files
- `scripts/summarize_benchmarks.py`, `scripts/build_submission_report.py` —
  grouped summaries and combined strict/loose tables

## Figures and report

- `scripts/generate_figures.py` regenerates all figures into `figures/`
  (figures 7–9 are the head-to-head comparisons).
- `report/reproduction_report.tex` (and a GitHub-renderable
  `report/reproduction_report.md`) documents the reproduction and the
  comparison results.

## Future work

`future_work/` contains baseline-parity results against qiskit, a BQSKit
availability probe, and hardware-metric probing (requires authenticated
hardware access, logged as unavailable here).

The main open lever is **NISQ**: the exact engine is Clifford-only, and the
numeric NISQ pipeline still trails the paper (159 vs 107). The `--deep`
databases (3-wire depth 6) and larger per-circuit budgets are expected to
narrow this; the deep graphs need more memory than this laptop allows.
