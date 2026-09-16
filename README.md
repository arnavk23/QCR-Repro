# Optimization-Driven Quantum Circuit Reduction

An exact, tolerance-free compute-graph engine for local term-replacement
quantum circuit reduction, evaluated against Rosenhahn, Osborne & Hirche,
*"Optimization driven quantum circuit reduction,"* New J. Phys. **27**,
104509 (2025) [[doi]](https://doi.org/10.1088/1367-2630/ae0e40).

On the ion-trap (all-Clifford) gate set, a bit-exact signed-symplectic
engine plus a two-qubit-aware search objective beats the published result,
and also beats BQSKit L2/L3 by more than 2x on both total and two-qubit
gate count (`results/comparison_bqskit_ion/`, n=30). On the NISQ gate set
a gap remains against the published baseline itself; six candidate
explanations for it are tested and ruled out in a systematic diagnostic
study, and a seventh (a node-keying fragility fix, kept for its own sake)
and an eighth (a compute-graph reach extension, net negative once its
overhead is accounted for) are added and ruled out the same way
(`report/draft_paper.tex`). A ninth attempt, directly porting
`matlab_demo/QCOptimDemo`'s own reduction loop, *does* help: run as the
entire search budget ahead of this codebase's exhaustive sweep (the
default as of this pass), it improves NISQ's mean gate count by ~3%
(p < 10⁻¹⁶, n = 100) over the exhaustive sweep alone — a real, modest win
layered on top of the still-open gap, not a fix for it
(`src/reducer.py`'s `reduce_matlab_burst` / `burst_frac`).

The repository also includes extensions not in the original paper: a
dependency-graph based block-reordering pass (`src/dag.py`), a
disk-backed compute-graph backend for databases too large to fit in RAM,
an incremental per-window search construction that removes redundant
from-scratch rebuilds in the sweep's hot loop (`src/reducer.py`,
`src/exact_reducer.py`), the matlab_demo-faithful burst reducer above, and
a richer mixed Clifford/non-Clifford ion-trap pool (±π/2, ±π/4, ±π/8) that
stress-tests the whole approach on a larger gate set. That last one is an
honest negative result worth reading before assuming this scales freely:
the larger pool forces much shallower compute graphs under the same build
budget, and at that depth the method currently loses even to plain qiskit
(Section "A mixed Clifford/non-Clifford ion-trap pool: a boundary case" in
`report/draft_paper.tex`) — a real boundary condition of the technique,
not a bug.

## Results

4 qubits, length-300 circuits, 100 circuits per gate set, identical inputs
across methods (`results/comparison/`):

| Gate set | Published "Ours" | This work | Δ |
|---|---:|---:|---:|
| Ion trap (RX/RY/RZ/RXX) | 111 gates (43 RXX) | **66.7** gates (25.9 RXX) | −40% |
| NISQ (RX/RZ/CZ) | 107 gates (43 CZ) | 155.1 gates (48.9 CZ) | gap remains |

Both differences vs. the published baseline are statistically significant
(p < 10⁻⁴⁵ and p < 10⁻⁶⁰ respectively, one-sample t-test, n = 100). The NISQ
row already reflects the matlab_demo-faithful burst reducer (`burst_frac=1.0`,
now the default for this gate set); without it the same protocol (identical
seeds, identical hybrid setting) gives 160.5 gates — the burst mode's own
paired improvement over the exhaustive sweep is independently significant
too (−3.3%, p = 9×10⁻²², paired t-test, n = 100). Full method, protocol, and
the NISQ diagnostic study are in `report/draft_paper.tex`.

**Against BQSKit** (n=30 per gate set, `results/comparison_bqskit_ion/`,
`results/comparison_bqskit_nisq30/`):

| Gate set | This work (total, twq) | BQSKit L2 | BQSKit L3 |
|---|---:|---:|---:|
| Ion trap | 67.8, 25.4 | 155.1, 36.6 | 133.0, 30.3 |
| NISQ (burst default) | 156.1, 50.8 | 227.5, 65.2 | 211.4, **50.3** |

Ion-trap is a clean win on both metrics. NISQ wins on total gate count but
two-qubit count is a genuine near-tie (BQSKit L3 edges us slightly there,
consistently across two independent runs, now with the burst default
included too) — a real, currently-open limitation, not noise; see the NISQ
cost-aware diagnostic in `report/draft_paper.tex` for why. In this run
BQSKit L3 itself also failed unitary verification on 1/30 circuits
(pass rate 0.967) — an issue on BQSKit's side of the comparison, not ours,
noted here rather than silently dropped. BQSKit L4 is omitted: every circuit at
that level failed unitary verification in our setup (likely an internal
qubit-relabeling issue at that optimization level), so its numbers aren't
trustworthy as reported.

## Installation

```bash
python -m pip install -e .
python -m pip install -e ".[baselines]"   # qiskit / BQSKit, for baseline comparisons
python -m pip install -e ".[ml]"          # scikit-learn, for the RF-gated lookup lever
```

Requires Python ≥ 3.10.

## Usage

```bash
# smoke test
python scripts/benchmark_comparison.py --gateset ion_trap --num-circuits 2 --budget 5 --no-baselines

python scripts/benchmark_comparison.py --gateset ion_trap --num-circuits 100 --budget 30 --fast-sweep
python scripts/benchmark_comparison.py --gateset nisq --num-circuits 100 --budget 60

# regenerate all figures
python scripts/generate_figures.py
```

Outputs land in `results/comparison/`: per-circuit CSVs, a Markdown report
with per-type means, and verdicts against the published numbers.

## Repository structure

```
src/             qcr_repro package: gate/token models, compute-graph
                 database (RAM + disk-backed), exact symplectic engine,
                 reducers, QASM I/O
scripts/         benchmarks, verification checks, figure generation
results/         benchmark outputs by protocol
figures/         generated figures
report/          draft_paper.tex — method, results, NISQ diagnostic study
matlab_demo/     reference MATLAB demo from the original paper
```

## Citation

```bibtex
@article{Rosenhahn2025Optimization,
  author  = {Rosenhahn, Bodo and Osborne, Tobias J and Hirche, Christoph},
  title   = {Optimization driven quantum circuit reduction},
  journal = {New Journal of Physics},
  volume  = {27},
  number  = {10},
  pages   = {104509},
  year    = {2025},
  doi     = {10.1088/1367-2630/ae0e40}
}
```

## License

MIT — see [LICENSE](LICENSE).
