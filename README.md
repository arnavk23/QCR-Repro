# Optimization-Driven Quantum Circuit Reduction

An exact, tolerance-free compute-graph engine for local term-replacement
quantum circuit reduction, evaluated against Rosenhahn, Osborne & Hirche,
*"Optimization driven quantum circuit reduction,"* New J. Phys. **27**,
104509 (2025) [[doi]](https://doi.org/10.1088/1367-2630/ae0e40).

On the ion-trap (all-Clifford) gate set, a bit-exact signed-symplectic
engine plus a two-qubit-aware search objective beats the published result,
and also beats BQSKit L2/L3 by more than 2x on both total and two-qubit
gate count (`results/comparison_bqskit_ion/`, n=30). On NISQ, the same
method wins on total gate count against BQSKit too, though two-qubit count
there is a genuine near-tie rather than a clean win
(`results/comparison_bqskit_nisq30/`, n=30). On the NISQ gate set a gap
remains against the published baseline itself; six candidate explanations
for it are tested and ruled out in a systematic diagnostic study
(`report/draft_paper.tex`).

The repository also includes extensions not in the original paper: a
dependency-graph based block-reordering pass (`src/dag.py`), a
disk-backed compute-graph backend for databases too large to fit in RAM,
an incremental per-window search construction that removes redundant
from-scratch rebuilds in the sweep's hot loop (`src/reducer.py`,
`src/exact_reducer.py`), and a richer mixed Clifford/non-Clifford
ion-trap pool (±π/2, ±π/4, ±π/8) that stress-tests the whole approach on a
larger gate set. That last one is an honest negative result worth reading
before assuming this scales freely: the larger pool forces much shallower
compute graphs under the same build budget, and at that depth the method
currently loses even to plain qiskit (Section "A mixed
Clifford/non-Clifford ion-trap pool: a boundary case" in
`report/draft_paper.tex`) — a real boundary condition of the technique,
not a bug.

## Results

4 qubits, length-300 circuits, 100 circuits per gate set, identical inputs
across methods (`results/comparison/`):

| Gate set | Published "Ours" | This work | Δ |
|---|---:|---:|---:|
| Ion trap (RX/RY/RZ/RXX) | 111 gates (43 RXX) | **66.7** gates (25.9 RXX) | −40% |
| NISQ (RX/RZ/CZ) | 107 gates (43 CZ) | 160.5 gates (49.6 CZ) | gap remains |

Both differences are statistically significant (p < 10⁻⁴⁵ and p < 10⁻⁶⁰
respectively, one-sample t-test, n = 100). Full method, protocol, and the
NISQ diagnostic study are in `report/draft_paper.tex`.

**Against BQSKit** (n=30 per gate set, `results/comparison_bqskit_ion/`,
`results/comparison_bqskit_nisq30/`):

| Gate set | This work (total, twq) | BQSKit L2 | BQSKit L3 |
|---|---:|---:|---:|
| Ion trap | 67.8, 25.4 | 155.1, 36.6 | 133.0, 30.3 |
| NISQ | 160.5, 51.2 | 228.9, 65.2 | 213.2, **50.8** |

Ion-trap is a clean win on both metrics. NISQ wins on total gate count but
two-qubit count is a genuine near-tie (BQSKit L3 edges us slightly there,
consistently across two independent runs) — a real, currently-open
limitation, not noise; see the NISQ cost-aware diagnostic in
`report/draft_paper.tex` for why. BQSKit L4 is omitted: every circuit at
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
