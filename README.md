# Optimization-Driven Quantum Circuit Reduction

An exact, tolerance-free compute-graph engine for local term-replacement
quantum circuit reduction, evaluated against Rosenhahn, Osborne & Hirche,
*"Optimization driven quantum circuit reduction,"* New J. Phys. **27**,
104509 (2025) [[doi]](https://doi.org/10.1088/1367-2630/ae0e40).

On the ion-trap (all-Clifford) gate set, a bit-exact signed-symplectic
engine plus a two-qubit-aware search objective beats the published result.
On the NISQ gate set a gap remains; six candidate explanations for it are
tested and ruled out in a systematic diagnostic study (`report/draft_paper.tex`).
The repository also includes two extensions not in the original paper: a
dependency-graph based block-reordering pass (`src/dag.py`) and a
disk-backed compute-graph backend for databases too large to fit in RAM.

## Results

4 qubits, length-300 circuits, 100 circuits per gate set, identical inputs
across methods (`results/comparison/`):

| Gate set | Published "Ours" | This work | Δ |
|---|---:|---:|---:|
| Ion trap (RX/RY/RZ/RXX) | 111 gates (43 RXX) | **71.6** gates (27.2 RXX) | −36% |
| NISQ (RX/RZ/CZ) | 107 gates (43 CZ) | 160.5 gates (49.6 CZ) | gap remains |

Both differences are statistically significant (p < 10⁻⁴⁵ and p < 10⁻⁶⁰
respectively, one-sample t-test, n = 100). Full method, protocol, and the
NISQ diagnostic study are in `report/draft_paper.tex`.

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

# full protocol (builds and caches lookup databases on first run)
python scripts/benchmark_comparison.py --gateset ion_trap --num-circuits 100 --budget 30
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
