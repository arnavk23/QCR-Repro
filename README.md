# Optimization Driven Quantum Circuit Reduction — Implementation Workspace

This repository contains a working, reproducible implementation workflow for the paper:

**Bodo Rosenhahn et al., “Optimization driven quantum circuit reduction,” New J. Phys. 27 (2025) 104509**

The project now includes:
- extracted reference MATLAB demo code,
- Python reimplementation and experiment scripts,
- benchmark outputs (strict and loose validation settings),
- generated figures,
- an academic LaTeX report and compiled PDF,
- a concrete `future_work` folder for baseline parity and hardware-metric probing.

## 1) Reference code from paper [38]

The provided archive is extracted at:
- [paper_code/QCOptimDemo](paper_code/QCOptimDemo)

Key files used for mapping behavior:
- [paper_code/QCOptimDemo/optimCodeGMode1DComp.m](paper_code/QCOptimDemo/optimCodeGMode1DComp.m)
- [paper_code/QCOptimDemo/FullQCGraphIOn3_4.m](paper_code/QCOptimDemo/FullQCGraphIOn3_4.m)
- [paper_code/QCOptimDemo/QoperatorsIon1.m](paper_code/QCOptimDemo/QoperatorsIon1.m)

## 2) Python implementation

Core package:
- [src/qcr_repro](src/qcr_repro)

Implemented modules include:
- gate/unitary construction,
- tokenization and operator pools,
- global-phase-aware equivalence checks,
- compute-graph generation,
- MATLAB-demo-compatible QASM I/O,
- local replacement reducer.

## 3) Main scripts

- [scripts/run_matlab_demo_port.py](scripts/run_matlab_demo_port.py): run reducer on MATLAB demo QASM.
- [scripts/compare_qasm.py](scripts/compare_qasm.py): compare two QASM files up to global phase.
- [scripts/benchmark_reducer.py](scripts/benchmark_reducer.py): sweep depth/iterations/seeds.
- [scripts/summarize_benchmarks.py](scripts/summarize_benchmarks.py): grouped paper-style summary tables.
- [scripts/build_submission_report.py](scripts/build_submission_report.py): consolidated submission tables.
- [scripts/generate_paper_style_figures.py](scripts/generate_paper_style_figures.py): generate implementation figure set.
- [scripts/future_work_baselines.py](scripts/future_work_baselines.py): baseline parity + hardware probe artifacts.

## 4) Environment status

Python + virtual environment are already set up in this workspace:
- [.venv](.venv)

Use:

```powershell
$env:PYTHONPATH='src'
& .\.venv\Scripts\python.exe --version
```

If reinstalling dependencies is needed:

```powershell
$env:PYTHONPATH='src'
& .\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

## 5) Reproduced benchmark outputs

Strict tolerance sweep (`atol=1e-5`):
- [results/benchmark_reducer.csv](results/benchmark_reducer.csv)
- [results/benchmark_summary.txt](results/benchmark_summary.txt)

Loose tolerance sweep (`atol=1e-3`):
- [results_tol1e3/benchmark_reducer.csv](results_tol1e3/benchmark_reducer.csv)
- [results_tol1e3/benchmark_summary.txt](results_tol1e3/benchmark_summary.txt)

Combined submission-style tables:
- [results/submission_report_table.csv](results/submission_report_table.csv)
- [results/submission_report.md](results/submission_report.md)

## 6) Generated figures

Implementation figure set:
- [implementation/figures](implementation/figures)

Includes:
- [implementation/figures/figure1_motivation.png](implementation/figures/figure1_motivation.png)
- [implementation/figures/figure2_compute_graph_growth.png](implementation/figures/figure2_compute_graph_growth.png)
- [implementation/figures/figure3_pipeline.png](implementation/figures/figure3_pipeline.png)
- [implementation/figures/figure4_reduction_curve.png](implementation/figures/figure4_reduction_curve.png)
- [implementation/figures/figure5_runtime_vs_length.png](implementation/figures/figure5_runtime_vs_length.png)
- [implementation/figures/figure6_boxplot.png](implementation/figures/figure6_boxplot.png)

## 7) Final report

LaTeX source:
- [implementation/report/reproduction_report.tex](implementation/report/reproduction_report.tex)

Compiled PDF:
- [implementation/report/reproduction_report.pdf](implementation/report/reproduction_report.pdf)

## 8) Future work folder (implemented)

Artifacts generated in:
- [implementation/future_work](implementation/future_work)

Includes:
- [implementation/future_work/baseline_parity_results.csv](implementation/future_work/baseline_parity_results.csv)
- [implementation/future_work/hardware_metrics.json](implementation/future_work/hardware_metrics.json)
- [implementation/future_work/future_work_status.md](implementation/future_work/future_work_status.md)

Notes:
- Qiskit parity baselines were executed and recorded.
- BQSKit availability was verified and logged.
- Hardware-level metrics probing is implemented; real backend metrics require authenticated provider access.

## 9) Remaining gap to full paper parity

The workspace now demonstrates a complete implementation pipeline for the provided demo and generated artifacts. Full parity with all paper claims still depends on:
- complete official source release details beyond the demo subset,
- exact hyperparameter/protocol settings used for all paper experiments,
- hardware backend access/configuration for execution-level metrics.
