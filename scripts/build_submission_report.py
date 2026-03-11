from __future__ import annotations

import csv
from pathlib import Path


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def key(row: dict[str, str]) -> tuple[str, str]:
    return row["depth"], row["iterations"]


def main() -> None:
    root = Path("results")
    strict_group = load_csv(root / "paper_style_strict.csv")
    loose_group = load_csv(root / "paper_style_loose.csv")
    strict_runs = load_csv(root / "benchmark_reducer.csv")
    loose_runs = load_csv(Path("results_tol1e3") / "benchmark_reducer.csv")

    loose_map = {key(r): r for r in loose_group}

    combined_rows: list[dict[str, str]] = []
    for s in strict_group:
        k = key(s)
        l = loose_map[k]
        combined_rows.append(
            {
                "depth": s["depth"],
                "iterations": s["iterations"],
                "runs": s["runs"],
                "best_end_len_strict": s["best_end_len"],
                "mean_end_len_strict": s["mean_end_len"],
                "std_end_len_strict": s["std_end_len"],
                "mean_runtime_sec_strict": s["mean_runtime_sec"],
                "eq_pass_rate_strict": s["eq_pass_rate"],
                "best_end_len_loose": l["best_end_len"],
                "mean_end_len_loose": l["mean_end_len"],
                "std_end_len_loose": l["std_end_len"],
                "mean_runtime_sec_loose": l["mean_runtime_sec"],
                "eq_pass_rate_loose": l["eq_pass_rate"],
            }
        )

    out_csv = root / "submission_report_table.csv"
    with out_csv.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(combined_rows[0].keys()))
        writer.writeheader()
        writer.writerows(combined_rows)

    strict_valid = [r for r in strict_runs if r["equivalent"] == "True"]
    loose_valid = [r for r in loose_runs if r["equivalent"] == "True"]

    best_strict = min(strict_valid, key=lambda r: (int(r["end_len"]), float(r["runtime_sec"])))
    best_loose = min(loose_valid, key=lambda r: (int(r["end_len"]), float(r["runtime_sec"])))

    out_md = root / "submission_report.md"
    with out_md.open("w", encoding="utf-8") as f:
        f.write("# Reproduction Benchmark Report\n\n")
        f.write("## Best configurations\n")
        f.write(f"- Strict (1e-5): depth={best_strict['depth']}, iterations={best_strict['iterations']}, seed={best_strict['seed']}, end_len={best_strict['end_len']}, runtime_sec={best_strict['runtime_sec']}\\n\n")
        f.write(f"- Loose (1e-3): depth={best_loose['depth']}, iterations={best_loose['iterations']}, seed={best_loose['seed']}, end_len={best_loose['end_len']}, runtime_sec={best_loose['runtime_sec']}\\n\n")
        f.write("## Combined table\n\n")
        f.write("| depth | iterations | runs | best_strict | mean_strict | eq_strict | best_loose | mean_loose | eq_loose |\n")
        f.write("|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n")
        for r in combined_rows:
            f.write(
                f"| {r['depth']} | {r['iterations']} | {r['runs']} | {r['best_end_len_strict']} | {r['mean_end_len_strict']} | {r['eq_pass_rate_strict']} | {r['best_end_len_loose']} | {r['mean_end_len_loose']} | {r['eq_pass_rate_loose']} |\n"
            )

    print(f"Wrote {out_csv}")
    print(f"Wrote {out_md}")


if __name__ == "__main__":
    main()
