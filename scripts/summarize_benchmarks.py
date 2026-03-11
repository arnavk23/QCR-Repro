from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


def load_rows(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8", newline="") as file:
        reader = csv.DictReader(file)
        for row in reader:
            rows.append(
                {
                    "depth": int(row["depth"]),
                    "iterations": int(row["iterations"]),
                    "seed": int(row["seed"]),
                    "start_len": int(row["start_len"]),
                    "end_len": int(row["end_len"]),
                    "replacements": int(row["replacements"]),
                    "runtime_sec": float(row["runtime_sec"]),
                    "equivalent": row["equivalent"] == "True",
                    "reduction_ratio": float(row["reduction_ratio"]),
                }
            )
    return rows


def group_summary(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[int, int], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[(int(row["depth"]), int(row["iterations"]))].append(row)

    summary: list[dict[str, object]] = []
    for (depth, iterations), items in sorted(groups.items()):
        end_lens = [int(i["end_len"]) for i in items]
        runtimes = [float(i["runtime_sec"]) for i in items]
        eq_rate = sum(bool(i["equivalent"]) for i in items) / len(items)
        summary.append(
            {
                "depth": depth,
                "iterations": iterations,
                "runs": len(items),
                "best_end_len": min(end_lens),
                "mean_end_len": round(statistics.mean(end_lens), 2),
                "std_end_len": round(statistics.pstdev(end_lens), 2),
                "mean_runtime_sec": round(statistics.mean(runtimes), 3),
                "std_runtime_sec": round(statistics.pstdev(runtimes), 3),
                "eq_pass_rate": round(eq_rate, 3),
            }
        )
    return summary


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create paper-style benchmark summaries.")
    parser.add_argument("--strict-csv", required=True)
    parser.add_argument("--loose-csv", required=True)
    parser.add_argument("--out-dir", default="results")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    strict_rows = load_rows(Path(args.strict_csv))
    loose_rows = load_rows(Path(args.loose_csv))

    strict_group = group_summary(strict_rows)
    loose_group = group_summary(loose_rows)

    write_csv(out_dir / "paper_style_strict.csv", strict_group)
    write_csv(out_dir / "paper_style_loose.csv", loose_group)

    strict_best_valid = min((r for r in strict_rows if bool(r["equivalent"])), key=lambda r: (int(r["end_len"]), float(r["runtime_sec"])))
    loose_best_valid = min((r for r in loose_rows if bool(r["equivalent"])), key=lambda r: (int(r["end_len"]), float(r["runtime_sec"])))

    with (out_dir / "paper_style_comparison.txt").open("w", encoding="utf-8") as file:
        file.write("Paper-Style Comparison\n")
        file.write("======================\n")
        file.write(f"Strict CSV: {args.strict_csv}\n")
        file.write(f"Loose CSV : {args.loose_csv}\n\n")

        file.write("Best valid run (strict tolerance)\n")
        file.write(str(strict_best_valid) + "\n\n")

        file.write("Best valid run (1e-3 tolerance)\n")
        file.write(str(loose_best_valid) + "\n\n")

        file.write("Grouped tables saved as:\n")
        file.write("- paper_style_strict.csv\n")
        file.write("- paper_style_loose.csv\n")

    print(f"Wrote {out_dir / 'paper_style_strict.csv'}")
    print(f"Wrote {out_dir / 'paper_style_loose.csv'}")
    print(f"Wrote {out_dir / 'paper_style_comparison.txt'}")


if __name__ == "__main__":
    main()
