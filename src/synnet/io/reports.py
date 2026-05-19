#!/usr/bin/env python3

from collections import Counter
from pathlib import Path


def write_network_tsv(path: Path, rows: list[dict[str, str | int | float]], subgenomes: list[str]) -> None:
    header = subgenomes + [
        "network_score",
        "subgenome_count",
        "supported_pairs",
        "alignment_score",
        "supported_alignment_ids",
        "selection_status",
    ]
    with path.open("w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[col]) for col in header) + "\n")


def write_summary_txt(
    path: Path,
    rows: list[dict[str, str | int | float]],
    subgenome_count: int,
    metrics: dict[str, int],
) -> None:
    completeness = Counter(int(row["subgenome_count"]) for row in rows)
    scores = Counter(int(row["network_score"]) for row in rows)

    with path.open("w") as handle:
        handle.write(f"total rows:      {len(rows)}\n")
        for count in range(subgenome_count, 1, -1):
            label = ":".join(["1"] * count + ["0"] * (subgenome_count - count))
            handle.write(f"{label} rows:    {completeness[count]}\n")
        handle.write("\n")
        max_score = subgenome_count * (subgenome_count - 1) // 2
        for score in range(1, max_score + 1):
            handle.write(f"score {score}: {scores[score]}\n")
        handle.write("\n")
        handle.write(f"resolved conflict components: {metrics.get('resolved_multi_gene_same_subgenome_components', 0)}\n")
        handle.write(f"excluded sparse components: {metrics.get('excluded_too_sparse_components', 0)}\n")
