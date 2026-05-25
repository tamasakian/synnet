#!/usr/bin/env python3

from collections import Counter
from pathlib import Path


def write_network_tsv(path: Path, rows: list[dict[str, str | int | float]], ogus: list[str]) -> None:
    header = ogus + [
        "edge_count",
        "node_count",
        "supported_edges",
        "edge_alignment_score_sum",
        "supported_edge_alignment_ids",
        "selection_status",
    ]
    with path.open("w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[col]) for col in header) + "\n")


def write_summary_txt(
    path: Path,
    rows: list[dict[str, str | int | float]],
    ogu_count: int,
    metrics: dict[str, int],
) -> None:
    completeness = Counter(int(row["node_count"]) for row in rows)
    scores = Counter(int(row["edge_count"]) for row in rows)

    with path.open("w") as handle:
        handle.write(f"total rows:      {len(rows)}\n")
        for count in range(ogu_count, 1, -1):
            label = ":".join(["1"] * count + ["0"] * (ogu_count - count))
            handle.write(f"{label} rows:    {completeness[count]}\n")
        handle.write("\n")
        max_edges = ogu_count * (ogu_count - 1) // 2
        for edge_count in range(1, max_edges + 1):
            handle.write(f"edge_count {edge_count}: {scores[edge_count]}\n")
        handle.write("\n")
        handle.write(f"resolved conflict components: {metrics.get('resolved_multi_gene_same_ogu_components', 0)}\n")
        handle.write(f"excluded sparse components: {metrics.get('excluded_too_sparse_components', 0)}\n")
        handle.write(f"unresolved components: {metrics.get('unresolved_components', 0)}\n")
        handle.write(f"noncanonical edges: {metrics.get('noncanonical_edges', 0)}\n")
