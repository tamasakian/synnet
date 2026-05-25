#!/usr/bin/env python3

import re
from pathlib import Path

from synnet.core.network import canonical_edge, normalize_gene_id


def collect_collinearity_files(paths: list[Path]) -> list[Path]:
    files: list[Path] = []
    for path in paths:
        if path.is_dir():
            files.extend(sorted(path.rglob("*.collinearity")))
        elif path.is_file():
            files.append(path)
        else:
            raise FileNotFoundError(f"input path not found: {path}")
    return sorted(dict.fromkeys(files))


ALIGNMENT_RE = re.compile(r"^## Alignment\s+(\d+):\s+score=([^\s]+)")


def read_edges(path: Path, gene_to_ogu: dict[str, str]) -> dict[tuple[str, str], tuple[float, str]]:
    edges: dict[tuple[str, str], tuple[float, str]] = {}
    current_alignment_id = ""
    current_score = 0.0
    with path.open() as handle:
        for line in handle:
            if line.startswith("## Alignment"):
                match = ALIGNMENT_RE.match(line)
                if match:
                    current_alignment_id = match.group(1)
                    current_score = float(match.group(2))
                continue
            if not line.startswith(" "):
                continue
            cols = line.strip().split()
            if len(cols) < 4 or not cols[0].endswith(":"):
                continue
            gene_a = normalize_gene_id(cols[1])
            gene_b = normalize_gene_id(cols[2])
            if gene_a not in gene_to_ogu or gene_b not in gene_to_ogu:
                continue
            if gene_a == gene_b:
                continue
            edge = canonical_edge(gene_a, gene_b)
            previous_score = edges.get(edge, (0.0, ""))[0]
            if current_score >= previous_score:
                edges[edge] = (current_score, current_alignment_id)
    return edges
