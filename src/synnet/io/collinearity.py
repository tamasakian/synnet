#!/usr/bin/env python3

from pathlib import Path


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


def read_edges(path: Path, gene_to_subgenome: dict[str, str]) -> set[tuple[str, str]]:
    edges: set[tuple[str, str]] = set()
    with path.open() as handle:
        for line in handle:
            if not line.startswith(" "):
                continue
            cols = line.strip().split()
            if len(cols) < 4 or not cols[0].endswith(":"):
                continue
            gene_a, gene_b = cols[1], cols[2]
            if gene_a not in gene_to_subgenome or gene_b not in gene_to_subgenome:
                continue
            if gene_a == gene_b:
                continue
            edges.add(tuple(sorted((gene_a, gene_b))))
    return edges
