#!/usr/bin/env python3

from pathlib import Path


def read_gene_ids(path: Path) -> list[str]:
    genes: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            genes.append(cols[1])
    return genes
