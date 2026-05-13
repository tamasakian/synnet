#!/usr/bin/env python3

from collections import defaultdict
from itertools import combinations


class DisjointSet:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def add(self, item: str) -> None:
        self.parent.setdefault(item, item)

    def find(self, item: str) -> str:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        self.add(left)
        self.add(right)
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left != root_right:
            self.parent[root_right] = root_left


def build_components(edges: set[tuple[str, str]]) -> dict[str, list[str]]:
    dsu = DisjointSet()
    for left, right in edges:
        dsu.union(left, right)

    components: dict[str, list[str]] = defaultdict(list)
    for gene in dsu.parent:
        components[dsu.find(gene)].append(gene)
    return components


def build_ortholog_rows(
    edges: set[tuple[str, str]],
    gene_to_subgenome: dict[str, str],
    subgenomes: list[str],
) -> tuple[list[dict[str, str | int]], dict[str, int]]:
    edge_lookup = set(edges)
    components = build_components(edges)
    excluded_multi_subgenome = 0
    excluded_too_sparse = 0
    rows: list[dict[str, str | int]] = []

    for genes in components.values():
        by_subgenome: dict[str, list[str]] = defaultdict(list)
        for gene in genes:
            by_subgenome[gene_to_subgenome[gene]].append(gene)

        present = [sub for sub in subgenomes if sub in by_subgenome]
        if len(present) < 2:
            excluded_too_sparse += 1
            continue
        if any(len(by_subgenome[sub]) != 1 for sub in present):
            excluded_multi_subgenome += 1
            continue

        ordered_genes = {
            sub: by_subgenome[sub][0] if sub in by_subgenome else ""
            for sub in subgenomes
        }
        supported_pairs: list[str] = []
        for sub_a, sub_b in combinations(present, 2):
            gene_a = ordered_genes[sub_a]
            gene_b = ordered_genes[sub_b]
            if tuple(sorted((gene_a, gene_b))) in edge_lookup:
                supported_pairs.append(f"{sub_a}-{sub_b}")

        network_score = len(supported_pairs)
        min_connected_edges = len(present) - 1
        if network_score < min_connected_edges:
            excluded_too_sparse += 1
            continue

        row: dict[str, str | int] = {sub: ordered_genes[sub] for sub in subgenomes}
        row["network_score"] = network_score
        row["subgenome_count"] = len(present)
        row["supported_pairs"] = ",".join(supported_pairs)
        rows.append(row)

    rows.sort(
        key=lambda row: (
            -int(row["subgenome_count"]),
            -int(row["network_score"]),
            *(str(row[sub]) for sub in subgenomes),
        )
    )

    metrics = {
        "components": len(components),
        "excluded_multi_gene_same_subgenome_components": excluded_multi_subgenome,
        "excluded_too_sparse_components": excluded_too_sparse,
    }
    return rows, metrics
