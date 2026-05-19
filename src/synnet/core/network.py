#!/usr/bin/env python3

from collections import defaultdict
from itertools import combinations, product


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
    edges: dict[tuple[str, str], tuple[float, str]],
    gene_to_subgenome: dict[str, str],
    subgenomes: list[str],
) -> tuple[list[dict[str, str | int]], dict[str, int]]:
    edge_lookup = set(edges)
    components = build_components(edge_lookup)
    resolved_multi_subgenome = 0
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

        selected, selection_status = select_best_subset(
            by_subgenome=by_subgenome,
            edge_evidence=edges,
            subgenomes=subgenomes,
        )
        if not selected:
            continue
        if selection_status.startswith("resolved"):
            resolved_multi_subgenome += 1

        ordered_genes = {sub: selected.get(sub, "") for sub in subgenomes}
        present = [sub for sub in subgenomes if ordered_genes[sub]]
        supported_pairs: list[str] = []
        supported_alignment_ids: list[str] = []
        alignment_score = 0.0
        for sub_a, sub_b in combinations(present, 2):
            gene_a = ordered_genes[sub_a]
            gene_b = ordered_genes[sub_b]
            edge = tuple(sorted((gene_a, gene_b)))
            if edge in edge_lookup:
                supported_pairs.append(f"{sub_a}-{sub_b}")
                alignment_score += edges[edge][0]
                supported_alignment_ids.append(f"{sub_a}-{sub_b}:{edges[edge][1]}")

        network_score = len(supported_pairs)
        min_connected_edges = len(present) - 1
        if network_score < min_connected_edges:
            excluded_too_sparse += 1
            continue

        row: dict[str, str | int] = {sub: ordered_genes[sub] for sub in subgenomes}
        row["network_score"] = network_score
        row["subgenome_count"] = len(present)
        row["supported_pairs"] = ",".join(supported_pairs)
        row["alignment_score"] = round(alignment_score, 6)
        row["supported_alignment_ids"] = ",".join(supported_alignment_ids)
        row["selection_status"] = selection_status
        rows.append(row)

    rows.sort(
        key=lambda row: (
            -int(row["subgenome_count"]),
            -int(row["network_score"]),
            -float(row["alignment_score"]),
            *(str(row[sub]) for sub in subgenomes),
        )
    )

    metrics = {
        "components": len(components),
        "resolved_multi_gene_same_subgenome_components": resolved_multi_subgenome,
        "excluded_too_sparse_components": excluded_too_sparse,
    }
    return rows, metrics


def select_best_subset(
    by_subgenome: dict[str, list[str]],
    edge_evidence: dict[tuple[str, str], tuple[float, str]],
    subgenomes: list[str],
) -> tuple[dict[str, str], str]:
    present = [sub for sub in subgenomes if sub in by_subgenome]
    is_conflict = any(len(by_subgenome[sub]) != 1 for sub in present)
    if not is_conflict:
        return {sub: by_subgenome[sub][0] for sub in present}, "direct"

    best_subset: dict[str, str] = {}
    best_key: tuple[int, int, float, tuple[str, ...]] | None = None
    best_support = 0

    for size in range(len(present), 1, -1):
        for chosen_subgenomes in combinations(present, size):
            gene_lists = [sorted(by_subgenome[sub]) for sub in chosen_subgenomes]
            for genes in product(*gene_lists):
                gene_by_subgenome = dict(zip(chosen_subgenomes, genes))
                supported_edges = []
                alignment_score = 0.0
                for sub_a, sub_b in combinations(chosen_subgenomes, 2):
                    edge = tuple(sorted((gene_by_subgenome[sub_a], gene_by_subgenome[sub_b])))
                    if edge in edge_evidence:
                        supported_edges.append(edge)
                        alignment_score += edge_evidence[edge][0]

                network_score = len(supported_edges)
                if network_score < size - 1:
                    continue

                tie_breaker = tuple(gene_by_subgenome.get(sub, "") for sub in subgenomes)
                key = (size, network_score, alignment_score, tie_breaker)
                if best_key is None or key > best_key:
                    best_key = key
                    best_subset = gene_by_subgenome
                    best_support = network_score

    if not best_subset:
        return {}, "unresolved"
    if best_support == len(best_subset) * (len(best_subset) - 1) // 2:
        return best_subset, "resolved"
    return best_subset, "resolved_sparse"
