#!/usr/bin/env python3

from collections import defaultdict
from itertools import combinations, product
from typing import DefaultDict


def normalize_gene_id(gene: str) -> str:
    return str(gene).strip()


def canonical_edge(left: str, right: str) -> tuple[str, str]:
    left = normalize_gene_id(left)
    right = normalize_gene_id(right)
    return tuple(sorted((left, right)))


class DisjointSet:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def add(self, item: str) -> None:
        item = normalize_gene_id(item)
        self.parent.setdefault(item, item)

    def find(self, item: str) -> str:
        item = normalize_gene_id(item)
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        left = normalize_gene_id(left)
        right = normalize_gene_id(right)
        self.add(left)
        self.add(right)
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left != root_right:
            self.parent[root_right] = root_left


def build_components(edges: set[tuple[str, str]]) -> dict[str, list[str]]:
    dsu = DisjointSet()
    for left, right in edges:
        edge = canonical_edge(left, right)
        dsu.union(*edge)

    components: DefaultDict[str, list[str]] = defaultdict(list)
    for gene in dsu.parent:
        components[dsu.find(gene)].append(gene)
    return dict(components)


def validate_edges(edges: dict[tuple[str, str], tuple[float, str]]) -> dict[str, int]:
    noncanonical = 0
    for edge in edges:
        if edge != canonical_edge(*edge):
            noncanonical += 1
    return {"noncanonical_edges": noncanonical}


def build_ortholog_rows(
    edges: dict[tuple[str, str], tuple[float, str]],
    gene_to_ogu: dict[str, str],
    ogus: list[str],
) -> tuple[list[dict[str, str | int | float]], dict[str, int]]:
    validation_metrics = validate_edges(edges)
    edge_lookup = {canonical_edge(*edge) for edge in edges}
    components = build_components(edge_lookup)
    resolved_multi_ogu = 0
    excluded_too_sparse = 0
    unresolved_components = 0
    rows: list[dict[str, str | int | float]] = []

    for genes in components.values():
        by_ogu: DefaultDict[str, list[str]] = defaultdict(list)
        for gene in genes:
            gene = normalize_gene_id(gene)
            if gene not in gene_to_ogu:
                unresolved_components += 1
                continue
            by_ogu[gene_to_ogu[gene]].append(gene)

        present = [ogu for ogu in ogus if ogu in by_ogu]
        if len(present) < 2:
            excluded_too_sparse += 1
            continue

        selected, selection_status = select_best_subset(
            by_ogu=by_ogu,
            edge_evidence=edges,
            ogus=ogus,
        )
        if not selected:
            unresolved_components += 1
            continue
        if selection_status.startswith("resolved"):
            resolved_multi_ogu += 1

        ordered_genes = {ogu: selected.get(ogu, "") for ogu in ogus}
        present = [ogu for ogu in ogus if ordered_genes[ogu]]
        supported_edges: list[str] = []
        supported_edge_alignment_ids: list[str] = []
        alignment_score = 0.0
        for ogu_a, ogu_b in combinations(present, 2):
            gene_a = ordered_genes[ogu_a]
            gene_b = ordered_genes[ogu_b]
            edge = canonical_edge(gene_a, gene_b)
            if edge in edge_lookup:
                supported_edges.append(f"{ogu_a}-{ogu_b}")
                alignment_score += edges[edge][0]
                supported_edge_alignment_ids.append(f"{ogu_a}-{ogu_b}:{edges[edge][1]}")

        edge_count = len(supported_edges)
        min_connected_edges = len(present) - 1
        if edge_count < min_connected_edges:
            excluded_too_sparse += 1
            continue

        row: dict[str, str | int | float] = {ogu: ordered_genes[ogu] for ogu in ogus}
        row["edge_count"] = edge_count
        row["node_count"] = len(present)
        row["supported_edges"] = ",".join(supported_edges)
        row["edge_alignment_score_sum"] = round(alignment_score, 6)
        row["supported_edge_alignment_ids"] = ",".join(supported_edge_alignment_ids)
        row["selection_status"] = selection_status
        rows.append(row)

    rows.sort(
        key=lambda row: (
            -int(row["node_count"]),
            -int(row["edge_count"]),
            -float(row["edge_alignment_score_sum"]),
            *(str(row[ogu]) for ogu in ogus),
        )
    )

    metrics = {
        "components": len(components),
        "resolved_multi_gene_same_ogu_components": resolved_multi_ogu,
        "excluded_too_sparse_components": excluded_too_sparse,
        "unresolved_components": unresolved_components,
        **validation_metrics,
    }
    return rows, metrics


def select_best_subset(
    by_ogu: dict[str, list[str]],
    edge_evidence: dict[tuple[str, str], tuple[float, str]],
    ogus: list[str],
) -> tuple[dict[str, str], str]:
    present = [ogu for ogu in ogus if ogu in by_ogu]
    is_conflict = any(len(by_ogu[ogu]) != 1 for ogu in present)
    if not is_conflict:
        return {ogu: by_ogu[ogu][0] for ogu in present}, "direct"

    best_subset: dict[str, str] = {}
    best_key: tuple[int, int, float, tuple[str, ...]] | None = None
    best_support = 0

    for size in range(len(present), 1, -1):
        for chosen_ogus in combinations(present, size):
            gene_lists = [sorted(by_ogu[ogu]) for ogu in chosen_ogus]
            for genes in product(*gene_lists):
                gene_by_ogu = dict(zip(chosen_ogus, genes))
                supported_edges = []
                alignment_score = 0.0
                for ogu_a, ogu_b in combinations(chosen_ogus, 2):
                    edge = canonical_edge(gene_by_ogu[ogu_a], gene_by_ogu[ogu_b])
                    if edge in edge_evidence:
                        supported_edges.append(edge)
                        alignment_score += edge_evidence[edge][0]

                edge_count = len(supported_edges)
                if edge_count < size - 1:
                    continue

                tie_breaker = tuple(gene_by_ogu.get(ogu, "") for ogu in ogus)
                key = (size, edge_count, alignment_score, tie_breaker)
                if best_key is None or key > best_key:
                    best_key = key
                    best_subset = gene_by_ogu
                    best_support = edge_count

    if not best_subset:
        return {}, "unresolved"
    if best_support == len(best_subset) * (len(best_subset) - 1) // 2:
        return best_subset, "resolved"
    return best_subset, "resolved_sparse"
