#!/usr/bin/env python3

from collections import defaultdict
from itertools import combinations
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
    matched_edges = select_greedy_ogu_pair_matches(edges, gene_to_ogu, ogus)
    matched_edges, prune_metrics = prune_conflicting_components(matched_edges, gene_to_ogu, ogus)
    edge_lookup = set(matched_edges)
    components = build_components(edge_lookup)
    excluded_multi_ogu = 0
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

        if any(len(by_ogu[ogu]) != 1 for ogu in present):
            excluded_multi_ogu += 1
            continue

        ordered_genes = {ogu: by_ogu[ogu][0] if ogu in by_ogu else "" for ogu in ogus}
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
        row["selection_status"] = "matched"
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
        "input_edges": len(edges),
        "matched_edges": len(matched_edges),
        "components": len(components),
        "excluded_multi_gene_same_ogu_components": excluded_multi_ogu,
        "excluded_too_sparse_components": excluded_too_sparse,
        "unresolved_components": unresolved_components,
        **prune_metrics,
        **validation_metrics,
    }
    return rows, metrics


def component_has_ogu_conflict(
    genes: list[str],
    gene_to_ogu: dict[str, str],
    ogus: list[str],
) -> bool:
    counts = {ogu: 0 for ogu in ogus}
    for gene in genes:
        gene = normalize_gene_id(gene)
        if gene not in gene_to_ogu:
            continue
        ogu = gene_to_ogu[gene]
        if ogu not in counts:
            continue
        counts[ogu] += 1
        if counts[ogu] > 1:
            return True
    return False


def prune_conflicting_components(
    edges: dict[tuple[str, str], tuple[float, str]],
    gene_to_ogu: dict[str, str],
    ogus: list[str],
) -> tuple[dict[tuple[str, str], tuple[float, str]], dict[str, int]]:
    pruned = dict(edges)
    removed_edges = 0

    while True:
        components = build_components(set(pruned))
        edges_to_remove: set[tuple[str, str]] = set()

        for genes in components.values():
            if not component_has_ogu_conflict(genes, gene_to_ogu, ogus):
                continue

            gene_set = set(genes)
            component_edges = [
                edge
                for edge in pruned
                if edge[0] in gene_set and edge[1] in gene_set
            ]
            if not component_edges:
                continue

            weakest_edge = min(
                component_edges,
                key=lambda edge: (pruned[edge][0], edge[0], edge[1]),
            )
            edges_to_remove.add(weakest_edge)

        if not edges_to_remove:
            break

        for edge in edges_to_remove:
            pruned.pop(edge, None)
            removed_edges += 1

    return pruned, {"pruned_conflict_edges": removed_edges}


def select_greedy_ogu_pair_matches(
    edges: dict[tuple[str, str], tuple[float, str]],
    gene_to_ogu: dict[str, str],
    ogus: list[str],
) -> dict[tuple[str, str], tuple[float, str]]:
    ogu_rank = {ogu: index for index, ogu in enumerate(ogus)}
    used_by_ogu_pair: dict[tuple[str, str], set[str]] = defaultdict(set)
    matched: dict[tuple[str, str], tuple[float, str]] = {}

    sorted_edges = sorted(
        edges.items(),
        key=lambda item: (item[1][0], item[0][0], item[0][1]),
        reverse=True,
    )

    for edge, evidence in sorted_edges:
        left, right = canonical_edge(*edge)
        if left not in gene_to_ogu or right not in gene_to_ogu:
            continue

        ogu_left = gene_to_ogu[left]
        ogu_right = gene_to_ogu[right]
        if ogu_left == ogu_right:
            continue

        ogu_pair = tuple(sorted((ogu_left, ogu_right), key=ogu_rank.get))
        used = used_by_ogu_pair[ogu_pair]
        if left in used or right in used:
            continue

        matched[(left, right)] = evidence
        used.add(left)
        used.add(right)

    return matched
