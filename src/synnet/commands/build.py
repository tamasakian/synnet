#!/usr/bin/env python3

import argparse
from pathlib import Path

from synnet.core.network import build_ortholog_rows
from synnet.io.bed4 import read_gene_ids
from synnet.io.collinearity import collect_collinearity_files, read_edges
from synnet.io.reports import write_network_tsv, write_summary_txt


def parse_gff_argument(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("--gff must be formatted as SUBGENOME=PATH")
    subgenome, path = value.split("=", 1)
    if not subgenome:
        raise argparse.ArgumentTypeError("subgenome name is empty")
    return subgenome, Path(path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="synnet",
        description="Build 1:1 syntenic ortholog networks from MCScanX collinearity files.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("build", help="build a syntenic ortholog network")
    build.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="MCScanX .collinearity files or directories containing them",
    )
    build.add_argument(
        "--gff",
        action="append",
        required=True,
        type=parse_gff_argument,
        help="subgenome-to-MCScanX-bed4 mapping, formatted as SUBGENOME=PATH",
    )
    build.add_argument(
        "--subgenomes",
        required=True,
        help="comma-separated output order, for example BI,BII,HI,HII",
    )
    build.add_argument(
        "-o",
        "--output-dir",
        required=True,
        type=Path,
        help="output directory",
    )
    return parser


def run_build(args: argparse.Namespace) -> None:
    subgenomes = [item.strip() for item in args.subgenomes.split(",") if item.strip()]
    gff_paths = dict(args.gff)
    missing = [sub for sub in subgenomes if sub not in gff_paths]
    if missing:
        names = ", ".join(missing)
        raise SystemExit(f"missing --gff mapping for subgenome(s): {names}")

    gene_to_subgenome: dict[str, str] = {}
    for subgenome in subgenomes:
        for gene in read_gene_ids(gff_paths[subgenome]):
            gene_to_subgenome[gene] = subgenome

    collinearity_files = collect_collinearity_files(args.inputs)
    if not collinearity_files:
        raise SystemExit("no .collinearity files found")

    edges: set[tuple[str, str]] = set()
    for path in collinearity_files:
        edges.update(read_edges(path, gene_to_subgenome))

    rows, _metrics = build_ortholog_rows(edges, gene_to_subgenome, subgenomes)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    network_path = args.output_dir / "syntenic_orthologs_1to1_network.tsv"
    summary_path = args.output_dir / "summary.txt"
    write_network_tsv(network_path, rows, subgenomes)
    write_summary_txt(summary_path, rows, len(subgenomes))

    print(network_path)
    print(summary_path)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "build":
        run_build(args)


if __name__ == "__main__":
    main()
