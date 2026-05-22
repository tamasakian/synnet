#!/usr/bin/env python3

import argparse
from pathlib import Path

from synnet.core.network import build_ortholog_rows
from synnet.io.bed4 import read_gene_ids
from synnet.io.collinearity import collect_collinearity_files, read_edges
from synnet.io.reports import write_network_tsv, write_summary_txt


def parse_gff_argument(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("--gff must be formatted as OGU=PATH")
    ogu, path = value.split("=", 1)
    if not ogu:
        raise argparse.ArgumentTypeError("OGU name is empty")
    return ogu, Path(path)


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
        help="OGU-to-MCScanX-bed4 mapping, formatted as OGU=PATH",
    )
    build.add_argument(
        "--ogus",
        required=True,
        help="comma-separated output order of Operational Genomic Units, for example BI,BII,HI,HII",
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
    ogus = [item.strip() for item in args.ogus.split(",") if item.strip()]
    gff_paths = dict(args.gff)
    missing = [ogu for ogu in ogus if ogu not in gff_paths]
    if missing:
        names = ", ".join(missing)
        raise SystemExit(f"missing --gff mapping for OGU(s): {names}")

    gene_to_ogu: dict[str, str] = {}
    for ogu in ogus:
        for gene in read_gene_ids(gff_paths[ogu]):
            gene_to_ogu[gene] = ogu

    collinearity_files = collect_collinearity_files(args.inputs)
    if not collinearity_files:
        raise SystemExit("no .collinearity files found")

    edges: dict[tuple[str, str], tuple[float, str]] = {}
    for path in collinearity_files:
        for edge, score in read_edges(path, gene_to_ogu).items():
            if score[0] >= edges.get(edge, (0.0, ""))[0]:
                edges[edge] = score

    rows, metrics = build_ortholog_rows(edges, gene_to_ogu, ogus)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    network_path = args.output_dir / "syntenic_orthologs_1to1_network.tsv"
    summary_path = args.output_dir / "summary.txt"
    write_network_tsv(network_path, rows, ogus)
    write_summary_txt(summary_path, rows, len(ogus), metrics)

    print(network_path)
    print(summary_path)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "build":
        run_build(args)


if __name__ == "__main__":
    main()
