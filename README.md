# synnet

`synnet` builds one to one syntenic ortholog networks
from MCScanX `.collinearity` files.

It treats every collinear gene pair as a network edge,
joins edges across any number of input collinearity files,
and reports components where each Operational Genomic Unit (OGU) contributes at
most one gene.

When a connected component contains multiple genes from the same OGU,
`synnet` resolves the conflict by selecting the best 1:1 subset.
It reads the MCScanX alignment header, for example:

```text
## Alignment 1: score=1300.0 e_value=2.6e-13 N=130 A1&B1 plus
```

and uses the alignment `score` as the weight for every gene-pair edge in that
alignment block. The alignment ID is retained for each supported OGU pair.
Candidate subsets are prioritized by:

1. larger node count
2. larger edge count
3. larger summed alignment score

## Input

`synnet` needs:

- one or more MCScanX `.collinearity` files, or directories containing them
- one MCScanX-style bed4 coordinate file per OGU

The coordinate files must have four tab-separated columns:

```text
chrom	gene	start	end
```

## Usage

Run from a source checkout:

```bash
PYTHONPATH=src python -m synnet.cli build \
  /path/to/mcscanx/results \
  --gff A1=/path/to/sp1_A1.gff \
  --gff A2=/path/to/sp1_A2.gff \
  --gff B1=/path/to/sp2_B1.gff \
  --gff B2=/path/to/sp2_B2.gff \
  --ogus A1,A2,B1,B2 \
  --output-dir /path/to/output
```

After installation, the same command is available as:

```bash
synnet build ...
```

## Output

The command writes two main files:

- `syntenic_orthologs_1to1_network.tsv`
- `summary.txt`

The network TSV has these columns:

```text
A1	A2	B1	B2	edge_count	node_count	supported_edges	edge_alignment_score_sum	supported_edge_alignment_ids	selection_status
```

Each row is one 1:1 syntenic ortholog network. Missing OGUs are written as empty
cells. `edge_count` is the number of supported OGU-pair edges in the network, so
with four OGUs the value ranges from 1 to 6. `edge_alignment_score_sum` is the
sum of the selected supporting MCScanX alignment scores.
`supported_edge_alignment_ids` stores the MCScanX alignment ID used for each
supported OGU pair, such as `A1-B1:1`. `selection_status` is `direct` for
already 1:1 components, `resolved` for conflict components resolved to a fully
connected subset, or `resolved_sparse` for conflict components resolved to a
connected but not fully connected subset.

`summary.txt` reports total row counts by node completeness and edge count.
