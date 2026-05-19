# synnet

`synnet` builds 1:1 syntenic ortholog networks from MCScanX `.collinearity`
files.

It treats every collinear gene pair as a network edge, joins edges across any
number of input collinearity files, and reports components where each subgenome
contributes at most one gene.

When a connected component contains multiple genes from the same subgenome,
`synnet` resolves the conflict by selecting the best 1:1 subset. It reads the
MCScanX alignment header, for example:

```text
## Alignment 25: score=417.0 e_value=2.6e-12 N=10 BI02&BII05_2 plus
```

and uses the alignment `score` as the weight for every gene-pair edge in that
alignment block. The alignment ID is retained for each supported subgenome pair.
Candidate subsets are prioritized by:

1. larger subgenome count
2. larger network score
3. larger summed alignment score

## Input

`synnet` needs:

- one or more MCScanX `.collinearity` files, or directories containing them
- one MCScanX-style bed4 coordinate file per subgenome

The coordinate files must have four tab-separated columns:

```text
chrom	gene	start	end
```

## Usage

Run from a source checkout:

```bash
PYTHONPATH=src python -m synnet.cli build \
  /path/to/mcscanx/results \
  --gff BI=/path/to/Cuscuta_campestris_BI.gff \
  --gff BII=/path/to/Cuscuta_campestris_BII.gff \
  --gff HI=/path/to/Cuscuta_chinensis_HI.gff \
  --gff HII=/path/to/Cuscuta_chinensis_HII.gff \
  --subgenomes BI,BII,HI,HII \
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
BI	BII	HI	HII	network_score	subgenome_count	supported_pairs	alignment_score	supported_alignment_ids	selection_status
```

Each row is one 1:1 syntenic ortholog network. Missing subgenomes are written as
empty cells. `network_score` is the number of supported subgenome-pair edges in
the network, so with four subgenomes the score ranges from 1 to 6.
`alignment_score` is the sum of the selected supporting MCScanX alignment
scores. `supported_alignment_ids` stores the MCScanX alignment ID used for each
supported subgenome pair, such as `BI-HI:25`. `selection_status` is `direct` for
already 1:1 components, `resolved` for conflict components resolved to a fully
connected subset, or `resolved_sparse` for conflict components resolved to a
connected but not fully connected subset.

`summary.txt` reports total row counts by subgenome completeness and by network
score.

## Example Summary

```text
total rows:      6720
1:1:1:1 rows:    4214
1:1:1:0 rows:    1455
1:1:0:0 rows:    1051

score 1: 1051
score 2: 151
score 3: 1409
score 4: 386
score 5: 647
score 6: 3076
```
