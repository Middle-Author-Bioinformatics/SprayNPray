# Cohort summary and genome ANI map

`spraynpray_compare.py` is a standalone script that combines sample-named
SprayNPray V2 CSV files, finds the source assemblies represented by their
protein hits, computes all-vs-all genome ANI using FastANI, and produces a
map with **black sample nodes** and **gray reference nodes**.

## Install and run on TheBelly

```bash
conda create -n spraynpray-compare -c conda-forge -c bioconda \
  python=3.11 numpy scipy matplotlib fastani ncbi-datasets-cli
conda activate spraynpray-compare

python ~/bin/SprayNPray/spraynpray_compare.py \
  --results ~/PITT/Rapsinski/long/hybrid_results/spraynpray_combined \
  --genomes ~/PITT/Rapsinski/long/hybrid_results/assemblies_combined \
  --threads 32 \
  --outdir ~/PITT/Rapsinski/long/hybrid_results/cohort_ani
```

Adjust the script path to where you cloned/downloaded it. Inputs should pair
`sample_3_ONT002.csv` with `sample_3_ONT002.fa`. The default genome extension is
`.fa`; change it with `--genome-extension .fna` or `.fa.gz`.

The default includes **every distinct GCA/GCF reference accession** in the
hits. It automatically downloads missing reference **nucleotide genomes**
using NCBI Datasets, keeping exact accession versions. Downloads are retained
under the output directory for subsequent runs. No protein identity is
converted into ANI, and no reference ANI is fabricated from a hit frequency.

If reference FASTAs are already available locally, add:

```bash
--reference-dir ~/databases/pseudomonas/genomes
```

The folder may contain `.fa`, `.fna`, `.fasta` or their gzip-compressed
equivalents, with exact accessions in the filenames, for example
`GCF_000233495.1.fna` or `GCF_000233495.1_ASM23349v1_genomic.fna`. Nested
NCBI Datasets folders are also recognized. Use complete genomic FASTAs rather
than CDS or protein files. Ambiguous duplicate reference files are rejected.
`--no-download` uses only local/cached references. Missing exact reference
versions stop the full analysis and are listed in `missing_references.csv`;
`--allow-missing-references` explicitly permits them as isolated gray nodes
without ANI measurements. No accession version or GCA/GCF substitution occurs.

For a smaller initial reference set, use `--top-references-per-sample 3`.
This selects the three references with the most hits in each sample and uses
their union for ANI. All hits and references still appear in the summary
tables, with `selected_for_ani` marking the chosen references.

## Master summary: percentages and identities

The input folder must contain **one sample-named CSV per sample**. Both the
V2 `contigs.csv` format with one top hit per selected ORF and the V2
`top_hits.csv` format are accepted; copy/rename the complete source files to
sample filenames before combining them. Do not put the output directory's
summary tables among the inputs.

For each sample/reference assembly, `master_summary.csv` reports:

| Column | Meaning |
| --- | --- |
| `orf_hits` | Number of selected ORFs whose displayed top hit uses that source assembly. |
| `hit_percent` | `100 × orf_hits / matched_orfs`. Denominator: all ORFs with qualifying hits in that sample. |
| `selected_orf_percent` | `100 × orf_hits / selected_orfs`. Includes the searched no-hit ORFs in the denominator. |
| `mean_protein_identity` | Arithmetic mean of the individual ORF hit identities for that assembly. Every matched ORF has equal weight. |
| `min_protein_identity`, `max_protein_identity` | Range of those identities. |

No-hit ORFs are excluded from `hit_percent` but reported in
`sample_overview.csv`. Hits without an assembly mapping form an explicit
`unmapped` group, so hit percentages still sum to 100% for samples with hits.
Samples with no hits remain in the overview with zero matched ORFs and no
invented assembly-match row. The mean identity reflects the precision available
in the input: contig summaries contain identities rounded by the upstream
program; `top_hits.csv` retains the original numeric precision.

An old oversized summary listing many genomes per ORF is rejected. A summary
with omitted ORFs is also rejected rather than silently producing biased
percentages; use the complete `top_hits.csv` for that sample instead.

**Interpretation:** these are percentages of the selected **protein hits**,
not genomic abundance or the fraction of genome DNA belonging to each strain.
Shared `WP_` accessions can use a deterministically chosen representative
source assembly. Frequent hit labels therefore do not by themselves establish
the closest reference strain. Use the genome ANI comparisons to assess that.
First-25-per-contig sampling is also not a random whole-genome gene sample.

To generate just the summary tables with standard Python and no downloads:

```bash
python spraynpray_compare.py \
  --results /path/to/spraynpray_combined \
  --summary-only \
  --outdir /path/to/cohort_summary
```

## ANI calculation

The default `--ani-scope samples` compares **all sample/sample pairs and all
sample/reference pairs in both directions**, skipping reference/reference
comparisons. It uses a sample-query job against the full sequence set and a
reference-query job against the samples. This keeps the comparisons needed
for sample relationships and nearest references while avoiding quadratic work
among references. Use `--ani-scope all` for the previous full-union all-vs-all
analysis, including references against one another.

For S samples and R references, the requested directional comparisons drop
from `(S + R)^2` to `S^2 + 2SR`. For example, 40 samples and 1,000 references
require 81,600 instead of 1,081,600 directional comparisons. This is a work
count, not a guaranteed wall-time speedup. All sample/reference ANI evidence
is retained. For further reductions, `--top-references-per-sample 3` narrows
the reference set; nearest matches then apply only to that selected set.
All FASTAs are normalized into unambiguous staged paths, and duplicate
sequence content still retains separate biological/sample labels.

FastANI output reports query, reference, estimated nucleotide identity,
mapped fragments and total query fragments. The script retains both directions
and their aligned fractions (`mapped / total`). For a pair to inform the
default graph/matrix, **both directions must be reported** and **each must
have an aligned fraction of at least 0.5** (`--min-af 0.5`). The plotted ANI
is the mean of the two directional estimates. Raw directional values and
excluded comparisons remain in the output tables with a reason.

FastANI normally does not report pairs much below approximately 80% ANI.
An unreported pair is **missing**, not 0% ANI and not 100% ANI. Estimates
based on low aligned fractions may represent shared fragments rather than
the genome as a whole; these are excluded from the default graph. The
`--frag-len` default is 3000 bases. Very fragmented/short genomes may therefore
have few usable fragments and no reported comparisons. The default
`--min-ani 80` can be increased to focus the map on closer relationships.

The nucleotide sequences are required to compute reference ANI. CSV protein
identities alone cannot provide it. FastANI is an estimator; results are not
an exact BLAST-based ANI calculation or a variant-calling validation.

## Graph interpretation

The layout fits **Euclidean distances to `100 − mean ANI`**, using all
accepted pairwise measurements in a connected component. It starts from
classical multidimensional scaling, then minimizes the squared differences
between measured distances and plotted distances using only observed pairs.
Missing pairs introduce no ANI constraint. Shortest-path distances are used
only to initialize the coordinates, not as new ANI measurements.

This is a **2D approximation**. Arbitrary pairwise ANI distances cannot all be
preserved exactly on a flat page. Exact ANI and aligned fractions appear on
hover and in `ani_pairs.csv`; `layout_diagnostics.csv` gives normalized stress
and the maximum absolute distance error for each connected component. A tiny
distance floor of 0.001 ANI percentage points stabilizes nearly identical
genomes. Identical/near-identical genomes can overlap visually.

Disconnected components are placed side by side. **Distances between those
components have no ANI meaning.** Missing sequence references and genomes
with no accepted comparisons remain visible as isolated nodes. This map is
not a phylogenetic tree and carries no branch-length or ancestry inference.

For readability, the display uses each node's five nearest accepted neighbors
plus a minimum spanning forest to preserve connectivity. `--neighbors 0`
draws every accepted comparison. All accepted pairs still inform the layout
regardless of which lines are drawn.

`ani_network.html` is self-contained and works offline. Hover nodes and lines,
zoom with the mouse wheel, drag to pan, toggle labels, or find a genome in the
selector. Reference accession/organism labels are now spread into two outer
columns, with dashed leader lines linking labels to the unchanged genome
points. The labels have separate nonoverlapping boxes, even when genome dots
coincide. Dashed leaders are label connectors, not ANI edges. Tall label sets
can be scrolled; reference labels are visible by default and can be toggled.
PNG and SVG versions are also produced for export. Colors are fixed: black
samples, gray references.

## Concise nearest-reference result and helper

Every full run now writes `nearest_references.csv`: **one row per sample**,
with the closest included reference, estimated overall genome ANI, exact ANI
distance (`100 - ANI`), minimum reciprocal aligned fraction, tie count and
status. The closest reference is the one with the **highest accepted measured
ANI**, not the shortest approximate 2D line. Overall ANI here is the mean of
the reciprocal genome ANI estimates for that sample/reference pair; it is not
protein identity or an average over all references.

Equal ANI ties (within 0.000001 percentage points) are counted; the displayed
reference is chosen by highest minimum aligned fraction, then accession.
No accepted reference produces a blank ANI and `no_accepted_reference`.
The distance column uses `100 - ANI` without the tiny numerical floor used
in fitting the graph. Values are rounded to five decimals for readability;
the full observations remain in `ani_pairs.csv`.

To add this summary and improve labels on an **existing** completed output:

```bash
python spraynpray_helper.py --analysis-dir /path/to/cohort_ani
```

Keep the helper beside the updated `spraynpray_compare.py`. It reads the saved
`ani_pairs.csv`, `network_nodes.csv` and `layout_diagnostics.csv`, preserves
the saved genome coordinates and ANI measurements, and rewrites only the
nearest-reference result and plot outputs/plot metadata. It runs no gene
prediction, alignment, reference download or ANI calculation. Use
`--summary-only` for just the concise CSV without plotting dependencies.

See [README_pipeline.md](README_pipeline.md) for the genome-folder wrapper.

## Outputs and reuse

| File | Contents |
| --- | --- |
| `master_summary.csv` | Per-sample/per-reference hit percentages and mean protein identities. |
| `sample_overview.csv` | Selected, matched, unmapped and no-hit ORF counts per sample. |
| `orf_hits.csv` | Parsed individual ORF evidence for audit/reanalysis. |
| `references.csv` | All observed references, hit counts, sample counts and ANI inclusion flags. |
| `nearest_references.csv` | One closest accepted reference and overall estimated genome ANI per sample. |
| `ani_pairs.csv` | Every requested unordered pair, both directional estimates/AFs, mean ANI and acceptance status. Reference-only pairs are omitted in the default sample-focused scope. |
| `ani_directional.csv` | Every actual directional FastANI observation. |
| `ani_matrix.csv` | Symmetric mean-ANI matrix for accepted pairs; unmeasured/excluded cells are blank. |
| `aligned_fraction_matrix.csv` | Minimum reciprocal aligned fraction for accepted pairs. |
| `network_nodes.csv`, `network_edges.csv` | Node identities, coordinates and the displayed edges with ANI. |
| `layout_diagnostics.csv` | Embedding error by component. |
| `ani_network.html`, `.png`, `.svg` | Interactive/offline and static maps. |
| `work/` | Staged input genomes, downloaded references, raw FastANI output, logs and the ANI cache. |

Rerun the same command with the same output directory to reuse downloads and
completed ANI. The cache validates normalized genome SHA-256 digests,
fragment length, engine settings and the raw output digest. Changing plot
thresholds, edge counts or thread count can reuse ANI. Failed/partial tool
outputs are never marked complete. `--force-ani` deliberately recomputes ANI.
Existing version 0.1 full all-vs-all results are reusable by the new default:
they already contain every required sample-linked comparison. A sample-only
cache cannot satisfy a later full `--ani-scope all` request. The two jobs in
sample-focused mode are also cached independently, so an interruption during
the second job does not repeat the completed first job.
`--fastani-results` can explicitly reuse a five-column FastANI output whose
genome fields use the original absolute input paths, staged paths or node IDs.

## Validation and sources

```bash
python -m unittest discover -s tests -p 'test_compare.py' -v
```

Tests cover hit counts/percentages/means, no-hit/unmapped handling, truncated
inputs, exact accession versions, NCBI Datasets command construction and
archive extraction, directional ANI/coverage filtering, missing pairs,
cache reuse/invalidation and failed searches, an exact-distance embedding
fixture, plus a real FastANI run with synthetic related DNA and offline figure
generation. The real-tool test requires `fastANI` on PATH. Reference downloads
are tested with a simulated NCBI package; live NCBI availability is external.

Primary references:

- [FastANI usage and output definitions](https://github.com/ParBLiSS/FastANI)
- [FastANI methods paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC6269478/)
- [NCBI Datasets genome downloads](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/)
- [NCBI Datasets genome package structure](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/)
