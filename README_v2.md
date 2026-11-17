# SprayNPray V2: assembly-backed prototype

`spraynpray_v2.py` is a separate command-line prototype. It runs **Prodigal →
DIAMOND blastp → assembly metadata lookup → taxonomic inference**, retaining
SprayNPray's contig summary style. The existing `spray-and-pray.py` remains the
original entry point.

## Requirements and quick start

Python **3.8+**, Prodigal **2.6.x**, and DIAMOND **2.x** on PATH. No additional
Python packages are needed. For example, create a dedicated environment:

```bash
conda create -n spraynpray-v2 -c conda-forge -c bioconda python=3.11 prodigal diamond
conda activate spraynpray-v2

python spraynpray_v2.py \
  --contigs assembly.fasta \
  --diamond-db /path/to/pseudomonas.dmnd \
  --protein-to-genome ~/databases/pseudomonas/prots/protein_to_genome.tsv \
  --assembly-info ~/databases/ncbi_assembly_info.tsv \
  --taxdump ~/databases/taxdump \
  --meta \
  --threads 32 \
  --outdir spraynpray_v2_results
```

Omit `--meta` for Prodigal's single-genome training mode. Prodigal is
single-threaded; `--threads` controls DIAMOND. Small assemblies may need meta
mode because single-genome mode requires sufficient sequence for training.

By default, Prodigal predicts all genes, then only the **first 25 ORFs per
contig**, ordered by their leftmost genomic coordinate regardless of strand,
are searched. Contigs with fewer genes use all available genes. Change this
with `--orfs-per-contig 50`, or search every predicted gene with
`--orfs-per-contig all`. Prodigal still sees the full input for gene prediction.
The summary lists **one actual top protein hit per selected ORF**, with its
organism name, protein accession, representative source assembly and identity.
An ORF without a qualifying alignment is listed as `no hit`.

Rerun the same command with the **same output directory** to reuse completed
Prodigal predictions and DIAMOND alignments automatically and regenerate the
reports. Input contig IDs must be unique in their first whitespace-delimited
header field. FASTA, map and assembly TSV inputs may be gzip compressed (`.gz`).

## Reference inputs

1. **DIAMOND database:** built from the same protein collection described by
   the protein-to-genome map. Use protein accessions as the FASTA IDs. Bare IDs
   such as `WP_002553999.1` and conventional `ref|WP_002553999.1|` hit IDs are
   supported. Protein and assembly accession versions are retained exactly;
   there is no silent version substitution.
2. **Protein-to-genome TSV:** two TAB-separated columns, protein accession
   then genome assembly accession. A header is optional. **Keep all rows when
   a protein belongs to multiple genomes.** Duplicate identical pairs are
   deduplicated in memory, but distinct source genomes are preserved. The map
   is streamed once, retaining only accessions found in the DIAMOND hits.
3. **Assembly metadata:** the original `ncbi_assembly_info.tsv` / NCBI assembly
   summary file, including the `#assembly_accession` TAB-separated header.
   Required columns are `assembly_accession`, `taxid`, `organism_name`.
   All provided columns are copied into the detailed match report. Exact
   accession matches are preferred. A GCF/GCA lookup through `gbrs_paired_asm`
   is allowed only if `paired_asm_comp` explicitly says `identical`; the report
   distinguishes the requested source accession from the actual metadata row.
   There is no automatic GCA/GCF prefix or accession version replacement.
4. **Taxdump directory:** `nodes.dmp` and `names.dmp` are required; `merged.dmp`
   is used when present to resolve replaced taxids. Deleted/missing taxids are
   left unresolved (they are absent from the current nodes). Other dump files
   are not required. Scientific names are loaded for the relevant lineages.

Use the original assembly metadata file rather than terminal-copied text.
The parser uses named TAB-separated columns and rejects mismatched row widths.

If a DIAMOND database has not been built yet, this creates one in your protein
directory (the temporary concatenated FASTA requires extra disk space):

```bash
cd ~/databases/pseudomonas/prots
gzip -cd -- ./*.faa.gz > pseudomonas_reference.faa
diamond makedb --in pseudomonas_reference.faa --db pseudomonas --threads 32
```

The protein-to-genome TSV is supplied to this script, **not** DIAMOND's
`--taxonmap`: these are different formats. Taxonomy is resolved after alignment.

## How inference works

1. Search selected ORFs using `diamond blastp --sensitive`. Defaults mirror
   the original program's alignment filters: E-value ≤ `1e-15`, query coverage
   ≥ 50%, subject coverage ≥ 50%. Minimum identity is configurable and defaults
   to 0%. Coverage uses the aligned coordinate span of one HSP.
2. Retain hits within **5% of the highest qualifying bitscore per ORF**
   (`--top-percent`). DIAMOND's `--top` reports this score window instead of an
   arbitrary fixed hit count. The parser reapplies filters/window selection,
   so reused hit files follow the same inference rules. One best HSP is used
   per protein hit, and repeated query/subject rows are collapsed.
3. Expand each protein hit to **every source genome in the provided map**,
   obtaining organism name and taxid from assembly metadata. An ORF receives
   the **lowest common ancestor (LCA)** of all these candidate taxids. A
   shared `WP_` protein therefore does not become one arbitrarily chosen strain.
   If any retained candidate lacks its mapping, assembly metadata, or a valid
   taxid, that ORF is marked `incomplete_metadata` and casts no taxonomic vote.
   An LCA at the root also casts no vote. Such ORFs still count toward the
   selected-ORF denominator used in the evidence coverage check below.
4. Each usable ORF contributes **one vote**, irrespective of how many hits,
   source genomes, or duplicate references it has. Its vote supports its LCA
   and the LCA's ancestors, never more specific descendants. Choose the
   deepest taxon supported by at least **80% of usable ORF votes**
   (`--consensus 0.8`). This is a heuristic support fraction, not a probability
   or calibrated confidence score. There must also be at least **3 usable
   ORFs** (`--min-tax-orfs 3`), covering at least **50% of selected ORFs**
   (`--min-classified-fraction 0.5`). Otherwise the contig remains unclassified
   with an explicit reason. A short contig with one or two selected ORFs is
   therefore unclassified by default; use `--min-tax-orfs 1` to relax this.
5. Predictions are capped at species by default; `--max-resolution taxid`
   permits deeper assignments supported by the same rules. Old `superkingdom`
   and newer `domain` nodes both populate the `domain` lineage column.

Example: ORF LCAs of species A, species A, species A, species A and genus G
support species A at 4/5 = 80%. Three A votes and two G votes support genus G
at 5/5, but species A at only 60%, so the result is genus G.

These calls describe similarity within the provided reference set. A
**Pseudomonas-only reference cannot establish that an unknown contig belongs
to Pseudomonas** or rule out contamination from taxa absent from that database.
Species/strain labels are provisional; this method does not implement ANI,
phylogenomic placement, or validated identity thresholds for rank assignment.
First-N sampling also examines only the beginning of long contigs and can
miss downstream discordance. Use `all` when broader contig evidence matters.

## Outputs

| File | Contents |
| --- | --- |
| `contigs.csv` | One row per input contig, including those with no ORFs or hits. Starts with SprayNPray's original seven column names, followed by ORF counts, prediction, support, status and rank columns. |
| `top_hits.csv` | Exactly one row per selected ORF, in contig/coordinate order, including `no_hit` rows. One top protein and one representative source assembly with metadata per matching ORF. |
| `matches.csv` | One row per retained ORF/protein/source-genome combination. Starts with `orf,taxa,top_hit`, then alignment metrics, mapping status and **every assembly metadata column**, prefixed `assembly_`. Unknown mappings are explicitly listed. |
| `orf_taxonomy.csv` | One row per selected ORF, with coordinates, candidate count, ORF LCA and status. |
| `selected_orfs.faa` | Only the ORFs submitted to DIAMOND. |
| `prodigal.faa`, `prodigal.gff` | Full Prodigal predictions (unless `--proteins` reused an existing file). |
| `diamond.tsv` | Unmodified DIAMOND output in the field order below (unless an existing file was reused). |
| `input_contigs.fasta` | Normalized input used by Prodigal; header IDs match the summary. |
| `run.log`, tool logs | Progress, commands and tool diagnostics. |
| `run_manifest.json` | Version, arguments, exact commands, counts and mapping diagnostics; written after successful completion. |
| `stage_cache.json` | Completed-stage signatures and output checksums, saved immediately after each successful stage to support automatic resume. |

`closest_blast_hits` lists exactly one **actual top hit** for each selected ORF,
not the ORF LCA or a list of every candidate genome. For example:

```text
1_1: Pseudomonas aeruginosa [WP_002553999.1; GCF_000233495.1; 99.50% identity]; 1_2: no hit
```

The best qualifying protein hit is selected by highest bitscore, then lowest
E-value, highest identity, and protein accession in lexical order. All tie
breaks are deterministic and do not depend on DIAMOND output order. A shared
protein accession can belong to many source assemblies: the display uses the
lexically first source with assembly metadata (or first source if none has
metadata). **This is a representative source, not evidence for an exact source
strain.** `top_hits.csv` includes its candidate-genome count and explicitly
flags representative source selection. `matches.csv` retains the complete
candidate matches, and `orf_taxonomy.csv` retains the conservative ORF LCA
assignments used for contig inference.

For `all`-ORF searches, only the summary display is capped at 30,000 UTF-16
units, with an omission notice and `summary_orfs_listed`/`summary_orfs_omitted`
counts. `top_hits.csv` still contains one row for every selected ORF, and all
evidence still contributes to inference. Proper CSV quoting preserves commas.

`GC-content` is a fraction (0–1) among A/C/G/T bases; ambiguous bases are
excluded. `hits_per_kb` is the number of **searched ORFs with qualifying
hits** per contig kb, not the density of every predicted gene. `cov` is
`Unknown` because this prototype does not accept coverage inputs.
`Average_AAI` retains the legacy label but reports the mean identity of each
ORF's selected top protein match, then the mean across matched ORFs. It is
**not a genome-wide AAI estimate**.

`support_fraction = supporting_orfs / classified_orfs`;
`classified_fraction = classified_orfs / selected_orfs`. Reports use empty
taxid/rank fields and `predicted_taxon=unclassified` when no call is made.

## Reusing predictions and alignments

### Automatic resume in the same output directory

After updating the code, **rerun your original command with the same
`--outdir`**. The log prints `Reusing completed Prodigal predictions` and
`Reusing completed DIAMOND hits` when the completed stages are compatible.
The reports, including the new one-hit-per-ORF summary, are regenerated.

Completed V2 `2.0.0a1`/`2.0.0a2` runs are recognized using their successful
`run_manifest.json`, normalized input, selected ORFs, alignment settings,
database path, and file timestamps. They are adopted into the new cache
without an alignment rerun when compatible. If those provenance files are
missing, supply `--proteins` and `--diamond-results` explicitly as below.

New runs save `stage_cache.json` immediately after each successful stage.
Prodigal cache entries include the normalized-contig SHA-256, mode and genetic
code. DIAMOND entries include the selected-protein SHA-256, database resolved
path/size/modification time, alignment settings and output format. Completed
output files are checked against their SHA-256 digests. A changed input,
database identity, search setting or damaged output reruns the affected stage.
Changing taxonomy/metadata inputs, inference thresholds or thread count reuses
compatible alignments. Database identity uses file metadata to avoid hashing
an entire large reference database on every run; use `--force-reblast` after
replacing a database while deliberately preserving its path, size and timestamp.

External tools write `.partial` files. These are promoted to completed output
files and marked reusable only after successful exit and parsing. A partial
or unmarked alignment is never automatically reused; a completed zero-hit
alignment is reusable. Executables are required only for stages that actually
need to run. `run_manifest.json` records whether each stage ran, was reused,
was adopted from a previous version, or was explicitly provided.

Use `--force-reblast` to rerun DIAMOND deliberately. A nonempty unrelated
directory is rejected; an existing recognizable V2 output directory is
accepted for resuming and refreshing reports.

### Explicit reuse from another directory

To revise inference without running external tools again, supply the full
Prodigal protein output and the previous DIAMOND TSV with the same inputs:

```bash
python spraynpray_v2.py \
  -i assembly.fasta \
  --proteins spraynpray_v2_results/prodigal.faa \
  --diamond-results spraynpray_v2_results/diamond.tsv \
  --protein-to-genome ~/databases/pseudomonas/prots/protein_to_genome.tsv \
  --assembly-info ~/databases/ncbi_assembly_info.tsv \
  --taxdump ~/databases/taxdump \
  --consensus 1.0 \
  -o spraynpray_v2_strict
```

Required DIAMOND fields (TAB-separated, no header):

```text
qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen
```

The existing SprayNPray `.blast` format is different and cannot be passed
directly. Reused `.faa` must have Prodigal's original `contig_ordinal` IDs and
` # start # end # strand # ` coordinate headers. Changing the sampling, filters,
or score window may require a new alignment search: a reused TSV cannot
recover hits or unsearched ORFs that were excluded in the original search.

## Verification

```bash
python -m unittest discover -s tests -p 'test_v2*.py' -v
```

The suite covers shared accessions, incomplete metadata, merged taxids, rank
resolution, sampling, actual single-top-hit selection, missing-hit rows,
large candidate sets, resumable stages, legacy-run adoption, cache invalidation,
failed/partial alignments and completed empty alignments. It also includes an optional real-tool
test using the repository's toy contigs and a synthetic reference/metadata
mapping. All 19 tests passed with Prodigal 2.6.3 and DIAMOND 2.2.9; the real-tool
test runs both 25-ORF and all-ORF searches and repeats each command to verify
automatic reuse without changing either tool output. This checks pipeline execution and
evidence accounting, not biological accuracy on independent reference data.

Reference documentation: [DIAMOND command-line options](https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options),
[Prodigal output](https://github.com/hyattpd/Prodigal/wiki/Understanding-the-Prodigal-Output),
[NCBI assembly summary format](https://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt).
