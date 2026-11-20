# Archive cleanup and portable cohort dashboard

Two standalone Python 3.9+ scripts use only the standard library. They work on
an existing `hybrid_results` directory; they do not rerun assembly, ANI or
GToTree, and do not require the GToTree conda environment.

## Build the dashboard first

From your `hybrid_results` directory:

```bash
python ~/bin/SprayNPray/spraynpray_dashboard.py \
  --root . \
  --outdir cohort_dashboard \
  --tree cohort_phylogeny_fresh.BldVHy/GToTree_output/GToTree_output.tre \
  --title "Rapsinski isolate genome comparison"
```

Use your actual tree path if its directory name differs. Omitting `--tree`
auto-detects a tree in `cohort_phylogen*/GToTree_output/*.tre` or
`cohort_phylogen*/runs/attempt_*/*.tre`; multiple valid trees require an explicit
choice. A dashboard can still be built without a tree or ANI, with those
analyses marked unavailable.

The script reads:

- `assemblies_combined/*.{fa,fna,fasta}` (also `.gz`), plus `sample_*` directories.
  A final `unicycler/assembly.fasta` or `flye/assembly.fasta` is a fallback when
  the sample is absent from the combined set. This makes an assembled sample
  omitted from downstream analyses visible rather than calling it a failure.
- `sample_*/quast/report.tsv` and `unicycler.done` / `flye.done`. Completion
  markers may be empty. A valid FASTA with no marker is reported as “FASTA
  present,” and a marker alone does not imply successful assembly.
- `cohort_ani/ani_pairs.csv` and `network_nodes.csv` from the comparison helper.
  Only accepted reciprocal genome ANI is used for grouping and nearest references;
  missing directions, rejected coverage and absent comparisons remain unknown.
- The native Newick file, adjacent `Genomes_summary_info.tsv`, and the nearest
  ancestor `selected_genomes.csv` to associate GToTree labels with sample IDs.
  Repeated underscores, renamed labels and staged `gNNNNNN` identifiers are
  handled using the selection table. Keep this table with the native results.

It writes `cohort_dashboard/dashboard.html` and `cohort_dashboard.zip`. Unzip and
open `dashboard.html` in a modern browser; no server, external JavaScript, CDN,
Python installation or internet connection is needed by the recipient. Keep the
HTML, `data/` and `assemblies/` together. The ZIP contains only the generated
report files, summary CSVs, original Newick and compressed final **input**
assemblies, not raw reads, reference FASTAs or the full SprayNPray output.
The HTML embeds the report data, so its plots still work if shared alone, but
its file downloads need the accompanying directories. Only share the ZIP with
people authorized to receive the included genome sequences.

The report provides:

- Sample-by-sample assembly availability, assembler, length, contig count, N50,
  GC, circular-contig annotations and assembly downloads. Missing or invalid
  final FASTAs are explicit; “completed” indicates output/marker presence,
  not a biological quality threshold.
- Per-sample QUAST details and GToTree marker recovery/redundancy when available.
  Marker recovery is not CheckM completeness, and marker redundancy is not a
  contamination estimate. No absent QC result is fabricated.
- A sample ANI heatmap with exact values, aligned fractions, rejected/missing
  status, nearest accepted reference and an adjustable grouping threshold.
- A phylogeny with one label per row, samples shown by default, optional
  references, label search, click details, zoom, native Newick and SVG export.
  Pruned paths retain their summed branch lengths. When lengths are incomplete,
  the displayed tree uses and identifies topological spacing.
- A concise interpretation and a follow-up agenda for reconciling isolate
  origins, collection dates and culture history.

### Interpreting similarity groups

Default `--group-ani 99` is an exploratory display setting, **not a validated
strain or lineage cutoff**. Complete-linkage groups merge only if every
cross-pair has accepted ANI at or above the chosen threshold. This avoids
joining distant samples through a chain of similar neighbors. Missing pairs
prevent a merge; samples with no accepted sample comparisons are “Unresolved.”

The overview reports distinct genomic similarity groups only when measured
sample pairs fall below the current threshold. It gives the observed range and
missing-pair count. Separate groups caused solely by missing comparisons do not
support a multiple-lineage conclusion. Reference comparisons do not determine
sample groups. ANI does not establish ancestry, transmission or chronology;
interpret the marker tree, sequencing/QC results and isolate metadata together.
For context on genome ANI, see the [FastANI primary publication](https://doi.org/10.1038/s41467-018-07641-9).

An existing bundle can be refreshed in the same folder; unrelated nonempty
folders are rejected. Summary tables record the initial grouping threshold;
“Current groups CSV” exports the threshold currently selected in the HTML.

## Preview cleanup

The script targets known intermediate filenames under the specified results
root. It does not delete final assemblies, final GFA graphs, combined FASTAs,
QUAST/fastp reports, ANI tables/plots, detailed hit tables, native trees,
alignments, partitions, logs, manifests or dashboard bundles. It leaves unknown
files, symlinks, unfinished assembly runs and abandoned phylogeny staging
folders alone. It does not search outside the results root.

Start with a dry run:

```bash
python ~/bin/SprayNPray/spraynpray_cleanup.py --root .
```

The default `safe` profile selects numbered Unicycler intermediate graphs,
known assembler work directories and known Flye stage directories only when a
final FASTA, completion marker and final assembly graph are present.

For a more substantial reduction:

```bash
python ~/bin/SprayNPray/spraynpray_cleanup.py --root . --profile compact
```

`compact` also selects:

- `trimmed/{R1,R2,ONT}.fastq.gz` only for samples with a final assembly,
  completion marker and QUAST report. **Keep your original reads elsewhere.**
- SprayNPray `diamond.tsv`, input/protein/GFF intermediates only when
  `top_hits.csv`, `contigs.csv` and `run_manifest.json` are present. Detailed
  `matches.csv` and inference summaries remain.
- `cohort_ani/work/genomes` staging copies and a known reference download ZIP
  only when ANI tables and node metadata exist.
- `cohort_phylogen*/genomes` staging copies only when a native tree, alignment,
  completed GToTree log, selection CSV and labels file exist.

Downloaded reference FASTAs are retained unless you also pass
`--remove-reference-cache` with `compact`. Removing them may require future
reference downloads. Compact cleanup intentionally sacrifices cached inputs:
rerunning affected stages can require regenerating trimmed reads, ORFs, search
results or staging files. Preserve raw reads and databases separately.

The TSV plan lists each selected file and its logical size. The console also
estimates reclaimable disk blocks, accounting for hardlinks to retained files;
filesystem snapshots, reflinks and compression can make actual savings differ.
A directory listing alone cannot establish how much of 80 GB will be reclaimed.

## Apply the reviewed profile

After reviewing `cleanup_plan.tsv`, repeat the same options with `--apply`:

```bash
python ~/bin/SprayNPray/spraynpray_cleanup.py --root . --profile compact --apply
```

Each invocation builds a fresh plan. Run on an idle results directory. Before
unlinking, files are checked against their planned inode, size and modification
time; changed files are skipped. The script writes a JSON deletion audit and
removes only empty parent directories of files it deleted. Deletion is permanent;
there is no hidden trash folder retaining the disk usage.

To include the downloaded reference cache, preview and then apply explicitly:

```bash
python ~/bin/SprayNPray/spraynpray_cleanup.py --root . --profile compact --remove-reference-cache
python ~/bin/SprayNPray/spraynpray_cleanup.py --root . --profile compact --remove-reference-cache --apply
```

## Validation

```bash
python -m unittest discover -s tests -v
```

Tests cover cleanup boundaries, empty completion markers, symlinks, changed
files, retained hardlinks, dry-run/apply/audit, incomplete cohorts, Flye fallback,
ANI rejection, complete linkage, Newick parsing, native label mapping, portable
ZIP contents and HTML escaping. If Node.js is installed, the suite also executes
the generated dashboard JavaScript and checks filters, grouping, click details,
reference toggling, zoom and downloads against a small DOM harness. This checks
logic; it is not a browser layout test. Node.js is not a runtime dependency of
either Python script or the offline dashboard.
