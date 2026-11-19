# One-command genome-folder analysis

`spraynpray_pipeline.py` runs the existing resumable tools in sequence:

1. Normalize one genome FASTA per sample.
2. Run `spraynpray_v2.py` per sample (Prodigal, first-25-ORF selection, DIAMOND,
   top hits, assembly metadata and taxonomy).
3. Collect each sample's complete `top_hits.csv` under its sample name.
4. Run `spraynpray_compare.py` for hit summaries, reference genome acquisition,
   reciprocal genome ANI, `nearest_references.csv` and the expanded-label map.

## Dependencies

Use the same Python environment for the wrapper and both underlying scripts.
Keep `spraynpray_pipeline.py`, `spraynpray_v2.py` and `spraynpray_compare.py`
beside one another in the repository. The optional existing-result helper is
`spraynpray_helper.py`.

```bash
# From the repository checkout; see README.md for cloning instructions.
conda config --set channel_priority strict
mamba env create -f environment.yml
conda activate spraynpray
```

## Run

```bash
python spraynpray_pipeline.py \
  --genomes /path/to/assemblies_combined \
  --diamond-db /path/to/reference.dmnd \
  --protein-to-genome /path/to/protein_to_genome.tsv \
  --assembly-info /path/to/ncbi_assembly_info.tsv \
  --taxdump /path/to/taxdump \
  --threads 32 \
  --outdir /path/to/spraynpray_pipeline_results
```

Accepted input suffixes: `.fa`, `.fna`, `.fasta`, and their gzip-compressed
versions. The full suffix is removed to obtain the sample ID; duplicate IDs
such as `sample.fa` and `sample.fna.gz` are rejected. Inputs are discovered
directly in the folder rather than recursively.

Prodigal defaults to single-genome mode. Add `--meta` for metagenomic mode.
`--orfs-per-contig` defaults to `25` and also accepts a positive integer or
`all`. The wrapper exposes the V2 alignment/inference filters and comparison
controls; see `--help`. It collects complete top-hit tables, preserving numeric
identity precision and avoiding truncated contig-summary fields in all-ORF mode.

The default ANI scope skips reference/reference work, retaining all sample
pairs and reciprocal sample/reference comparisons. Use `--ani-scope all` if
reference/reference ANI is required. Reference downloads can be replaced with
local genomes using `--reference-dir`; `--no-download` and
`--allow-missing-references` are passed through. The optional
`--top-references-per-sample N` limits references for a smaller first run.

The map shows up to **10 closest reference connections per sample**, ranked by
accepted reciprocal genome ANI. Shared references are drawn once; all input
samples remain visible. `--network-references-per-sample N` changes this display
limit; `0` restores all reference nodes. This does not restrict ANI computation,
tables, or the nearest-reference summary. The separate
`--top-references-per-sample` option above restricts analysis inputs by protein
hit count. The interactive HTML supports click-to-select genome details, including
strain/isolate metadata when present. The wrapper passes `--assembly-info` to
the comparison stage automatically; input samples are distinguished from the
metadata of their closest reference.

## CPU allocation and resuming

`--threads` is the total CPU budget. Per-sample V2 jobs run sequentially by
default. `--sample-jobs 4 --threads 32` runs up to four sample jobs with eight
DIAMOND threads apiece, followed by the cohort ANI stage with 32 threads.
Prodigal itself remains single-threaded. Concurrent V2 jobs each load their
reference/taxonomy resources, so the default avoids multiplying those memory
requirements. The number of sample jobs cannot exceed the thread budget.

Rerun the same command with the same output directory to resume. The underlying
V2 code reuses compatible completed gene predictions and DIAMOND output, and
the comparison code reuses reference downloads and completed ANI. Metadata,
summary and plot outputs are refreshed. Use `--force-reblast` or `--force-ani`
to deliberately recompute those stages.

A sample failure stops the cohort comparison after current sample jobs finish;
the logs and successful sample outputs are retained for a retry. Cohort inputs
are copied only after the sample's complete top-hit table is validated. When
input samples are removed, their obsolete combined CSVs are removed only if
recorded as wrapper outputs; historical sample analysis folders remain intact.
Unrecognized extra CSVs in the combined output folder are rejected.

## Output layout

| Location | Contents |
| --- | --- |
| `samples/SAMPLE/` | Complete SprayNPray V2 results and prediction/alignment caches for each sample. |
| `assemblies_combined/` | Normalized, sample-named genome FASTAs used for the cohort analysis. |
| `spraynpray_combined/` | One complete top-hit CSV per active sample. |
| `cohort_ani/nearest_references.csv` | Concise closest-reference and overall ANI summary. |
| `cohort_ani/ani_network.html` | Interactive map with black samples, gray references and separated reference labels. |
| `cohort_ani/` | All other summary tables, ANI measurements/matrices, plots and reference/ANI caches. |
| `logs/`, `pipeline.log` | Sample and cohort command logs. |
| `pipeline_manifest.json` | Active sample list, completion/failure status and run arguments. |

For an existing cohort analysis, use the helper directly rather than starting
a new wrapper run:

```bash
python spraynpray_helper.py --analysis-dir /path/to/existing/cohort_ani
```

This updates the concise summary and reference labels without recalculating ANI.
Nearest-reference and graph interpretation are documented in
[README_compare.md](README_compare.md).
