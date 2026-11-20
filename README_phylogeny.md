# Phylogenomic trees from a completed ANI analysis

`spraynpray_phylogeny.py` takes a completed `spraynpray_compare` output folder
and submits **every input sample plus its five closest accepted references**
to [GToTree](https://github.com/AstrobioMike/GToTree). References are ranked by
mean reciprocal genome ANI, with higher minimum aligned fraction and then
accession breaking ties. Shared references are included once. This matches
the network's new default top-five reference selection; it is a per-sample
limit, so a cohort can still have more than five unique references.

The helper uses local nucleotide FASTAs already used by the comparison. It
runs no ANI, DIAMOND or reference downloads. GToTree identifies marker genes
across the supplied genomes independently of the first-25-ORF SprayNPray
sampling, aligns/trims them, and builds a concatenated-marker phylogeny.

## Installation

Keep this script beside the updated `spraynpray_compare.py`. Run it in an
existing GToTree environment, or install GToTree separately using its
[author's installation instructions](https://github.com/AstrobioMike/GToTree/wiki/installation):

```bash
conda create -n gtotree -c astrobiomike -c conda-forge -c bioconda gtotree
conda activate gtotree
```

GToTree is optional and is not added to the main SprayNPray environment. Use
`gtt-hmms` to list the marker sets available in your GToTree installation.

## Run for the Pseudomonas cohort

From the SprayNPray checkout, with the GToTree environment active:

```bash
python spraynpray_phylogeny.py \
  --analysis-dir ~/MAB/PITT/Rapsinski/long/hybrid_results/cohort_ani \
  --outdir ~/MAB/PITT/Rapsinski/long/hybrid_results/cohort_phylogeny \
  --hmm-set Gammaproteobacteria \
  --references-per-sample 5 \
  --threads 24
```

For wrapper outputs, point `--analysis-dir` to the wrapper's `cohort_ani`
subfolder. At least three genomes must be selected. The default marker set is
`Bacteria`; use `Gammaproteobacteria` when that set fits the entire cohort, as
in a Pseudomonas-focused analysis. A custom HMM path is also accepted. The
output directory and custom HMM path must have no whitespace because of
GToTree's shell workflow. Original source FASTA paths can contain spaces: the
helper stages them under safe, unique filenames.

The default is an amino-acid marker tree. Add `--nucleotide` for GToTree's
nucleotide marker alignment/tree (`-z`), an option to consider for closely
related strains. Neither tree is an ANI-distance network or a recombination-
filtered variant phylogeny. Default GToTree marker/genome filtering is retained;
inspect its marker-completeness reports and `tree_membership.csv` for genomes
excluded during filtering, and inspect the alignment/tree before interpretation.
No rooting/outgroup is imposed by this helper.

`--tree-program IQ-TREE` selects GToTree's IQ-TREE mode; otherwise GToTree
chooses its default engine. `--threads` is passed as parallel jobs (`-j`), with
HMMER (`-n`) and MUSCLE (`-M`) limited to one thread per job. Tree-engine
parallelism follows GToTree's own implementation.

## Genome sources, labels and reuse

The helper preferentially uses normalized cached FASTAs whose checksums match
the saved ANI cache, then the original `source_file` paths or cached reference
files. If originals moved, provide replacement folders:

```bash
--genomes /path/to/assemblies_combined \
--reference-dir /path/to/reference_genomes
```

Sample files must use their sample names and a supported FASTA suffix;
reference files must contain their exact versioned assembly accessions, as in
the comparison tool. Replacements must match the original sequences when ANI
cache checksums are available. Missing or changed selected genomes stop the
run with an explicit message rather than being silently omitted.

Reference tree labels include their accession, organism, strain and isolate
when saved metadata provides those values; sample labels preserve their IDs
with a sample prefix. Labels are sanitized to safe Newick characters and repeated
underscores are collapsed, matching GToTree's final label normalization. Use
`--assembly-info /path/to/ncbi_assembly_info.tsv` to enrich older outputs before
building labels. This does not assign reference strain names to input samples.

`--references-per-sample N` changes the tree's reference selection; `0` includes
every reference with an accepted sample/reference ANI comparison. Samples with
no accepted reference are still submitted. This tree option does not redraw
the network; the display's separate control is
`--network-references-per-sample N` on the comparison/helper/wrapper.

`--prepare-only` stages the selected FASTAs and writes a GToTree command script
without requiring GToTree to be installed or running it. Repeating the normal
helper command reuses a compatible completed tree whose checksum still matches.
New inputs/settings or `--force` create a new attempt under `runs/`, retaining
previous attempt outputs and logs. Use `--force` after upgrading GToTree or
replacing a named built-in marker set: these installed resource versions are
not automatically fingerprinted. Custom HMM files are checksummed.

## Recover an already finished GToTree run

GToTree renames its final tree to the output-directory basename, for example
`runs/attempt_001/attempt_001.tre`. Earlier versions of this helper looked only
for the intermediate `Aligned_SCGs*.tre` names and could incorrectly mark a
successful native run as failed. The corrected helper recognizes both naming
styles and checks existing attempts before launching another run.

After updating the code, repeat your original settings with **the same output
folder** and add `--reuse-only` to guarantee GToTree will not be launched:

```bash
python spraynpray_phylogeny.py \
  --analysis-dir /path/to/cohort_ani \
  --outdir /path/to/existing/cohort_phylogeny \
  --hmm-set Gammaproteobacteria \
  --references-per-sample 5 \
  --threads 24 \
  --reuse-only
```

GToTree also collapses repeated underscores in labels (`sample__name` becomes
`sample_name`). Earlier helper versions generated double underscores and then
rejected these valid normalized output labels. Recovery and membership checks
now accept that exact normalization, and new runs generate normalized labels.
Original manifests with double-underscore labels remain compatible. Sequence
checksums, genome IDs and marker/tree settings must still match; arbitrary
label differences and normalization collisions are not accepted.

Recovery requires matching genome/settings signatures in an existing attempt
manifest, a native GToTree completion log, a nonempty alignment and a final
tree whose leaves match the submitted labels and alignment after this normalization. The helper repairs
the completion manifest and regenerates `tree_membership.csv`; the original
tree/alignment are preserved. A later incomplete attempt does not prevent
recovery of an earlier compatible finished one. Missing, partial or incompatible
outputs produce an error with `--reuse-only`, rather than a new analysis.
Rejection warnings now identify the specific mismatched settings or missing/
incompatible artifact and print the absolute output path being checked. Confirm
that `--outdir` points to the folder holding the original completed attempt;
a similarly named folder elsewhere will not contain that run.
Normal runs also attempt this recovery automatically; `--force` deliberately
bypasses it.

## Outputs

| Output | Contents |
| --- | --- |
| `selected_genomes.csv` | Genome IDs, sample/reference status, metadata, tree labels, source paths and normalized sequence hashes. |
| `reference_selection.csv` | Each sample's selected references, rank, ANI and aligned fraction. |
| `genomes/`, `genomes.txt` | Staged DNA FASTAs and the GToTree input list. |
| `labels.tsv` | Headerless GToTree mapping from staged FASTA basenames to informative tree labels. |
| `run_gtotree.sh` | Exact GToTree command, run from this output directory. |
| `runs/attempt_NNN/attempt_NNN.tre` | Final informatively labeled Newick tree. Older `Aligned_SCGs*.tre` output names are also recognized. |
| `runs/attempt_NNN/` | GToTree alignment, marker/filtering reports, citations and other outputs. |
| `tree_membership.csv` | Whether each submitted genome appears in the final tree; filtered genomes are explicitly listed. |
| `gtotree_attempt_NNN.log` | Full stdout/stderr for that attempt. |
| `run_manifest.json`, `selection_manifest.json`, `runs/*.manifest.json` | Settings, genome selection, completion/cache status and attempt history. |

Open the Newick tree in a tree viewer such as FigTree or iTOL. Cite GToTree and
the dependencies listed in its generated `citations.txt` for the actual run.

GToTree flag/label/output handling was checked against its
[upstream implementation](https://github.com/AstrobioMike/GToTree/blob/f93319ca0702efa13620e6cbf1bb912ff0181a1a/bin/GToTree).
The tests exercise selection/staging, GToTree command construction with a
simulated runner, output validation, resume/failure handling and membership
reporting. A full GToTree marker search/tree run requires its separately
installed toolchain and marker sets.
