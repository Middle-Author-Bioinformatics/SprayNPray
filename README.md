# SprayNPray

SprayNPray profiles genome and metagenome contigs. The assembly-backed V2
workflow runs **Prodigal → DIAMOND → taxonomic inference**, then optionally
combines samples for genome ANI analysis and an interactive reference map.
The original `spray-and-pray.py` workflow and examples are retained below.

| Task | Script | Detailed guide |
| --- | --- | --- |
| Start with a folder of genome FASTAs and run every step | `spraynpray_pipeline.py` | [Genome-folder pipeline](README_pipeline.md) |
| Profile a single assembly | `spraynpray_v2.py` | [V2 inputs, inference and outputs](README_v2.md) |
| Compare existing sample results and genomes | `spraynpray_compare.py` | [Cohort summary and ANI map](README_compare.md) |
| Build a GToTree phylogeny from samples and their closest references | `spraynpray_phylogeny.py` | [Phylogeny helper](README_phylogeny.md) |
| Refresh nearest-reference summaries and labels without rerunning ANI | `spraynpray_helper.py` | [Existing-result helper](README_compare.md#concise-nearest-reference-result-and-helper) |
| Archive intermediates and share an offline assembly/ANI/tree report | `spraynpray_cleanup.py`, `spraynpray_dashboard.py` | [Cleanup and portable dashboard](README_dashboard.md) |

## Rapid and simple taxonomic profiling of genome and metagenome contigs 
#### - Bin validation
#### - Metagenome profiling
#### - Contaminant identification and filtering
#### - Endosymbiont identification
#### - Identification of bacteria-to-eukaryote horizontal gene transfers

## Citing SprayNPray
If you found this software useful to your research please cite as follows:

Garber, A. I., Armbruster, C. R., Lee, S. E., Cooper, V. S., Bomberger, J. M., & McAllister, S. M. (2021). SprayNPray: user-friendly taxonomic profiling of genome and metagenome contigs. In bioRxiv (p. 2021.07.17.452725). https://doi.org/10.1101/2021.07.17.452725

Please also cite the various dependencies used by spray-and-pray: [DIAMOND](https://pubmed.ncbi.nlm.nih.gov/25402007/), [Prodigal](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2848648/), [Metabat2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6662567/), and [Biopython](https://biopython.org/)


## Installation for V2 and cohort analysis

Clone this repository and create the supplied dependency environment using
[mamba](https://mamba.readthedocs.io/) or conda:

```bash
git clone https://github.com/Middle-Author-Bioinformatics/SprayNPray.git
cd SprayNPray
conda config --set channel_priority strict
mamba env create -f environment.yml
conda activate spraynpray
```

With conda, replace the environment creation command with
`conda env create -f environment.yml`. To update an existing environment after
updating the checkout, run `mamba env update -n spraynpray -f environment.yml`.
The file specifies conda-forge before Bioconda and excludes default channels.
Run the following commands from this repository directory; keep the Python
scripts together. From another directory, use their full paths.

The environment installs Python 3.11, Prodigal, DIAMOND, FastANI, NCBI Datasets,
NumPy, SciPy and Matplotlib. Reference databases are supplied separately.
`setup.sh` is retained for the original workflow; use `environment.yml` for
the V2/cohort tools. The environment file installs dependencies; packaging
SprayNPray itself for `mamba install spraynpray` is a separate publishing step,
described in [CONDA_PACKAGING.md](CONDA_PACKAGING.md).

## Required reference inputs

Prepare these once and reuse them across samples:

- A DIAMOND `.dmnd` database built from reference proteins with accession IDs.
- `protein_to_genome.tsv`: two tab-separated columns, protein accession and
  versioned source assembly accession. Keep every distinct protein/assembly
  pair when an accession is shared between assemblies.
- `ncbi_assembly_info.tsv`: the original NCBI assembly summary, including its
  `#assembly_accession` header, taxids and organism names.
- An NCBI taxdump directory containing `nodes.dmp` and `names.dmp`, with
  `merged.dmp` when available.

The protein database and protein-to-genome map must describe the same
collection. See [V2 reference inputs](README_v2.md#reference-inputs) for formats
and database construction. Genome ANI also needs reference **nucleotide**
FASTAs: missing references are downloaded using NCBI Datasets, or supplied via
`--reference-dir /path/to/reference_genomes` with accession-named files such as
`GCF_000233495.1.fa`. Protein percent identity is not genome ANI.

## Quick start: a folder of genomes

Put one genome per file in a folder; `.fa`, `.fna`, `.fasta` and gzip-compressed
versions are accepted. The filename without its FASTA suffix becomes the
sample name. For example, `sample_104_ONT013.fa` becomes `sample_104_ONT013`.

```bash
python spraynpray_pipeline.py \
  --genomes /path/to/assemblies_combined \
  --diamond-db /path/to/pseudomonas.dmnd \
  --protein-to-genome /path/to/protein_to_genome.tsv \
  --assembly-info /path/to/ncbi_assembly_info.tsv \
  --taxdump /path/to/taxdump \
  --threads 32 \
  --outdir /path/to/spraynpray_pipeline_results
```

Prodigal uses single-genome mode by default; add `--meta` for metagenomic mode.
Only the **first 25 ORFs per contig**, ordered by genomic coordinate, are
searched by default. Change this with `--orfs-per-contig 50` or
`--orfs-per-contig all`. The contig report lists **one top protein hit per
selected ORF**, including explicit no-hit entries. A protein shared by multiple
assemblies has a representative source in the display; the full candidate set
is retained for taxonomic inference and detailed reporting.

The default ANI scope computes all sample/sample and reciprocal
sample/reference comparisons, skipping reference/reference pairs to save time.
Use `--ani-scope all` for the full union all-vs-all calculation. Optionally use
`--top-references-per-sample 5` to reduce the included reference set; the nearest
reference is then the closest among that restricted set. `--threads` controls
the total CPU budget; `--sample-jobs 4 --threads 32` permits four concurrent
V2 sample runs with eight DIAMOND threads each. See the
[pipeline guide](README_pipeline.md) for memory considerations and all options.

The map defaults to each sample's **5 closest references by accepted genome
ANI**, with up to 5 reference connections per sample. Shared references appear
once; every input sample remains visible. Change the display with
`--network-references-per-sample 20`, or use `0` for all reference nodes. This
controls the plot only and preserves the full ANI results and closest-reference
summary. It differs from the hit-count-based analysis limit above. Click any
dot in the HTML to view genome details, including available strain/isolate
metadata and ANI context.

Rerun with the same output directory to reuse compatible completed Prodigal,
DIAMOND and ANI stages. Use `--force-reblast` or `--force-ani` only when you
want to recompute those stages.

| Result | Location under the wrapper output directory |
| --- | --- |
| Per-contig reports and one-hit-per-ORF tables | `samples/SAMPLE/contigs.csv`, `samples/SAMPLE/top_hits.csv` |
| Per-sample reference hit percentages and mean protein identity | `cohort_ani/master_summary.csv` |
| One nearest reference and estimated overall genome ANI per sample | `cohort_ani/nearest_references.csv` |
| Pairwise genome ANI and reciprocal aligned fractions | `cohort_ani/ani_pairs.csv` |
| Interactive map with black samples, gray references and separated reference labels | `cohort_ani/ani_network.html` |

The nearest reference is selected by the highest accepted reciprocal genome
ANI, with distance `100 - ANI`. Map line lengths approximate those distances;
use the tables or hover values for measurements. Protein-hit percentages
summarize the representative sources of the sampled ORFs, not genome-wide
strain abundance. See the [comparison guide](README_compare.md) for acceptance
thresholds, ties, missing comparisons and map interpretation.

## Run individual stages or refresh existing results

For a single assembly, run the profiling stage directly:

```bash
python spraynpray_v2.py \
  --contigs sample_104_ONT013.fa \
  --diamond-db /path/to/pseudomonas.dmnd \
  --protein-to-genome /path/to/protein_to_genome.tsv \
  --assembly-info /path/to/ncbi_assembly_info.tsv \
  --taxdump /path/to/taxdump \
  --threads 32 \
  --outdir sample_104_ONT013_spraynpray
```

To compare existing results, place each sample's complete `top_hits.csv` in a
folder under its sample name, for example `sample_104_ONT013.csv`, with the
matching `sample_104_ONT013.fa` in the genome folder:

```bash
python spraynpray_compare.py \
  --results /path/to/spraynpray_combined \
  --genomes /path/to/assemblies_combined \
  --threads 32 \
  --outdir /path/to/cohort_ani
```

V2 `contigs.csv` files renamed by sample are also supported; complete
`top_hits.csv` files preserve every ORF and full numeric identity precision.
To add the concise summary and expand reference labels in an existing
comparison output, without recalculating ANI:

```bash
python spraynpray_helper.py --analysis-dir /path/to/cohort_ani
```

For older outputs, add `--assembly-info /path/to/ncbi_assembly_info.tsv` to
enrich reference strain/isolate details in the click panel. New wrapper runs
pass this metadata automatically.

The helper also applies the top-5 reference display without recalculating ANI;
pass `--network-references-per-sample 0` to restore all reference nodes.
Add `--summary-only` to the helper for just `nearest_references.csv`.
Each script documents its flags with `--help`.

## Build a phylogenomic tree

The optional [GToTree helper](README_phylogeny.md) uses every input sample and
its **five closest accepted ANI references**, including shared references once:

```bash
# Run in an environment containing GToTree.
python spraynpray_phylogeny.py \
  --analysis-dir /path/to/cohort_ani \
  --outdir /path/to/cohort_phylogeny \
  --hmm-set Gammaproteobacteria \
  --references-per-sample 5 \
  --threads 24
```

Use `Bacteria` (the default) or another marker set if Gammaproteobacteria does
not fit the full cohort. The helper reuses local genome FASTAs and runs no ANI
or DIAMOND. It keeps accession/organism/strain labels when available, retains
GToTree's alignment and reports, and records genomes excluded by its filtering.
See the [phylogeny guide](README_phylogeny.md) for installation, nucleotide mode,
reference limits, preparation-only mode and completed-tree reuse.

## Original SprayNPray workflow

The following installation, profiling and contig-filtering examples use the
original `spray-and-pray.py`. Its `sprayandpray` environment and database setup
are separate from the V2/cohort environment above.


Provide either raw contigs or ORFs (amino acid) in FASTA format.
A reference protein dataset also needs to be provided. Ideally, this should be NCBI's RefSeq or nr database from NCBI, 
which can be downloaded via (or follow the install instructions):
    
    wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz

This program will then DIAMOND BLAST each ORF against the reference database
(This takes about 30 minutes for a typical bacteria genome of ~4Mbp).

The program will then output a file in CSV format summarizing the dominant taxanomic groups matched to each contig.
This allows the user to visually inspect the data, seeing what the closest taxonomic group to each contig is.
Additionally, the user can specify a set of criteria (e.g. GC-content, read coverage, coding density, closest taxonomic hits) to re-write the provided contigs into a new FASTA file.


### Original installation with conda:
  
    git clone https://github.com/Middle-Author-Bioinformatics/SprayNPray.git
    cd SprayNPray
    bash setup.sh
    bash dlDB.sh /path/to/preferred/directory/for/db
    conda activate sprayandpray

Please also place spray-and-pray.py into the PATH variable: e.g. ~/bin/SprayNPray. You can update the PATH variable in a file called '.bash_profile' or '.profile' in your home directory.

Do not worry about the dependencies after conda installation. Just enter `conda deactivate` when finished using the program.


### Original installation without conda:

    git clone https://github.com/Middle-Author-Bioinformatics/SprayNPray.git
    export PATH=$PATH:$(pwd)

(be sure to put spray-and-pray.py into your $PATH in your bash profile)

### Install Dependencies:

* [DIAMOND](https://github.com/bbuchfink/diamond)
* [Prodigal](https://github.com/hyattpd/Prodigal)
* [Metabat](https://bitbucket.org/berkeleylab/metabat)
* [Python3](https://www.python.org/download/releases/3.0/)
* [Biopython3](https://biopython.org/)
* [Joblib](https://joblib.readthedocs.io/en/latest/)

### Test run to check dependencies or to check conda environment:

    bash testRun.sh

### Original usage

### quick-start

    spray-and-pray.py -g genomeContigs.fa -out genomeContigs -ref /path/to/nr/nr.faa

### Decontaminating a Pseudomonas assembly

    spray-and-pray.py -a pseudomonas_crude.fa -out pseudomonas_clean.fa -Genus Pseudomonas -species aeruginosa -perc 50 --fa -ref /path/to/nr/nr.faa

In the above command, we require that at least 50% of contig's genes have a top DIAMOND hit to Pseudomonas aruginosa. Contigs matching these parameters will be written to a new FASTA file: pseudomonas_clean.fa


### Pulling out endosymbiont genomes from an assembly of the mealybug Maconellicoccus hirsutus

    spray-and-pray.py -a M_hirsutus_assembly.fa -out endosymbionts.fa -cd 0.5 -L 1000000 -perc 50 --fa -ref /path/to/nr/nr.faa -Domain Bacteria

In the above comomand, we require a gene density of 0.5 genes per kb, maximum length of 1 Mb, top DIAMOND hit to be to a bacterial gene, and that at least 50% of the contig's genes to be of bacterial origin to be written to endosymbionts.fa


### Identifying putative HGTs in an assembly of the mealybug Planococcus citri

    spray-and-pray.py -a P_citri.fa -out putative_hgts.csv --hgt -ref /path/to/nr/nr.faa


