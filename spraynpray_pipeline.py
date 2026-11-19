#!/usr/bin/env python3
"""Genome folder -> per-sample SprayNPray V2 -> cohort summary and ANI map.

Keep this wrapper beside spraynpray_v2.py and spraynpray_compare.py. Existing
completed gene prediction, alignment, reference download and ANI stages are
reused by those tools. Default sample processing is sequential.
"""

from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import spraynpray_compare as compare

LOG = logging.getLogger("spraynpray_pipeline")
ROOT = Path(__file__).resolve().parent


def discover_genomes(folder):
    samples = {}
    for path in sorted(folder.iterdir()):
        suffix = next((s for s in sorted(compare.FASTA_SUFFIXES, key=len, reverse=True) if path.name.endswith(s)), None)
        if not path.is_file() or suffix is None:
            continue
        name = path.name[:-len(suffix)]
        if not name or name in (".", "..") or name in samples:
            raise ValueError("Empty or duplicate sample ID in genome filenames: " + name)
        samples[name] = path.resolve()
    if not samples:
        raise ValueError("No .fa/.fna/.fasta genomes (optionally .gz) in " + str(folder))
    return samples


def run_sample(sample, source, args, threads):
    genome = args.outdir / "assemblies_combined" / (sample + ".fa")
    compare.stage_fasta(source, genome)
    output = args.outdir / "samples" / sample
    command = [sys.executable, str(ROOT / "spraynpray_v2.py"), "--contigs", str(genome),
               "--diamond-db", str(args.diamond_db.resolve()), "--protein-to-genome", str(args.protein_to_genome.resolve()),
               "--assembly-info", str(args.assembly_info.resolve()), "--taxdump", str(args.taxdump.resolve()),
               "--outdir", str(output), "--threads", str(threads), "--orfs-per-contig", args.orfs_per_contig]
    for name in ("genetic_code", "evalue", "min_identity", "query_cover", "subject_cover", "top_percent", "sensitivity",
                 "consensus", "min_tax_orfs", "min_classified_fraction", "max_resolution", "prodigal", "diamond"):
        command.extend(["--" + name.replace("_", "-"), str(getattr(args, name))])
    if args.meta:
        command.append("--meta")
    if args.force_reblast:
        command.append("--force-reblast")
    log_path = args.outdir / "logs" / (sample + ".log")
    LOG.info("SprayNPray: %s (%s DIAMOND threads)", sample, threads)
    compare.run_command(command, log_path)
    top_hits = output / "top_hits.csv"
    # Validate before replacing this sample's cohort input. Top hits retain
    # complete ORF coverage and original identity precision, even in all mode.
    compare.load_sample(top_hits)
    target = args.outdir / "spraynpray_combined" / (sample + ".csv")
    temporary = target.with_name(target.name + ".tmp")
    shutil.copyfile(top_hits, temporary)
    temporary.replace(target)
    return dict(sample=sample, source_file=str(source), staged_genome=str(genome), result_csv=str(target), status="complete")


def compare_command(args):
    command = [sys.executable, str(ROOT / "spraynpray_compare.py"),
               "--results", str(args.outdir / "spraynpray_combined"),
               "--genomes", str(args.outdir / "assemblies_combined"),
               "--outdir", str(args.outdir / "cohort_ani"), "--assembly-info", str(args.assembly_info.resolve())]
    for name in ("threads", "fastani", "datasets", "frag_len", "min_af", "min_ani", "neighbors", "ani_scope", "top_references_per_sample", "network_references_per_sample"):
        command.extend(["--" + name.replace("_", "-"), str(getattr(args, name))])
    if args.reference_dir is not None:
        command.extend(["--reference-dir", str(args.reference_dir.resolve())])
    for name in ("no_download", "allow_missing_references", "force_ani"):
        if getattr(args, name):
            command.append("--" + name.replace("_", "-"))
    return command


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--genomes", type=Path, required=True, help="Folder of one genome FASTA per sample")
    p.add_argument("--diamond-db", type=Path, required=True)
    p.add_argument("--protein-to-genome", type=Path, required=True)
    p.add_argument("--assembly-info", type=Path, required=True)
    p.add_argument("--taxdump", type=Path, required=True)
    p.add_argument("-o", "--outdir", type=Path, required=True)
    p.add_argument("--threads", type=int, default=32, help="Total CPU budget; split across concurrent sample jobs")
    p.add_argument("--sample-jobs", type=int, default=1, help="Concurrent per-sample V2 jobs (each also loads its reference/taxonomy data)")
    p.add_argument("--meta", action="store_true")
    p.add_argument("--orfs-per-contig", default="25", help="Positive integer or all")
    p.add_argument("--genetic-code", type=int, default=11)
    p.add_argument("--evalue", type=float, default=1e-15)
    p.add_argument("--min-identity", type=float, default=0)
    p.add_argument("--query-cover", type=float, default=50)
    p.add_argument("--subject-cover", type=float, default=50)
    p.add_argument("--top-percent", type=float, default=5)
    p.add_argument("--sensitivity", choices=("default", "fast", "mid-sensitive", "sensitive", "very-sensitive", "ultra-sensitive"), default="sensitive")
    p.add_argument("--consensus", type=float, default=.8)
    p.add_argument("--min-tax-orfs", type=int, default=3)
    p.add_argument("--min-classified-fraction", type=float, default=.5)
    p.add_argument("--max-resolution", choices=("species", "taxid"), default="species")
    p.add_argument("--prodigal", default="prodigal")
    p.add_argument("--diamond", default="diamond")
    p.add_argument("--force-reblast", action="store_true")
    p.add_argument("--reference-dir", type=Path)
    p.add_argument("--no-download", action="store_true")
    p.add_argument("--allow-missing-references", action="store_true")
    p.add_argument("--top-references-per-sample", type=int, default=0)
    p.add_argument("--fastani", default="fastANI")
    p.add_argument("--datasets", default="datasets")
    p.add_argument("--frag-len", type=int, default=3000)
    p.add_argument("--min-af", type=float, default=.5)
    p.add_argument("--min-ani", type=float, default=80)
    p.add_argument("--network-references-per-sample", type=int, default=5,
                   help="Top references shown per sample ranked by accepted genome ANI; 0 shows all reference nodes")
    p.add_argument("--neighbors", type=int, default=5)
    p.add_argument("--ani-scope", choices=("samples", "all"), default="samples")
    p.add_argument("--force-ani", action="store_true")
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    if not args.genomes.is_dir() or args.threads < 1 or args.sample_jobs < 1:
        p.error("Provide an existing genome folder and positive thread/job counts")
    if args.sample_jobs > args.threads:
        p.error("--sample-jobs cannot exceed the total --threads budget")
    if args.reference_dir is not None and not args.reference_dir.is_dir():
        p.error("--reference-dir must be an existing directory")
    if not 0 <= args.min_af <= 1 or not 0 <= args.min_ani <= 100 or args.frag_len < 100 or args.neighbors < 0 or args.top_references_per_sample < 0 or args.network_references_per_sample < 0:
        p.error("Invalid ANI threshold, fragment length, neighbor count or reference limit")
    for path in (args.protein_to_genome, args.assembly_info, args.taxdump / "nodes.dmp", args.taxdump / "names.dmp"):
        if not path.is_file():
            p.error("File not found: " + str(path))
    if not args.diamond_db.is_file() and not Path(str(args.diamond_db) + ".dmnd").is_file():
        p.error("DIAMOND database not found")
    try:
        samples = discover_genomes(args.genomes)
    except ValueError as exc:
        p.error(str(exc))
    args.outdir = args.outdir.resolve()
    args.outdir.mkdir(parents=True, exist_ok=True)
    for name in ("samples", "assemblies_combined", "spraynpray_combined", "logs"):
        (args.outdir / name).mkdir(exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(args.outdir / "pipeline.log")])
    manifest_path = args.outdir / "pipeline_manifest.json"
    old = compare.read_json(manifest_path)
    previous_ids = set(old.get("sample_ids", []))
    if any(not isinstance(n, str) or not n or n in (".", "..") or Path(n).name != n for n in previous_ids):
        p.error("Invalid sample IDs in the existing pipeline manifest")
    combined = args.outdir / "spraynpray_combined"
    existing_ids = {path.stem for path in combined.glob("*.csv")}
    unknown = existing_ids - set(samples) - previous_ids
    if unknown:
        p.error("Unrecognized CSV files in wrapper output folder: " + ", ".join(sorted(unknown)))
    # Remove only obsolete cohort inputs recorded as this wrapper's outputs.
    # Historical per-sample analysis folders and staged FASTAs are retained.
    for name in previous_ids - set(samples):
        (combined / (name + ".csv")).unlink(missing_ok=True)
    manifest = dict(version=compare.VERSION, status="running", sample_ids=list(samples), samples=[],
                    arguments={k: str(v.resolve()) if isinstance(v, Path) else v for k, v in vars(args).items()})
    compare.write_json(manifest_path, manifest)
    jobs = min(args.sample_jobs, len(samples))
    per_sample_threads = max(1, args.threads // jobs)
    try:
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = {pool.submit(run_sample, name, source, args, per_sample_threads): name for name, source in samples.items()}
            failures = []
            for future in as_completed(futures):
                name = futures[future]
                try:
                    manifest["samples"].append(future.result())
                except (ValueError, OSError, subprocess.CalledProcessError) as exc:
                    failures.append(name)
                    manifest["samples"].append(dict(sample=name, status="failed", error=str(exc)))
                    LOG.error("Sample %s failed; see logs/%s.log", name, name)
                compare.write_json(manifest_path, manifest)
        if failures:
            raise ValueError("Sample analysis failed for: " + ", ".join(sorted(failures)) + ". Cohort comparison was not run; repeat this command to resume.")
        manifest["status"] = "comparing"
        compare.write_json(manifest_path, manifest)
        compare.run_command(compare_command(args), args.outdir / "logs" / "compare.log")
        manifest["status"] = "complete"
        manifest["nearest_references"] = str(args.outdir / "cohort_ani" / "nearest_references.csv")
        compare.write_json(manifest_path, manifest)
        LOG.info("Finished: %s", args.outdir / "cohort_ani" / "ani_network.html")
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        manifest["status"] = "failed"
        manifest["error"] = str(exc)
        compare.write_json(manifest_path, manifest)
        LOG.error("%s", exc)
        LOG.error("See pipeline.log and logs/*.log. Rerun the same command to resume completed stages.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
