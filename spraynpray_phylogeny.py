#!/usr/bin/env python3
"""Build a GToTree phylogeny from ANI results: all samples plus their closest references.

Keep this helper beside spraynpray_compare.py. No ANI, DIAMOND, or genome
reference downloads are performed. Requires a completed cohort comparison.
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import shlex
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import spraynpray_compare as compare

LOG = logging.getLogger("spraynpray_phylogeny")


def select_genomes(nodes, pairs, limit):
    if limit < 0:
        raise ValueError("Reference limit must be nonnegative")
    kinds = {n["id"]: n["kind"] for n in nodes}
    if len(kinds) != len(nodes) or not nodes or any(k not in ("sample", "reference") for k in kinds.values()):
        raise ValueError("Expected unique sample/reference IDs in network_nodes.csv")
    candidates = defaultdict(list)
    seen = set()
    for pair in pairs:
        a, b = pair["genome_a"], pair["genome_b"]
        if a not in kinds or b not in kinds:
            raise ValueError("ANI pair contains an unknown genome")
        key = tuple(sorted((a, b)))
        if key in seen:
            raise ValueError("Duplicate ANI pair: " + str(key))
        seen.add(key)
        if pair["status"] != "accepted" or {kinds[a], kinds[b]} != {"sample", "reference"}:
            continue
        sample, reference = (a, b) if kinds[a] == "sample" else (b, a)
        ani, af = compare.identity(pair["ani_used"]), float(pair["min_aligned_fraction"])
        if not 0 <= af <= 1:
            raise ValueError("Invalid aligned fraction")
        candidates[sample].append((reference, ani, af))
    rankings, references = [], set()
    samples = [n for n in nodes if n["kind"] == "sample"]
    if not samples:
        raise ValueError("No input samples in comparison")
    for node in samples:
        sample = node["id"]
        ordered = sorted(candidates[sample], key=lambda value: (-value[1], -value[2], value[0]))
        if not ordered:
            LOG.warning("%s has no accepted reference ANI; retaining the sample without references", sample)
        for rank, (reference, ani, af) in enumerate(ordered[:limit] if limit else ordered, 1):
            references.add(reference)
            rankings.append(dict(sample=sample, rank=rank, reference=reference, ani=ani, min_aligned_fraction=af))
    return [n for n in nodes if n["kind"] == "sample" or n["id"] in references], rankings


def locate_fasta(node, analysis, genomes=None, reference_dir=None):
    """Use original exact paths, or user-specified replacements after moving a run."""
    folder = genomes if node["kind"] == "sample" else reference_dir
    if folder is not None:
        if node["kind"] == "reference":
            path = compare.reference_files(folder).get(node["id"])
            if path:
                return path
        else:
            paths = [folder / (node["id"] + suffix) for suffix in compare.FASTA_SUFFIXES]
            matches = [p for p in paths if p.is_file()]
            if len(matches) > 1:
                raise ValueError("Multiple sample FASTAs for " + node["id"])
            if matches:
                return matches[0].resolve()
        raise ValueError("No replacement FASTA for " + node["id"] + " in " + str(folder))
    # Old network tables lack staged_file, but prepare_nodes used this stable
    # index. Only trust that fallback when its checksum agrees with the ANI cache.
    if node.get("ani_sha256"):
        cached = analysis / "work" / "genomes" / ("g%06d.fa" % node["cached_index"])
        if cached.is_file() and compare.sha256(cached) == node["ani_sha256"]:
            return cached.resolve()
    source = node.get("source_file", "")
    if source:
        path = Path(source)
        path = path if path.is_absolute() else analysis / path
        if path.is_file():
            return path.resolve()
    if node["kind"] == "reference":
        path = compare.reference_files(analysis / "work" / "references").get(node["id"])
        if path:
            return path
    raise ValueError("FASTA unavailable for %s; provide --genomes or --reference-dir with the exact original genomes" % node["id"])


def label_for(node):
    parts = [node["kind"], node["id"]]
    if node["kind"] == "reference":
        parts.extend(node.get(k, "") for k in ("organism", "strain", "isolate"))
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", "__".join(p for p in parts if p)).strip("_")


def prepare_inputs(args):
    _, nodes = compare.read_csv(args.analysis_dir / "network_nodes.csv")
    _, pairs = compare.read_csv(args.analysis_dir / "ani_pairs.csv")
    cache = compare.read_json(args.analysis_dir / "work" / "ani_cache.json")
    digests = {r["id"]: r["sha256"] for r in cache.get("signature", {}).get("genomes", [])}
    for index, node in enumerate(nodes):
        node.update(cached_index=index, ani_sha256=digests.get(node["id"]))
    compare.enrich_reference_nodes(nodes, args.assembly_info)
    selected, rankings = select_genomes(nodes, pairs, args.references_per_sample)
    if len(selected) < 3:
        raise ValueError("At least three selected genomes are required for this phylogeny helper")
    # Resolve all sources before writing inputs, so missing genomes are never silently skipped.
    sources = [locate_fasta(n, args.analysis_dir, args.genomes, args.reference_dir) for n in selected]
    inputs = args.outdir / "genomes"
    inputs.mkdir(parents=True, exist_ok=True)
    labels, records, used_labels = [], [], set()
    for index, (node, source) in enumerate(zip(selected, sources)):
        filename = "g%06d.fa" % index
        target = inputs / filename
        stats = compare.stage_fasta(source, target)
        if node.get("ani_sha256") and node["ani_sha256"] != stats["sha256"]:
            raise ValueError("Genome changed since ANI for %s; supply the original FASTA or recompute ANI" % node["id"])
        label = label_for(node)
        if label in used_labels:
            raise ValueError("Sanitized tree labels collide; rename the input sample: " + node["id"])
        used_labels.add(label)
        labels.append((filename, label))
        records.append(dict(node, tree_label=label, input_fasta="genomes/" + filename,
                            original_fasta=str(source), sha256=stats["sha256"]))
    (args.outdir / "genomes.txt").write_text("".join(r["input_fasta"] + "\n" for r in records))
    with (args.outdir / "labels.tsv").open("w", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(labels)
    compare.write_csv(args.outdir / "selected_genomes.csv",
                      ["id", "kind", "organism", "strain", "isolate", "tree_label", "original_fasta", "input_fasta", "sha256"], records)
    compare.write_csv(args.outdir / "reference_selection.csv",
                      ["sample", "rank", "reference", "ani", "min_aligned_fraction"], rankings)
    return records


def report_tree_membership(tree, records, outdir):
    # Labels generated here contain only safe unquoted Newick characters.
    leaves = set(re.findall(r"[,(]\s*([A-Za-z0-9_.-]+)\s*(?=[:),;])", tree.read_text()))
    if not leaves.intersection(r["tree_label"] for r in records):
        raise ValueError("Tree labels do not match the submitted genome labels")
    rows = [dict(id=r["id"], kind=r["kind"], tree_label=r["tree_label"],
                 present_in_tree=r["tree_label"] in leaves) for r in records]
    compare.write_csv(outdir / "tree_membership.csv", ["id", "kind", "tree_label", "present_in_tree"], rows)
    missing = [r["id"] for r in rows if not r["present_in_tree"]]
    if missing:
        LOG.warning("GToTree excluded %s selected genomes; inspect tree_membership.csv and its marker/filtering reports", len(missing))
    return missing


def find_tree(output):
    """GToTree renames the final tree to the output directory's basename."""
    for name in (output.name + ".tre", "Aligned_SCGs_mod_names.tre", "Aligned_SCGs.tre"):
        path = output / name
        if path.is_file():
            text = path.read_text().strip()
            if "(" in text and text.endswith(";"):
                return path
    return None


def completed_attempt(output, records):
    """Require native completion, an alignment, and compatible final leaf labels."""
    tree = find_tree(output)
    log = output / "gtotree-runlog.txt"
    if tree is None or not log.is_file():
        return None
    text = log.read_text()
    if "Done!!" not in text or "Total process runtime:" not in text:
        return None
    alignment = next((output / name for name in ("Aligned_SCGs_mod_names.faa", "Aligned_SCGs_mod_names.fa",
                     "Aligned_SCGs.faa", "Aligned_SCGs.fa")
                     if (output / name).is_file() and (output / name).stat().st_size), None)
    if alignment is None:
        return None
    leaves = set(re.findall(r"[,(]\s*([A-Za-z0-9_.-]+)\s*(?=[:),;])", tree.read_text()))
    with alignment.open() as handle:
        headers = {line[1:].strip().split()[0] for line in handle if line.startswith(">") and line[1:].strip()}
    if not leaves or not leaves.issubset({r["tree_label"] for r in records}) or not leaves.issubset(headers):
        return None
    return tree


def recover_completed_attempt(outdir, signature, records):
    """Adopt a completed native run that an earlier helper could not recognize."""
    paths = [outdir / "run_manifest.json"] + sorted((outdir / "runs").glob("attempt_*.manifest.json"), reverse=True)
    checked = set()
    for path in paths:
        manifest = compare.read_json(path)
        command = manifest.get("command", [])
        if manifest.get("signature") != signature or "-o" not in command:
            continue
        try:
            output = (outdir / command[command.index("-o") + 1]).resolve()
        except (IndexError, TypeError):
            continue
        if output in checked or output.parent != (outdir / "runs").resolve():
            continue
        checked.add(output)
        tree = completed_attempt(output, records)
        if tree is None or (manifest.get("tree_sha256") and compare.sha256(tree) != manifest["tree_sha256"]):
            continue
        excluded = report_tree_membership(tree, records, outdir)
        manifest.update(recovered_from_status=manifest.get("status"), recovered_existing_attempt=True,
                        status="complete", tree=str(tree.relative_to(outdir)), tree_sha256=compare.sha256(tree),
                        excluded_genomes=excluded)
        compare.write_json(outdir / "run_manifest.json", manifest)
        compare.write_json(output.parent / (output.name + ".manifest.json"), manifest)
        return manifest
    return None


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--analysis-dir", type=Path, required=True, help="Completed spraynpray_compare output folder (cohort_ani)")
    p.add_argument("--outdir", type=Path, required=True, help="Separate phylogeny output folder; use a path without whitespace")
    p.add_argument("--references-per-sample", type=int, default=5, help="Top accepted ANI references per sample; 0 includes all accepted sample references")
    p.add_argument("--genomes", type=Path, help="Replacement sample FASTA folder if original paths moved")
    p.add_argument("--reference-dir", type=Path, help="Replacement exact-accession reference FASTA folder")
    p.add_argument("--assembly-info", type=Path, help="Optional NCBI metadata TSV to enrich reference labels")
    p.add_argument("--hmm-set", default="Bacteria", help="GToTree marker-set name or custom HMM path; e.g. Gammaproteobacteria for Pseudomonas")
    p.add_argument("--threads", type=int, default=8, help="GToTree parallel jobs; HMM and MUSCLE each use one thread per job")
    p.add_argument("--gtotree", default="GToTree", help="GToTree executable in the active environment")
    p.add_argument("--tree-program", choices=("FastTree", "FastTreeMP", "IQ-TREE"), help="Omit to use GToTree's default tree engine")
    p.add_argument("--nucleotide", action="store_true", help="GToTree -z: build a nucleotide rather than amino-acid marker alignment/tree")
    p.add_argument("--prepare-only", action="store_true", help="Stage genomes and write the GToTree command without running it")
    p.add_argument("--reuse-only", action="store_true", help="Reuse/recover a compatible completed tree; never launch GToTree")
    p.add_argument("--force", action="store_true", help="Run a new GToTree attempt even if compatible completed outputs exist")
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    for key in ("analysis_dir", "outdir", "genomes", "reference_dir", "assembly_info"):
        value = getattr(args, key)
        if value is not None:
            setattr(args, key, value.resolve())
    if any(c.isspace() for c in args.hmm_set):
        p.error("Use an HMM set name or path without whitespace for GToTree")
    if args.reuse_only and (args.force or args.prepare_only):
        p.error("--reuse-only cannot be combined with --force or --prepare-only")
    if args.threads < 1 or args.references_per_sample < 0:
        p.error("Threads must be positive; reference limit must be nonnegative")
    if not args.analysis_dir.is_dir() or args.analysis_dir == args.outdir or args.outdir in args.analysis_dir.parents:
        p.error("Provide an existing comparison folder and a separate phylogeny output folder")
    if any(c.isspace() for c in str(args.outdir)):
        p.error("GToTree requires this helper's output directory to have no whitespace")
    for folder in (args.genomes, args.reference_dir):
        if folder is not None and not folder.is_dir():
            p.error("Genome folders must exist")
    selection_path = args.outdir / "selection_manifest.json"
    if (args.outdir.exists() and any(args.outdir.iterdir())
            and compare.read_json(selection_path).get("tool") != "spraynpray_phylogeny"):
        p.error("Output directory is not a recognized phylogeny-helper run; choose a new folder")
    args.outdir.mkdir(parents=True, exist_ok=True)
    compare.write_json(selection_path, dict(tool="spraynpray_phylogeny", status="preparing"))
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    manifest_path = args.outdir / "run_manifest.json"
    try:
        records = prepare_inputs(args)
        hmm_path = Path(args.hmm_set).expanduser()
        if hmm_path.is_file():
            args.hmm_set = str(hmm_path.resolve())
        signature = dict(genomes=[dict(id=r["id"], label=r["tree_label"], sha256=r["sha256"]) for r in records],
                         hmm_set=args.hmm_set, hmm_sha256=compare.sha256(hmm_path) if hmm_path.is_file() else None,
                         nucleotide=args.nucleotide, tree_program=args.tree_program, executable=args.gtotree)
        previous = compare.read_json(manifest_path)
        tree = args.outdir / previous.get("tree", "missing")
        if (not args.prepare_only and not args.force and previous.get("status") == "complete"
                and previous.get("signature") == signature and tree.is_file()
                and compare.sha256(tree) == previous.get("tree_sha256")):
            report_tree_membership(tree, records, args.outdir)
            compare.write_json(selection_path, dict(tool="spraynpray_phylogeny", status="reused",
                               signature=signature, command=previous.get("command")))
            LOG.info("Reusing completed tree: %s", tree)
            return 0
        if not args.prepare_only and not args.force:
            recovered = recover_completed_attempt(args.outdir, signature, records)
            if recovered:
                compare.write_json(selection_path, dict(tool="spraynpray_phylogeny", status="recovered",
                                   signature=signature, command=recovered["command"]))
                LOG.info("Recovered completed GToTree output without rerunning: %s", args.outdir / recovered["tree"])
                return 0
        if args.reuse_only:
            raise ValueError("No compatible completed tree could be verified; --reuse-only prevents a new GToTree run")
        runs = args.outdir / "runs"
        runs.mkdir(exist_ok=True)
        attempt = next(i for i in range(1, 1000000) if not (runs / ("attempt_%03d" % i)).exists())
        output = "runs/attempt_%03d" % attempt
        command = [args.gtotree, "-f", "genomes.txt", "-m", "labels.tsv", "-H", args.hmm_set,
                   "-j", str(args.threads), "-n", "1", "-M", "1", "-o", output]
        if args.tree_program:
            command.extend(["-T", args.tree_program])
        if args.nucleotide:
            command.append("-z")
        script = "#!/usr/bin/env bash\nset -euo pipefail\ncd " + shlex.quote(str(args.outdir)) + "\nexec " + shlex.join(command) + "\n"
        (args.outdir / "run_gtotree.sh").write_text(script)
        compare.write_json(selection_path, dict(tool="spraynpray_phylogeny", signature=signature, command=command,
                           samples=sum(r["kind"] == "sample" for r in records), references=sum(r["kind"] == "reference" for r in records)))
        LOG.info("Selected %s input samples and %s unique references", sum(r["kind"] == "sample" for r in records),
                 sum(r["kind"] == "reference" for r in records))
        if args.prepare_only:
            LOG.info("Prepared inputs and command: %s", args.outdir / "run_gtotree.sh")
            return 0
        if shutil.which(args.gtotree) is None:
            raise ValueError("GToTree not found; activate its conda environment or use --prepare-only")
        log = args.outdir / ("gtotree_attempt_%03d.log" % attempt)
        manifest = dict(status="running", signature=signature, command=command, log=log.name)
        compare.write_json(manifest_path, manifest)
        try:
            compare.run_command(command, log, cwd=args.outdir)
            tree = completed_attempt(args.outdir / output, records)
            if tree is None:
                raise ValueError("GToTree did not produce a verifiable completed tree/alignment; inspect " + str(log))
            excluded = report_tree_membership(tree, records, args.outdir)
            manifest.update(status="complete", tree=str(tree.relative_to(args.outdir)), tree_sha256=compare.sha256(tree),
                            excluded_genomes=excluded)
            LOG.info("Tree: %s", tree)
        except (OSError, ValueError, subprocess.CalledProcessError):
            manifest["status"] = "failed"
            raise
        finally:
            (args.outdir / output).mkdir(parents=True, exist_ok=True)
            compare.write_json(manifest_path, manifest)
            compare.write_json(runs / ("attempt_%03d.manifest.json" % attempt), manifest)
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        LOG.error("%s", exc)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
