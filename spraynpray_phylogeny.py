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

VERSION = "0.5.0"
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



def normalize_tree_label(label):
    """Match GToTree's final `tr -s '_'` normalization of safe input labels."""
    return re.sub(r"_+", "_", label)


def compatible_signatures(saved, current):
    """Only repeated-underscore label changes are equivalent; sequence/settings stay exact."""
    if not isinstance(saved, dict) or not isinstance(current, dict):
        return False
    def normalized(signature):
        return dict(signature, genomes=[dict(genome, label=normalize_tree_label(genome.get("label", "")))
                                        for genome in signature.get("genomes", [])])
    return normalized(saved) == normalized(current)


def label_for(node):
    parts = [node["kind"], node["id"]]
    if node["kind"] == "reference":
        parts.extend(node.get(k, "") for k in ("organism", "strain", "isolate"))
    return normalize_tree_label(re.sub(r"[^A-Za-z0-9_.-]+", "_", "__".join(p for p in parts if p)).strip("_"))


def load_selection(args):
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
    return selected, rankings


def prepare_inputs(args):
    selected, rankings = load_selection(args)
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
    actual_labels = {normalize_tree_label(label): label for label in leaves}
    requested = [normalize_tree_label(r["tree_label"]) for r in records]
    if len(set(requested)) != len(requested):
        raise ValueError("Submitted labels collide after GToTree underscore normalization")
    if not set(actual_labels).intersection(requested):
        raise ValueError("Tree labels do not match the submitted genome labels")
    rows = [dict(id=r["id"], kind=r["kind"], tree_label=actual_labels.get(label, r["tree_label"]),
                 present_in_tree=label in actual_labels) for r, label in zip(records, requested)]
    compare.write_csv(outdir / "tree_membership.csv", ["id", "kind", "tree_label", "present_in_tree"], rows)
    missing = [r["id"] for r in rows if not r["present_in_tree"]]
    if missing:
        LOG.warning("GToTree excluded %s selected genomes; inspect tree_membership.csv and its marker/filtering reports", len(missing))
    return missing


def find_tree(output, explicit=None):
    """GToTree renames the final tree to the output directory's basename."""
    paths = [explicit] if explicit is not None else [output / name for name in
             (output.name + ".tre", "Aligned_SCGs_mod_names.tre", "Aligned_SCGs.tre")]
    for path in paths:
        if path.is_file():
            text = path.read_text().strip()
            if "(" in text and text.endswith(";"):
                return path
    return None


def completed_attempt(output, records, reasons=None, tree_path=None):
    """Require native completion, an alignment, and compatible final leaf labels."""
    def rejected(reason):
        if reasons is not None:
            reasons.append(reason)
        return None
    tree = find_tree(output, tree_path)
    log = output / "gtotree-runlog.txt"
    if tree is None:
        return rejected("no nonempty final Newick tree found")
    if not log.is_file():
        return rejected("native gtotree-runlog.txt is missing")
    text = log.read_text()
    if "Done!!" not in text or "Total process runtime:" not in text:
        return rejected("native log lacks the final completion/runtime markers")
    alignment = next((output / name for name in ("Aligned_SCGs_mod_names.faa", "Aligned_SCGs_mod_names.fa",
                     "Aligned_SCGs.faa", "Aligned_SCGs.fa")
                     if (output / name).is_file() and (output / name).stat().st_size), None)
    if alignment is None:
        return rejected("nonempty marker alignment is missing")
    leaves = set(re.findall(r"[,(]\s*([A-Za-z0-9_.-]+)\s*(?=[:),;])", tree.read_text()))
    with alignment.open() as handle:
        headers = {line[1:].strip().split()[0] for line in handle if line.startswith(">") and line[1:].strip()}
    expected = {normalize_tree_label(r["tree_label"]) for r in records}
    normalized_leaves = {normalize_tree_label(label) for label in leaves}
    normalized_headers = {normalize_tree_label(label) for label in headers}
    if len(expected) != len(records):
        return rejected("submitted labels collide after GToTree normalization")
    if not normalized_leaves:
        return rejected("no leaf labels could be read from the tree")
    if not normalized_leaves.issubset(expected):
        return rejected("unrecognized tree labels: " + ", ".join(sorted(normalized_leaves - expected)[:3]))
    if not normalized_leaves.issubset(normalized_headers):
        return rejected("tree labels absent from alignment: " + ", ".join(sorted(normalized_leaves - normalized_headers)[:3]))
    return tree


def recover_completed_attempt(outdir, signature, records, wanted_tree=None):
    """Adopt a completed native run that an earlier helper could not recognize."""
    paths = [outdir / "run_manifest.json"] + sorted((outdir / "runs").glob("attempt_*.manifest.json"), reverse=True)
    checked = set()
    for path in paths:
        manifest = compare.read_json(path)
        command = manifest.get("command", [])
        if not path.is_file():
            continue
        if not compatible_signatures(manifest.get("signature"), signature):
            fields = sorted(k for k in set(manifest.get("signature", {}) or {}) | set(signature)
                            if (manifest.get("signature", {}) or {}).get(k) != signature.get(k))
            LOG.warning("Cannot reuse %s: signature differs in %s", path, ", ".join(fields) or "missing signature")
            continue
        if "-o" not in command:
            LOG.warning("Cannot reuse %s: no recorded GToTree output argument", path)
            continue
        try:
            output = (outdir / command[command.index("-o") + 1]).resolve()
        except (IndexError, TypeError):
            continue
        if output in checked:
            continue
        if output.parent != (outdir / "runs").resolve():
            LOG.warning("Cannot reuse %s: recorded output %s is outside %s", path, output, outdir / "runs")
            continue
        checked.add(output)
        reasons = []
        tree = completed_attempt(output, records, reasons)
        if tree is None:
            LOG.warning("Cannot reuse %s: %s", output, "; ".join(reasons))
            continue
        if wanted_tree is not None and tree.resolve() != wanted_tree:
            continue
        if manifest.get("tree_sha256") and compare.sha256(tree) != manifest["tree_sha256"]:
            LOG.warning("Cannot reuse %s: tree checksum differs from the completed manifest", output)
            continue
        excluded = report_tree_membership(tree, records, outdir)
        manifest.update(recovered_from_status=manifest.get("status"), recovered_existing_attempt=True,
                        status="complete", tree=str(tree.relative_to(outdir)), tree_sha256=compare.sha256(tree),
                        excluded_genomes=excluded)
        compare.write_json(outdir / "run_manifest.json", manifest)
        compare.write_json(output.parent / (output.name + ".manifest.json"), manifest)
        return manifest
    return None



def recovery_records(args):
    """Inspect saved selection without reopening genomes or rewriting any inputs."""
    selected, _ = load_selection(args)
    path = args.outdir / "selected_genomes.csv"
    saved = {}
    if path.is_file():
        _, rows = compare.read_csv(path)
        saved = {r["id"]: r for r in rows}
        if len(saved) != len(rows) or set(saved) != {n["id"] for n in selected}:
            raise ValueError("Saved genome selection differs from the requested cohort/reference limit")
    records = []
    for node in selected:
        row = saved.get(node["id"], {})
        if row.get("kind", node["kind"]) != node["kind"]:
            raise ValueError("Saved sample/reference kind differs for " + node["id"])
        if node.get("ani_sha256") and row.get("sha256") and node["ani_sha256"] != row["sha256"]:
            raise ValueError("Saved genome fingerprint differs from ANI for " + node["id"])
        records.append(dict(node, tree_label=row.get("tree_label") or label_for(node),
                            sha256=row.get("sha256") or node.get("ani_sha256") or ""))
    return records


def run_signature(args, records):
    hmm = Path(args.hmm_set).expanduser()
    return dict(genomes=[dict(id=r["id"], label=r["tree_label"], sha256=r["sha256"]) for r in records],
                hmm_set=str(hmm.resolve()) if hmm.is_file() else args.hmm_set,
                hmm_sha256=compare.sha256(hmm) if hmm.is_file() else None,
                nucleotide=args.nucleotide, tree_program=args.tree_program, executable=args.gtotree)


def native_attempts(outdir, explicit_tree=None):
    if explicit_tree is not None:
        return [explicit_tree.parent]
    paths = sorted((outdir / "runs").glob("attempt_*"), reverse=True)
    return [p for p in paths if p.is_dir()] + ([outdir] if find_tree(outdir) else [])


def recover_native_outputs(args, signature, records):
    """Inspect native outputs even when no helper attempt manifest survived.

    Without an attempt-bound signature we report the existing artifact, never
    manufacture a sequence/settings-verified cache entry for it.
    """
    paths = [args.outdir / "run_manifest.json"] + list((args.outdir / "runs").glob("attempt_*.manifest.json"))
    if any(compare.read_json(path).get("signature") for path in paths):
        return None  # an explicitly incompatible recorded run must not be bypassed
    attempts = native_attempts(args.outdir, args.tree)
    LOG.info("No attempt-bound fingerprints found; inspecting %s native output directories under %s", len(attempts), args.outdir)
    for output in attempts:
        reasons = []
        tree = completed_attempt(output, records, reasons, args.tree)
        if tree is None:
            LOG.warning("Cannot recover native output %s: %s", output, "; ".join(reasons))
            continue
        if args.tree is not None and tree.resolve() != args.tree:
            LOG.warning("Native final tree %s differs from explicitly requested %s", tree, args.tree)
            continue
        log = (output / "gtotree-runlog.txt").read_text()
        counts = re.search(r"Overall,\s*(\d+) genomes of the input (\d+) were retained", log)
        all_retained = re.search(r"All (\d+) input genomes were retained", log)
        retained, submitted = (tuple(map(int, counts.groups())) if counts else
                              (int(all_retained[1]), int(all_retained[1])) if all_retained else (None, None))
        leaves = set(re.findall(r"[,(]\s*([A-Za-z0-9_.-]+)\s*(?=[:),;])", tree.read_text()))
        if submitted != len(records) or retained != len(leaves) or retained < 3:
            LOG.warning("Cannot recover native output %s: native submitted/retained counts %s/%s do not match requested %s genomes and %s tree leaves",
                        output, submitted, retained, len(records), len(leaves))
            continue
        excluded = report_tree_membership(tree, records, args.outdir)
        manifest = dict(status="recovered_native", tree=str(tree.resolve()), tree_sha256=compare.sha256(tree),
                        verification="native_completion_alignment_and_genome_labels",
                        input_fingerprints_verified=False, requested_signature=signature, excluded_genomes=excluded)
        compare.write_json(args.outdir / "recovered_native_manifest.json", manifest)
        LOG.warning("Native tree recovered by completion/alignment/label checks; original genome fingerprints and tree settings cannot be verified without an attempt manifest")
        LOG.info("Existing tree: %s", tree.resolve())
        return manifest
    return None


def explain_missing_outputs(args):
    attempts = native_attempts(args.outdir, args.tree)
    LOG.error("Searched %s native output directories and helper manifests under %s", len(attempts), args.outdir)
    if not attempts:
        roots = {args.outdir.parent, args.analysis_dir.parent}
        found = set()
        for root in roots:
            for tree in root.glob("*phylogen*/runs/attempt_*/*.tre"):
                if tree.is_file():
                    found.add(str(tree.resolve()))
            for tree in root.glob("runs/attempt_*/*.tre"):
                if tree.is_file():
                    found.add(str(tree.resolve()))
        for tree in sorted(found)[:10]:
            LOG.error("Tree found nearby: %s (use the containing phylogeny folder as --outdir)", tree)


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--version", action="version", version=VERSION)
    p.add_argument("--tree", type=Path, help="Explicit final .tre file for --reuse-only when recovering native outputs")
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
    for key in ("analysis_dir", "outdir", "genomes", "reference_dir", "assembly_info", "tree"):
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
    if args.tree is not None and not args.reuse_only:
        p.error("--tree is only used with --reuse-only")
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    LOG.info("SprayNPray phylogeny helper %s; output folder: %s", VERSION, args.outdir)
    selection_path = args.outdir / "selection_manifest.json"
    manifest_path = args.outdir / "run_manifest.json"
    try:
        # Recovery precedes ownership markers, staging, and executable checks.
        # Failed recovery leaves saved inputs and manifests unchanged.
        if not args.force and not args.prepare_only:
            records = recovery_records(args)
            signature = run_signature(args, records)
            previous = compare.read_json(manifest_path)
            tree = args.outdir / previous.get("tree", "missing")
            if (previous.get("status") == "complete" and compatible_signatures(previous.get("signature"), signature)
                    and tree.is_file() and compare.sha256(tree) == previous.get("tree_sha256")):
                if args.tree is None or tree.resolve() == args.tree:
                    report_tree_membership(tree, records, args.outdir)
                    LOG.info("Reusing completed tree: %s", tree)
                    return 0
            recovered = recover_completed_attempt(args.outdir, signature, records, args.tree)
            if recovered:
                LOG.info("Recovered completed GToTree output without rerunning: %s", args.outdir / recovered["tree"])
                return 0
            if recover_native_outputs(args, signature, records):
                return 0
            if args.reuse_only:
                explain_missing_outputs(args)
                raise ValueError("No completed tree could be verified; --reuse-only ran no tools and changed no input files")
            if any(find_tree(output) for output in native_attempts(args.outdir)):
                raise ValueError("Existing trees were found but recovery checks rejected them. No new run started; inspect the warnings or use --force deliberately")
        if not args.prepare_only and shutil.which(args.gtotree) is None:
            raise ValueError("No completed tree was recovered and GToTree is not installed in this environment. No inputs were changed; activate its environment only to start a new tree")
        if (args.outdir.exists() and any(args.outdir.iterdir())
                and compare.read_json(selection_path).get("tool") != "spraynpray_phylogeny"):
            raise ValueError("Output directory is not a recognized phylogeny-helper run; choose a new folder")
        args.outdir.mkdir(parents=True, exist_ok=True)
        if not selection_path.is_file():
            compare.write_json(selection_path, dict(tool="spraynpray_phylogeny", status="preparing"))
        records = prepare_inputs(args)
        signature = run_signature(args, records)
        if Path(args.hmm_set).expanduser().is_file():
            args.hmm_set = str(Path(args.hmm_set).expanduser().resolve())
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
        log = args.outdir / ("gtotree_attempt_%03d.log" % attempt)
        manifest = dict(status="running", signature=signature, command=command, log=log.name)
        compare.write_json(manifest_path, manifest)
        try:
            compare.run_command(command, log, cwd=args.outdir)
            reasons = []
            tree = completed_attempt(args.outdir / output, records, reasons)
            if tree is None:
                raise ValueError("GToTree output could not be verified: %s; inspect %s" % ("; ".join(reasons), log))
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
