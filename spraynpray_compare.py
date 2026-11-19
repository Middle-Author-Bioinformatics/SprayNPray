#!/usr/bin/env python3
"""Summarize SprayNPray hits, calculate all-vs-all FastANI, and map genomes.

Standalone Python 3.9+ script. Plotting requires numpy, scipy and matplotlib.
FastANI and (for missing reference downloads) NCBI datasets are external tools.
See README_compare.md for methods, interpretation and example commands.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import itertools
import json
import logging
import math
import re
import shutil
import statistics
import subprocess
import sys
import textwrap
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

VERSION = "0.3.0"
LOG = logging.getLogger("spraynpray_compare")
ACCESSION = re.compile(r"GC[AF]_\d+\.\d+")
ENTRY = re.compile(
    r"(?P<orf>\S+): (?:(?P<nohit>no hit)|(?P<organism>.*?) "
    r"\[(?P<protein>[^;\[\]]+); (?P<assembly>GC[AF]_\d+\.\d+|unmapped); "
    r"(?P<identity>\d+(?:\.\d+)?)% identity\])")
FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(value, indent=2) + "\n")
    temp.replace(path)


def read_json(path):
    try:
        value = json.loads(path.read_text())
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def write_csv(path, fields, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            break
        except OverflowError:
            limit //= 10
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("Missing CSV header: " + str(path))
        rows = list(reader)
        if any(None in r or any(v is None for v in r.values()) for r in rows):
            raise ValueError("Malformed CSV row: " + str(path))
        return reader.fieldnames, rows


def identity(value):
    number = float(value)
    if not math.isfinite(number) or not 0 <= number <= 100:
        raise ValueError("Identity must be finite and between 0 and 100: " + str(value))
    return number


def parse_summary_entries(text):
    position = 0
    while position < len(text):
        match = ENTRY.match(text, position)
        if match is None:
            raise ValueError("Unsupported hit summary near: " + text[position:position + 120])
        row = match.groupdict()
        yield dict(orf=row["orf"], assembly="" if row["nohit"] else row["assembly"],
                   organism=row["organism"] or "", protein=row["protein"] or "",
                   identity=None if row["nohit"] else identity(row["identity"]))
        position = match.end()
        if position == len(text):
            break
        if text[position:position + 2] != "; ":
            raise ValueError("Unexpected text after ORF match: " + text[position:position + 120])
        position += 2
        if position == len(text):
            raise ValueError("Trailing separator in ORF match summary")


ANNOTATION_FIELDS = ["strain", "isolate", "infraspecific_name", "taxid", "bioproject", "biosample",
                     "assembly_level", "metadata_assembly_accession", "metadata_match"]


def annotation_values(row, prefix=""):
    def clean(field):
        value = row.get(prefix + field, "").strip()
        return "" if value.lower() in ("na", "n/a", "none", "not available", "-") else value
    infra = clean("infraspecific_name")
    return dict(strain=infra.split("=", 1)[1] if infra.lower().startswith("strain=") else "",
                isolate=clean("isolate"), infraspecific_name=infra, taxid=clean("taxid"),
                bioproject=clean("bioproject"), biosample=clean("biosample"), assembly_level=clean("assembly_level"),
                metadata_assembly_accession=clean("assembly_accession"))


def enrich_reference_nodes(nodes, assembly_info=None):
    """Stream optional NCBI metadata, matching exact versioned accessions only."""
    if assembly_info is None:
        return
    wanted = {n["id"]: n for n in nodes if n["kind"] == "reference"}
    opener = gzip.open if str(assembly_info).endswith(".gz") else open
    header, seen = None, {}
    with opener(assembly_info, "rt") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("##"):
                continue
            if line.lstrip("# ").startswith("assembly_accession\t"):
                header = line.lstrip("# ").rstrip("\r\n").split("\t")
                if len(header) != len(set(header)) or "organism_name" not in header:
                    raise ValueError("Invalid assembly metadata header")
                continue
            if line.startswith("#"):
                continue
            if header is None:
                raise ValueError("Missing TAB-separated assembly_accession header in metadata")
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(header):
                raise ValueError("Assembly metadata row width mismatch at line %s" % number)
            row = dict(zip(header, fields))
            accession = row["assembly_accession"]
            if accession in wanted:
                if accession in seen and seen[accession] != row:
                    raise ValueError("Conflicting assembly metadata for " + accession)
                seen[accession] = row
                node = wanted[accession]
                node.update({k: v for k, v in annotation_values(row).items() if v})
                if row.get("organism_name", "").strip().lower() not in ("", "na"):
                    node["organism"] = row["organism_name"]
                node["metadata_match"] = "exact"
    if header is None:
        raise ValueError("Missing assembly metadata header")


def load_sample(path):
    fields, rows = read_csv(path)
    hits, selected = [], 0
    if "closest_blast_hits" in fields:
        seen_contigs = set()
        for row in rows:
            if row["contig"] in seen_contigs:
                raise ValueError("Duplicate contig in " + str(path))
            seen_contigs.add(row["contig"])
            if int(row.get("summary_orfs_omitted", "0") or 0):
                raise ValueError("Truncated summary in %s. Use the run's complete top_hits.csv, renamed to the sample ID." % path)
            entries = list(parse_summary_entries(row["closest_blast_hits"]))
            expected = int(row.get("selected_orfs", len(entries)))
            if len(entries) != expected:
                raise ValueError("Selected ORF count does not match listed entries in %s, contig %s" % (path, row["contig"]))
            selected += expected
            hits.extend(dict(entry, contig=row["contig"]) for entry in entries)
    elif {"orf", "protein_accession", "source_genome_accession", "pident"}.issubset(fields):
        selected = len(rows)
        for row in rows:
            no_hit = row.get("status") == "no_hit"
            assembly = row["source_genome_accession"] or "unmapped"
            if not no_hit and assembly != "unmapped" and ACCESSION.fullmatch(assembly) is None:
                raise ValueError("Invalid source assembly in " + str(path))
            hits.append(dict(orf=row["orf"], contig=row.get("contig", ""),
                             assembly="" if no_hit else assembly, organism=row.get("taxa", ""),
                             protein=row["protein_accession"], identity=None if no_hit else identity(row["pident"]),
                             **dict(annotation_values(row, "assembly_"), metadata_match=row.get("metadata_match", ""))))
    else:
        raise ValueError("Expected V2 contig summary or top_hits.csv columns: " + str(path))
    if len({h["orf"] for h in hits}) != len(hits):
        raise ValueError("Duplicate ORF IDs in " + str(path))
    matched = [h for h in hits if h["identity"] is not None]
    grouped = defaultdict(list)
    for hit in matched:
        grouped[hit["assembly"]].append(hit)
    summary = []
    for assembly, group in sorted(grouped.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        values = [h["identity"] for h in group]
        summary.append(dict(sample=path.stem, reference_assembly=assembly,
                            organism_names=" | ".join(sorted({h["organism"] for h in group})),
                            orf_hits=len(group), matched_orfs=len(matched), selected_orfs=selected,
                            hit_percent=100 * len(group) / len(matched),
                            selected_orf_percent=100 * len(group) / selected if selected else 0,
                            mean_protein_identity=statistics.mean(values),
                            min_protein_identity=min(values), max_protein_identity=max(values),
                            **{field: " | ".join(sorted({h.get(field, "") for h in group if h.get(field)}))
                               for field in ANNOTATION_FIELDS}))
    overview = dict(sample=path.stem, source_csv=str(path.resolve()), selected_orfs=selected,
                    matched_orfs=len(matched), no_hit_orfs=selected - len(matched),
                    unmapped_hit_orfs=sum(h["assembly"] == "unmapped" for h in matched),
                    reference_assemblies=sum(a != "unmapped" for a in grouped),
                    mean_protein_identity=statistics.mean(h["identity"] for h in matched) if matched else "")
    return summary, overview, hits


def summarize(folder, outdir, top_references=0):
    files = sorted(folder.glob("*.csv"))
    if not files:
        raise ValueError("No .csv files found in " + str(folder))
    master, overviews, hits = [], [], []
    selected_refs = set()
    for path in files:
        summary, overview, sample_hits = load_sample(path)
        master.extend(summary)
        overviews.append(overview)
        hits.extend(dict(h, sample=path.stem) for h in sample_hits)
        references = [r for r in summary if r["reference_assembly"] != "unmapped"]
        if top_references:
            references = references[:top_references]
        selected_refs.update(r["reference_assembly"] for r in references)
    fields = ["sample", "reference_assembly", "organism_names", "orf_hits", "matched_orfs", "selected_orfs",
              "hit_percent", "selected_orf_percent", "mean_protein_identity", "min_protein_identity", "max_protein_identity"] + ANNOTATION_FIELDS
    write_csv(outdir / "master_summary.csv", fields, master)
    write_csv(outdir / "sample_overview.csv", list(overviews[0]), overviews)
    write_csv(outdir / "orf_hits.csv", ["sample", "contig", "orf", "assembly", "organism", "protein", "identity"] + ANNOTATION_FIELDS, hits)
    references = sorted({r["reference_assembly"] for r in master if r["reference_assembly"] != "unmapped"})
    write_csv(outdir / "references.csv", ["assembly", "organism_names", "total_orf_hits", "samples_with_hits", "selected_for_ani"], [
        dict(assembly=a, organism_names=" | ".join(sorted({r["organism_names"] for r in master if r["reference_assembly"] == a})),
             total_orf_hits=sum(r["orf_hits"] for r in master if r["reference_assembly"] == a),
             samples_with_hits=sum(r["reference_assembly"] == a for r in master), selected_for_ani=a in selected_refs)
        for a in references])
    return master, overviews, sorted(selected_refs)


def fasta_records(path):
    header, parts = None, []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    if not parts:
                        raise ValueError("Empty FASTA record in " + str(path))
                    yield header, "".join(parts)
                header, parts = line[1:], []
                if not header:
                    raise ValueError("Empty FASTA header in " + str(path))
            else:
                if header is None:
                    raise ValueError("Sequence before FASTA header in " + str(path))
                sequence = "".join(line.split()).upper()
                if re.fullmatch("[ACGTRYSWKMBDHVN]+", sequence) is None:
                    raise ValueError("Non-DNA characters in " + str(path))
                parts.append(sequence)
    if header is not None:
        if not parts:
            raise ValueError("Empty FASTA record in " + str(path))
        yield header, "".join(parts)


def stage_fasta(source, target):
    lengths = []
    temp = target.with_name(target.name + ".tmp")
    with open(temp, "w") as handle:
        for header, sequence in fasta_records(source):
            lengths.append(len(sequence))
            handle.write(">%s\n%s\n" % (header, sequence))
    if not lengths:
        raise ValueError("No sequences in " + str(source))
    digest = sha256(temp)
    if target.is_file() and sha256(target) == digest:
        temp.unlink()
    else:
        temp.replace(target)
    total, cumulative, n50 = sum(lengths), 0, 0
    for length in sorted(lengths, reverse=True):
        cumulative += length
        if cumulative >= total / 2:
            n50 = length
            break
    return dict(length_bp=total, contigs=len(lengths), n50=n50, sha256=digest)


def reference_files(folder):
    result = {}
    if folder is None or not folder.exists():
        return result
    for path in sorted(folder.rglob("*")):
        if not path.is_file() or not str(path).endswith(FASTA_SUFFIXES):
            continue
        match = ACCESSION.search(path.name)
        if match is None:
            match = ACCESSION.fullmatch(path.parent.name)
        if match:
            accession = match.group()
            if accession in result:
                raise ValueError("Multiple reference FASTAs for %s in %s" % (accession, folder))
            result[accession] = path.resolve()
    return result


def run_command(command, log_path, cwd=None):
    LOG.info("Running %s", " ".join(map(str, command)))
    with open(log_path, "w") as handle:
        subprocess.run(list(map(str, command)), stdout=handle, stderr=subprocess.STDOUT, cwd=cwd, check=True)


def extract_reference_package(archive_path, destination, requested):
    """Extract exact requested versions only; never extract untrusted zip paths."""
    extracted = {}
    with zipfile.ZipFile(archive_path) as archive:
        for member in archive.infolist():
            parts = Path(member.filename).parts
            if len(parts) < 2 or not parts[-1].endswith("_genomic.fna"):
                continue
            accession = parts[-2]
            if accession not in requested or "_cds_from_" in parts[-1]:
                continue
            if accession in extracted:
                raise ValueError("Multiple genomic FASTAs in download for " + accession)
            target = destination / (accession + ".fna")
            temporary = target.with_name(target.name + ".partial")
            with archive.open(member) as source, open(temporary, "wb") as out:
                shutil.copyfileobj(source, out)
            # Check DNA and complete records before accepting the reference.
            if not any(True for _ in fasta_records(temporary)):
                raise ValueError("Empty downloaded reference: " + accession)
            temporary.replace(target)
            extracted[accession] = target.resolve()
    return extracted


def resolve_references(accessions, args, work):
    destination = work / "references"
    destination.mkdir(exist_ok=True)
    local = reference_files(args.reference_dir)
    cached = reference_files(destination)
    resolved = {a: local.get(a, cached.get(a)) for a in accessions}
    missing = [a for a in accessions if resolved[a] is None]
    if missing and not args.no_download:
        if shutil.which(args.datasets) is None:
            raise ValueError("Missing reference FASTAs. Install NCBI datasets or supply --reference-dir. Missing: " + ", ".join(missing[:10]))
        for start in range(0, len(missing), 50):
            batch = missing[start:start + 50]
            list_path = work / "reference_download_list.txt"
            list_path.write_text("\n".join(batch) + "\n")
            archive = work / "reference_download.zip"
            archive.unlink(missing_ok=True)
            try:
                run_command([args.datasets, "download", "genome", "accession", "--inputfile", list_path,
                             "--include", "genome", "--filename", archive, "--no-progressbar"],
                            work / ("reference_download_%s.log" % start))
                resolved.update(extract_reference_package(archive, destination, set(batch)))
            except (subprocess.CalledProcessError, zipfile.BadZipFile) as exc:
                if not args.allow_missing_references:
                    raise
                LOG.warning("Reference download batch failed; its missing references will lack ANI: %s", exc)
        missing = [a for a in accessions if resolved.get(a) is None]
    if missing:
        write_csv(args.outdir / "missing_references.csv", ["assembly"], [dict(assembly=a) for a in missing])
        if not args.allow_missing_references:
            raise ValueError("Missing exact reference versions; see missing_references.csv. Supply local FASTAs or use --allow-missing-references to retain them without ANI.")
        LOG.warning("%s references have no sequence; retained in tables but excluded from the default ANI-ranked map", len(missing))
    elif (args.outdir / "missing_references.csv").exists():
        (args.outdir / "missing_references.csv").unlink()
    return resolved


def prepare_nodes(overviews, references, master, args, work):
    nodes = []
    for overview in overviews:
        sample = overview["sample"]
        path = args.genomes / (sample + args.genome_extension)
        if not path.is_file():
            raise ValueError("Sample FASTA not found: " + str(path))
        nodes.append(dict(id=sample, kind="sample", organism="", source_file=str(path.resolve()),
                          matched_orfs=overview["matched_orfs"], selected_orfs=overview["selected_orfs"]))
    reference_paths = resolve_references(references, args, work)
    for reference in references:
        groups = [r for r in master if r["reference_assembly"] == reference]
        nodes.append(dict(id=reference, kind="reference", organism=" | ".join(sorted({r["organism_names"] for r in groups})),
                          source_file=str(reference_paths[reference]) if reference_paths[reference] else "",
                          orf_hits=sum(r["orf_hits"] for r in groups),
                          **{field: " | ".join(sorted({r.get(field, "") for r in groups if r.get(field)}))
                             for field in ANNOTATION_FIELDS}))
    enrich_reference_nodes(nodes, getattr(args, "assembly_info", None))
    if len({n["id"] for n in nodes}) != len(nodes):
        raise ValueError("A sample filename collides with a reference accession; rename that sample")
    staged = work / "genomes"
    staged.mkdir(exist_ok=True)
    for i, node in enumerate(nodes):
        if node["source_file"]:
            relative = "genomes/g%06d.fa" % i
            node.update(stage_fasta(Path(node["source_file"]), work / relative))
            node.update(staged_file=relative, sequence_status="available")
        else:
            node.update(staged_file="", sequence_status="missing", length_bp="", contigs="", n50="")
    return nodes


def parse_fastani(path, nodes, work):
    aliases = {}
    for node in nodes:
        if node["staged_file"]:
            for alias in (node["id"], node["source_file"], node["staged_file"], str((work / node["staged_file"]).resolve())):
                if alias in aliases and aliases[alias] != node["id"]:
                    raise ValueError("Ambiguous genome path in ANI inputs: " + alias)
                aliases[alias] = node["id"]
    result = {}
    with open(path) as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 5:
                fields = raw.split()
            if len(fields) != 5:
                raise ValueError("FastANI line %s must have five fields" % number)
            query, reference, ani, mapped, fragments = fields
            if query not in aliases or reference not in aliases:
                raise ValueError("Unknown genome path in FastANI line %s" % number)
            mapped, fragments = int(mapped), int(fragments)
            if not 0 <= mapped <= fragments or fragments <= 0:
                raise ValueError("Invalid fragment counts in FastANI line %s" % number)
            value = dict(ani=identity(ani), mapped_fragments=mapped, query_fragments=fragments,
                         aligned_fraction=mapped / fragments)
            key = aliases[query], aliases[reference]
            if key in result and result[key] != value:
                raise ValueError("Conflicting duplicate FastANI comparison: " + str(key))
            result[key] = value
    return result


def calculate_ani(nodes, args, work):
    scope = getattr(args, "ani_scope", "all")
    base_signature = dict(genomes=[dict(id=n["id"], sha256=n["sha256"]) for n in nodes if n["staged_file"]],
                          fragment_length=args.frag_len, min_fraction=0.2, executable=args.fastani)
    signature = dict(base_signature, scope=scope)
    raw_path, cache_path = work / "fastani.tsv", work / "ani_cache.json"
    cache = read_json(cache_path)
    if args.fastani_results:
        return parse_fastani(args.fastani_results, nodes, work), "provided"
    compatible = [signature, dict(base_signature, scope="all"), base_signature] if scope == "samples" else [signature, base_signature]
    if (not args.force_ani and cache.get("signature") in compatible and raw_path.is_file()
            and cache.get("output_sha256") == sha256(raw_path)):
        LOG.info("Reusing completed ANI results (including compatible older all-vs-all runs)")
        return parse_fastani(raw_path, nodes, work), "reused"
    if shutil.which(args.fastani) is None:
        raise ValueError("FastANI executable not found: " + args.fastani)
    available = [n for n in nodes if n["staged_file"]]
    if scope == "samples":
        samples = [n for n in available if n["kind"] == "sample"]
        references = [n for n in available if n["kind"] == "reference"]
        requested = len(samples) * len(available) + len(references) * len(samples)
        LOG.info("Requesting %s directional comparisons instead of %s full-union comparisons", requested, len(available) ** 2)
        jobs = [("samples", samples, available)]
        if references:
            jobs.append(("references", references, samples))
        outputs = []
        for name, queries, targets in jobs:
            job_signature = dict(base_signature, queries=[n["id"] for n in queries], references=[n["id"] for n in targets])
            job_output = work / ("fastani_%s.tsv" % name)
            job_cache_path = work / ("fastani_%s_cache.json" % name)
            job_cache = read_json(job_cache_path)
            if (not args.force_ani and job_cache.get("signature") == job_signature and job_output.is_file()
                    and job_cache.get("output_sha256") == sha256(job_output)):
                LOG.info("Reusing completed ANI job: %s", name)
                parse_fastani(job_output, nodes, work)
            else:
                qlist, rlist = work / (name + "_queries.txt"), work / (name + "_references.txt")
                qlist.write_text("\n".join(n["staged_file"] for n in queries) + "\n")
                rlist.write_text("\n".join(n["staged_file"] for n in targets) + "\n")
                partial = work / (job_output.name + ".partial")
                partial.unlink(missing_ok=True)
                run_command([args.fastani, "--ql", qlist.name, "--rl", rlist.name, "-o", partial.name,
                             "-t", args.threads, "--fragLen", args.frag_len, "--minFraction", "0.2"],
                            work / ("fastani_%s.log" % name), cwd=work)
                parse_fastani(partial, nodes, work)
                partial.replace(job_output)
                write_json(job_cache_path, dict(signature=job_signature, output_sha256=sha256(job_output)))
            outputs.append(job_output)
        combined = work / "fastani.tsv.partial"
        with open(combined, "wb") as destination:
            for output in outputs:
                with open(output, "rb") as source:
                    shutil.copyfileobj(source, destination)
                destination.write(b"\n")
        values = parse_fastani(combined, nodes, work)
        combined.replace(raw_path)
        write_json(cache_path, dict(signature=signature, output_sha256=sha256(raw_path)))
        return values, "ran"
    list_path = work / "genomes.txt"
    list_path.write_text("\n".join(n["staged_file"] for n in nodes if n["staged_file"]) + "\n")
    partial = work / "fastani.tsv.partial"
    partial.unlink(missing_ok=True)
    run_command([args.fastani, "--ql", "genomes.txt", "--rl", "genomes.txt", "-o", partial.name,
                 "-t", args.threads, "--fragLen", args.frag_len, "--minFraction", "0.2"], work / "fastani.log", cwd=work)
    values = parse_fastani(partial, nodes, work)
    partial.replace(raw_path)
    write_json(cache_path, dict(signature=signature, output_sha256=sha256(raw_path)))
    return values, "ran"


PAIR_FIELDS = ["genome_a", "genome_b", "kind_a", "kind_b", "ani_a_to_b", "ani_b_to_a", "ani_mean",
               "aligned_fraction_a_to_b", "aligned_fraction_b_to_a", "min_aligned_fraction", "ani_used", "status"]


def pairwise_tables(nodes, directional, args):
    pairs, edges = [], []
    for a, b in itertools.combinations(nodes, 2):
        if getattr(args, "ani_scope", "all") == "samples" and a["kind"] == b["kind"] == "reference":
            continue
        forward, reverse = directional.get((a["id"], b["id"])), directional.get((b["id"], a["id"]))
        values = [v["ani"] for v in (forward, reverse) if v]
        mean = statistics.mean(values) if values else ""
        fraction = min(forward["aligned_fraction"], reverse["aligned_fraction"]) if forward and reverse else ""
        status = ("missing_sequence" if not a["staged_file"] or not b["staged_file"] else
                  "not_reported" if not values else "missing_direction" if not forward or not reverse else
                  "low_aligned_fraction" if fraction < args.min_af else
                  "below_min_ani" if mean < args.min_ani else "accepted")
        record = dict(genome_a=a["id"], genome_b=b["id"], kind_a=a["kind"], kind_b=b["kind"],
                      ani_a_to_b=forward["ani"] if forward else "", ani_b_to_a=reverse["ani"] if reverse else "",
                      ani_mean=mean, aligned_fraction_a_to_b=forward["aligned_fraction"] if forward else "",
                      aligned_fraction_b_to_a=reverse["aligned_fraction"] if reverse else "",
                      min_aligned_fraction=fraction, ani_used=mean if status == "accepted" else "", status=status)
        pairs.append(record)
        if status == "accepted":
            edges.append(dict(record, distance=max(100 - mean, 0.001)))
    write_csv(args.outdir / "ani_pairs.csv", PAIR_FIELDS, pairs)
    write_csv(args.outdir / "ani_directional.csv", ["query", "reference", "ani", "mapped_fragments", "query_fragments", "aligned_fraction"],
              [dict(query=q, reference=r, **value) for (q, r), value in sorted(directional.items())])
    lookup = {(e["genome_a"], e["genome_b"]): e for e in edges}
    matrix = []
    af_matrix = []
    for a in nodes:
        row, af_row = {"genome": a["id"]}, {"genome": a["id"]}
        for b in nodes:
            edge = lookup.get((a["id"], b["id"]), lookup.get((b["id"], a["id"])))
            row[b["id"]] = 100 if a == b and a["staged_file"] else edge["ani_used"] if edge else ""
            af_row[b["id"]] = 1 if a == b and a["staged_file"] else edge["min_aligned_fraction"] if edge else ""
        matrix.append(row)
        af_matrix.append(af_row)
    fields = ["genome"] + [n["id"] for n in nodes]
    write_csv(args.outdir / "ani_matrix.csv", fields, matrix)
    write_csv(args.outdir / "aligned_fraction_matrix.csv", fields, af_matrix)
    return edges


def nearest_reference_summary(nodes, pairs, outdir):
    """Rank accepted genome ANI, never approximate 2D display distances."""
    kinds = {n["id"]: n["kind"] for n in nodes}
    candidates = defaultdict(list)
    for pair in pairs:
        if pair["status"] != "accepted":
            continue
        a, b = pair["genome_a"], pair["genome_b"]
        if {kinds.get(a), kinds.get(b)} != {"sample", "reference"}:
            continue
        sample, reference = (a, b) if kinds[a] == "sample" else (b, a)
        candidates[sample].append((identity(pair["ani_used"]), float(pair["min_aligned_fraction"]), reference))
    result = []
    for node in nodes:
        if node["kind"] != "sample":
            continue
        matches = sorted(candidates[node["id"]], key=lambda item: (-item[0], -item[1], item[2]))
        if matches:
            ani, fraction, reference = matches[0]
            ties = sum(abs(value[0] - ani) <= 1e-6 for value in matches)
            result.append(dict(sample=node["id"], nearest_reference=reference, estimated_ani=round(ani, 5),
                               ani_distance=round(max(0.0, 100 - ani), 5), min_aligned_fraction=round(fraction, 6),
                               tied_reference_count=ties, status="matched"))
        else:
            result.append(dict(sample=node["id"], nearest_reference="", estimated_ani="", ani_distance="",
                               min_aligned_fraction="", tied_reference_count=0, status="no_accepted_reference"))
    write_csv(outdir / "nearest_references.csv", ["sample", "nearest_reference", "estimated_ani", "ani_distance",
              "min_aligned_fraction", "tied_reference_count", "status"], result)
    return result


def connected_components(ids, edges):
    adjacent = {n: set() for n in ids}
    for edge in edges:
        a, b = edge["genome_a"], edge["genome_b"]
        adjacent[a].add(b)
        adjacent[b].add(a)
    remaining, result = set(ids), []
    for first in ids:
        if first not in remaining:
            continue
        pending, component = [first], []
        remaining.remove(first)
        while pending:
            node = pending.pop()
            component.append(node)
            for other in sorted(adjacent[node]):
                if other in remaining:
                    pending.append(other)
                    remaining.remove(other)
        result.append(component)
    return result


def embed_network(nodes, edges):
    import numpy as np
    from scipy.optimize import least_squares
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import shortest_path

    positions, diagnostics, offset = {}, [], 0.0
    for number, component in enumerate(connected_components([n["id"] for n in nodes], edges), 1):
        ids = {name: i for i, name in enumerate(component)}
        constraints = [e for e in edges if e["genome_a"] in ids and e["genome_b"] in ids]
        n = len(ids)
        if n == 1:
            coords = np.zeros((1, 2))
            stress, error = 0.0, 0.0
        else:
            i = np.array([ids[e["genome_a"]] for e in constraints])
            j = np.array([ids[e["genome_b"]] for e in constraints])
            target = np.array([e["distance"] for e in constraints])
            distances = np.full((n, n), np.inf)
            np.fill_diagonal(distances, 0)
            distances[i, j] = distances[j, i] = target
            # Missing distances are used only for initialization, never as ANI observations.
            initial = shortest_path(distances, directed=False)
            centering = np.eye(n) - np.ones((n, n)) / n
            gram = -0.5 * centering @ (initial ** 2) @ centering
            values, vectors = np.linalg.eigh(gram)
            order = np.argsort(values)[::-1][:2]
            coords = np.zeros((n, 2))
            coords[:, :len(order)] = vectors[:, order] * np.sqrt(np.maximum(values[order], 0))
            rng = np.random.default_rng(42)
            coords += rng.normal(0, max(float(np.max(target)), 0.001) * 1e-6, coords.shape)

            def residual(flat):
                x = flat.reshape(n, 2)
                return np.linalg.norm(x[i] - x[j], axis=1) - target

            def jacobian(flat):
                x = flat.reshape(n, 2)
                delta = x[i] - x[j]
                derivative = delta / np.maximum(np.linalg.norm(delta, axis=1), 1e-12)[:, None]
                rows = np.repeat(np.arange(len(i)), 4)
                columns = np.column_stack((2 * i, 2 * i + 1, 2 * j, 2 * j + 1)).ravel()
                data = np.column_stack((derivative, -derivative)).ravel()
                return csr_matrix((data, (rows, columns)), shape=(len(i), 2 * n))

            fit = least_squares(residual, coords.ravel(), jac=jacobian, max_nfev=300,
                                ftol=1e-7, xtol=1e-7, gtol=1e-7)
            coords = fit.x.reshape(n, 2)
            difference = residual(fit.x)
            stress = math.sqrt(float(np.sum(difference ** 2) / np.sum(target ** 2)))
            error = float(np.max(np.abs(difference)))
            if not fit.success:
                LOG.warning("Component %s layout reached its optimization limit", number)
        coords -= coords.mean(axis=0)
        width = float(np.ptp(coords[:, 0]))
        coords[:, 0] += offset - float(coords[:, 0].min())
        offset += max(width, 1.0) + 2.0
        for name, position in zip(component, coords):
            positions[name] = dict(x=float(position[0]), y=float(position[1]), component=number)
        diagnostics.append(dict(component=number, nodes=n, measured_pairs=len(constraints),
                                normalized_stress=stress, max_distance_error=error))
    return positions, diagnostics


def display_edges(nodes, edges, neighbors):
    if neighbors == 0:
        return edges
    selected = set()
    adjacent = defaultdict(list)
    for index, edge in enumerate(edges):
        adjacent[edge["genome_a"]].append((edge["distance"], index))
        adjacent[edge["genome_b"]].append((edge["distance"], index))
    for node in nodes:
        selected.update(index for _, index in sorted(adjacent[node["id"]])[:neighbors])
    # Add a minimum spanning forest so sparsifying displayed edges never disconnects a component.
    parent = {n["id"]: n["id"] for n in nodes}

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for index, edge in sorted(enumerate(edges), key=lambda item: (item[1]["distance"], item[0])):
        a, b = root(edge["genome_a"]), root(edge["genome_b"])
        if a != b:
            parent[a] = b
            selected.add(index)
    return [edges[index] for index in sorted(selected)]


def select_network_references(nodes, edges, limit=10, neighbors=5):
    """Keep each sample's top ANI references and sparse sample/sample links.

    Reference selection is per sample, so a shared reference selected by
    another sample cannot add a lower-ranked edge to this sample. A zero
    limit restores the original unrestricted graph/edge selection.
    """
    if limit < 0:
        raise ValueError("Network reference limit must be nonnegative")
    if limit == 0:
        return nodes, display_edges(nodes, edges, neighbors)
    kinds = {n["id"]: n["kind"] for n in nodes}
    candidates = defaultdict(list)
    sample_edges = []
    for edge in edges:
        if edge.get("status", "accepted") != "accepted":
            continue
        a, b = edge["genome_a"], edge["genome_b"]
        if kinds[a] == kinds[b] == "sample":
            sample_edges.append(edge)
        elif {kinds[a], kinds[b]} == {"sample", "reference"}:
            sample, reference = (a, b) if kinds[a] == "sample" else (b, a)
            candidates[sample].append((edge, reference))
    selected, reference_ids = [], set()
    for sample in sorted(candidates):
        ranked = sorted(candidates[sample], key=lambda item: (
            -float(item[0]["ani_used"]), -float(item[0]["min_aligned_fraction"]), item[1]))
        for edge, reference in ranked[:limit]:
            selected.append(edge)
            reference_ids.add(reference)
    samples = [n for n in nodes if n["kind"] == "sample"]
    selected.extend(display_edges(samples, sample_edges, neighbors))
    visible = [n for n in nodes if n["kind"] == "sample" or n["id"] in reference_ids]
    return visible, selected


def expanded_reference_layout(nodes, positions):
    """Put full reference labels in two nonoverlapping outer columns.

    Only labels move. The same isotropic transform is applied to all points,
    so spreading labels cannot change genome distances or measured ANI.
    """
    references = sorted((n for n in nodes if n["kind"] == "reference"),
                        key=lambda n: (positions[n["id"]]["x"], n["id"]))
    middle = (len(references) + 1) // 2
    groups = [references[:middle], references[middle:]]
    lines = {n["id"]: [n["id"]] + textwrap.wrap(n.get("organism", ""), width=44,
              break_long_words=True, break_on_hyphens=False) for n in references}
    heights = {key: max(18, 14 * len(value)) for key, value in lines.items()}
    column_width = max([270] + [max(map(len, value)) * 6.5 + 30 for value in lines.values()])
    width = 1080 + 2 * column_width + 120
    height = max([800] + [80 + sum(heights[n["id"]] + 8 for n in group) for group in groups])
    xs = [positions[n["id"]]["x"] for n in nodes]
    ys = [positions[n["id"]]["y"] for n in nodes]
    scale = min(1080 / max(max(xs) - min(xs), .1), (height - 120) / max(max(ys) - min(ys), .1))
    xy = {n["id"]: (column_width + 60 + (positions[n["id"]]["x"] - min(xs)) * scale,
                     60 + (max(ys) - positions[n["id"]]["y"]) * scale) for n in nodes}
    labels = []
    for side, group in enumerate(groups):
        group = sorted(group, key=lambda n: (xy[n["id"]][1], n["id"]))
        tops = []
        for i, node in enumerate(group):
            name = node["id"]
            top = max(40, min(xy[name][1] - heights[name] / 2, height - heights[name] - 40))
            if tops:
                previous = group[i - 1]["id"]
                top = max(top, tops[-1] + heights[previous] + 8)
            tops.append(top)
        if group:
            tops[-1] = min(tops[-1], height - heights[group[-1]["id"]] - 40)
            for i in range(len(tops) - 2, -1, -1):
                tops[i] = min(tops[i], tops[i + 1] - heights[group[i]["id"]] - 8)
        for node, top in zip(group, tops):
            labels.append(dict(id=node["id"], x=20 if side == 0 else column_width + 1160,
                               y=top, height=heights[node["id"]], lines=lines[node["id"]], side=side))
    return xy, labels, width, height


def genome_details(node, nodes, edges):
    """Plain text for the click panel; reference strain never becomes sample strain."""
    lines = ["Genome: " + node["id"], "Type: " + node["kind"],
             "Organism: " + (node.get("organism") or "Not provided"),
             "Strain: " + (node.get("strain") or "Not provided"),
             "Isolate: " + (node.get("isolate") or "Not provided")]
    for field, label in (("taxid", "NCBI taxid"), ("bioproject", "BioProject"), ("biosample", "BioSample"),
                         ("assembly_level", "Assembly level"), ("infraspecific_name", "Infraspecific name"),
                         ("metadata_assembly_accession", "Metadata assembly"), ("metadata_match", "Metadata match"),
                         ("length_bp", "Length (bp)"), ("contigs", "Contigs"), ("n50", "N50 (bp)"),
                         ("source_file", "Source FASTA"), ("sequence_status", "Sequence"), ("component", "Layout component")):
        if node.get(field) not in (None, ""):
            lines.append(label + ": " + str(node[field]))
    by_id = {n["id"]: n for n in nodes}
    related = []
    for edge in edges:
        if edge.get("status", "accepted") != "accepted" or node["id"] not in (edge["genome_a"], edge["genome_b"]):
            continue
        other = edge["genome_b"] if node["id"] == edge["genome_a"] else edge["genome_a"]
        if by_id[other]["kind"] != node["kind"]:
            related.append((edge, by_id[other]))
    related.sort(key=lambda item: (-float(item[0]["ani_used"]), -float(item[0]["min_aligned_fraction"]), item[1]["id"]))
    if node["kind"] == "sample":
        for field, label in (("selected_orfs", "Selected ORFs"), ("matched_orfs", "Matched ORFs")):
            if node.get(field) not in (None, ""):
                lines.append(label + ": " + str(node[field]))
        if related:
            edge, reference = related[0]
            lines.extend(["Nearest reference: " + reference["id"],
                          "Reference organism: " + (reference.get("organism") or "Not provided"),
                          "Reference strain: " + (reference.get("strain") or "Not provided"),
                          "Reference isolate: " + (reference.get("isolate") or "Not provided")])
        else:
            lines.append("Nearest reference: No accepted genome ANI comparison")
    else:
        lines.append("Protein hits: %s (representative assembly labels)" % node.get("orf_hits", "Not provided"))
        if related:
            edge, sample = related[0]
            lines.append("Closest input sample: " + sample["id"])
    if related:
        edge = related[0][0]
        lines.extend(["Mean reciprocal genome ANI: %.5f%%" % float(edge["ani_used"]),
                      "Minimum reciprocal aligned fraction: %.3f" % float(edge["min_aligned_fraction"])])
    return "\n".join(lines)


def draw_network(nodes, edges, positions, diagnostics, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    limit = getattr(args, "network_references_per_sample", 10)
    visible, shown = select_network_references(nodes, edges, limit, args.neighbors)
    visible_ids = {n["id"] for n in visible}
    for node in nodes:
        node.update(positions[node["id"]])
        node["displayed"] = node["id"] in visible_ids
    # Keep all coordinates so the helper can expand the reference limit later.
    write_csv(args.outdir / "network_nodes.csv", ["id", "kind", "organism", "source_file", "sequence_status",
              "length_bp", "contigs", "n50", "orf_hits", "matched_orfs", "selected_orfs",
              "component", "x", "y", "displayed"] + ANNOTATION_FIELDS, nodes)
    adjacent = defaultdict(list)
    for edge in edges:
        adjacent[edge["genome_a"]].append(edge)
        adjacent[edge["genome_b"]].append(edge)
    details = {n["id"]: genome_details(n, nodes, adjacent[n["id"]]) for n in visible}
    nodes = visible
    write_csv(args.outdir / "network_edges.csv", PAIR_FIELDS + ["distance"], shown)
    write_csv(args.outdir / "layout_diagnostics.csv", list(diagnostics[0]), diagnostics)
    span_x = max(n["x"] for n in nodes) - min(n["x"] for n in nodes)
    span_y = max(n["y"] for n in nodes) - min(n["y"] for n in nodes)
    figure_height = min(10, max(5, 14 * span_y / max(span_x, .1)))
    figure, axis = plt.subplots(figsize=(15, figure_height))
    for edge in shown:
        a, b = positions[edge["genome_a"]], positions[edge["genome_b"]]
        axis.plot([a["x"], b["x"]], [a["y"], b["y"]], color="#c7ccd1", linewidth=.6, zorder=1)
    # Larger gray markers remain visible behind coincident black sample nodes.
    for kind, color, size in (("reference", "#b0b0b0", 65), ("sample", "#111111", 32)):
        group = [n for n in nodes if n["kind"] == kind]
        axis.scatter([n["x"] for n in group], [n["y"] for n in group], c=color, s=size, label=kind, zorder=3)
        if kind == "sample":
            for node in group:
                axis.annotate(node["id"], (node["x"], node["y"]), xytext=(4, 4), textcoords="offset points", fontsize=6)
    axis.set_aspect("equal", adjustable="datalim")
    axis.set_axis_off()
    axis.legend(loc="upper right", frameon=False)
    figure.suptitle("Genome ANI map", fontsize=18, x=.06, ha="left")
    figure.text(.06, .04, "Distances approximate 100 − ANI within each component. Gaps between disconnected components have no ANI meaning.\n"
                "Exact ANI and aligned fractions: ani_pairs.csv and the interactive HTML. Missing ANI is never treated as zero.", fontsize=9)
    figure.savefig(args.outdir / "ani_network.png", dpi=200, bbox_inches="tight")
    plt.close(figure)

    xy, labels, canvas_width, canvas_height = expanded_reference_layout(nodes, positions)
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %.3f %.3f" role="img" aria-label="Genome ANI map"><rect width="%.3f" height="%.3f" fill="white" pointer-events="none"/><g id="drawing">' %
             (canvas_width, canvas_height, canvas_width, canvas_height)]
    for edge in shown:
        x1, y1 = xy[edge["genome_a"]]
        x2, y2 = xy[edge["genome_b"]]
        title = "%s ↔ %s | ANI %.5f%% | aligned fractions %.3f / %.3f" % (
            edge["genome_a"], edge["genome_b"], edge["ani_used"], edge["aligned_fraction_a_to_b"], edge["aligned_fraction_b_to_a"])
        parts.append('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="#bfc6cc" stroke-width="1" pointer-events="stroke"><title>%s</title></line>' %
                     (x1, y1, x2, y2, html.escape(title)))
    for label in labels:
        x, y = xy[label["id"]]
        end_x = (label["x"] + max(map(len, label["lines"])) * 6.5 + 8 if label["side"] == 0 else label["x"] - 8)
        parts.append('<g class="reference-label" pointer-events="none"><line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="#9fa7ae" stroke-width=".7" stroke-dasharray="3 3"/>' %
                     (x, y, end_x, label["y"] + 8))
        parts.append('<text x="%.3f" y="%.3f" font-size="11" fill="#434b53">' % (label["x"], label["y"] + 11))
        for i, line in enumerate(label["lines"]):
            parts.append('<tspan x="%.3f" dy="%s">%s</tspan>' % (label["x"], "0" if i == 0 else "14", html.escape(line)))
        parts.append('</text></g>')
    for node in sorted(nodes, key=lambda n: n["kind"] == "sample"):
        x, y = xy[node["id"]]
        color = "#111111" if node["kind"] == "sample" else "#b0b0b0"
        detail = details[node["id"]]
        parts.append('<circle tabindex="0" role="button" aria-label="Select genome %s" data-id="%s" data-info="%s" cx="%.3f" cy="%.3f" r="%s" fill="%s" stroke="white" stroke-width="1"><title>%s</title></circle>' % (
            html.escape(node["id"], quote=True), html.escape(node["id"], quote=True), html.escape(detail, quote=True), x, y,
            4 if node["kind"] == "sample" else 6, color, html.escape(detail)))
        if node["kind"] == "sample":
            parts.append('<text class="sample-label" x="%.3f" y="%.3f" font-size="10" fill="#333">%s</text>' % (
                x + 6, y - 5, html.escape(node["id"])))
    parts.append("</g></svg>")
    svg = "".join(parts)
    (args.outdir / "ani_network.svg").write_text(svg)
    options = '<option value="">Find a genome…</option>' + "".join(
        '<option value="%s">%s (%s)</option>' % (html.escape(n["id"], quote=True), html.escape(n["id"]), n["kind"]) for n in nodes)
    page = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Genome ANI map</title><style>body{font:15px system-ui;margin:28px;color:#20242a;background:#fafafa}h1{font-size:28px;font-weight:600;margin-bottom:8px}p{max-width:1050px;line-height:1.5}.controls{display:flex;gap:18px;align-items:center;flex-wrap:wrap}select,button{padding:8px;border:1px solid #ccd1d6;border-radius:5px;background:white}.map-scroll{max-height:75vh;overflow:auto;margin-top:16px;border:1px solid #e3e6e8;background:white}svg{width:100%;min-width:1200px;height:auto;background:white;touch-action:none}text{pointer-events:none}circle{cursor:pointer}circle:focus{outline:none;stroke:#197a9a;stroke-width:3}#info{white-space:pre-line;min-height:45px;max-height:230px;overflow:auto;font-size:13px;margin-top:12px;padding:14px;border:1px solid #d5dce3;border-radius:6px;background:white;color:#334}</style>
<h1>Genome ANI map</h1><p>Black dots: your samples. Gray dots: SprayNPray reference assemblies. Shorter distances approximate higher genome ANI <strong>within each connected component</strong>. Click a dot to show genome details, including strain/isolate when provided. Hover lines for exact ANI and aligned fractions. Scroll to zoom and drag the background to pan.</p>
<div class="controls"><select id="find">__OPTIONS__</select><label><input id="samples" type="checkbox" checked>Sample labels</label><label><input id="references" type="checkbox" checked>Reference labels</label><button id="reset">Reset view</button></div><section id="info" aria-label="Genome details" aria-live="polite">Select a genome dot or use Find a genome to view its details.</section><div class="map-scroll">__SVG__</div>
<p>Reference names are spread into outer columns with dashed leader lines; genome dots and ANI distances are unchanged. Dashed leaders identify labels and do not represent ANI.</p>
<p>Distances target <code>100 − mean reciprocal ANI</code>; 2D geometry is an approximation, not a phylogenetic tree. Separate components are placed side by side for display: their gaps are arbitrary. Missing or low-coverage ANI comparisons create no distance constraints. Only pairs with both directions and aligned fractions ≥ MINAF are used. REFERENCE_POLICY NEIGHBOR_POLICY All accepted pairs inform the layout; full ANI measurements remain in ani_pairs.csv. Identical genomes can overlap. Protein-hit percentages indicate representative assembly labels, not strain abundance.</p>
<script>const svg=document.querySelector('svg'),drawing=document.getElementById('drawing'),original=[0,0,__WIDTH__,__HEIGHT__];let box=original.slice(),drag=null;
function point(e,m){return new DOMPoint(e.clientX,e.clientY).matrixTransform(m||svg.getScreenCTM().inverse())}function view(){svg.setAttribute('viewBox',box.join(' '))}svg.addEventListener('wheel',e=>{e.preventDefault();const p=point(e),k=e.deltaY>0?1.15:1/1.15;box=[p.x-(p.x-box[0])*k,p.y-(p.y-box[1])*k,box[2]*k,box[3]*k];view()},{passive:false});
svg.addEventListener('pointerdown',e=>{if(e.target.closest('circle'))return;const m=svg.getScreenCTM().inverse();drag={p:point(e,m),m,box:box.slice()};svg.setPointerCapture(e.pointerId)});svg.addEventListener('pointermove',e=>{if(!drag)return;const p=point(e,drag.m);box=[drag.box[0]-(p.x-drag.p.x),drag.box[1]-(p.y-drag.p.y),drag.box[2],drag.box[3]];view()});svg.addEventListener('pointerup',()=>drag=null);svg.addEventListener('pointercancel',()=>drag=null);
document.getElementById('reset').onclick=()=>{box=original.slice();view()};for(const [control,kind] of [['samples','sample'],['references','reference']])document.getElementById(control).onchange=e=>document.querySelectorAll('.'+kind+'-label').forEach(t=>t.setAttribute('display',e.target.checked?'':'none'));
function select(id){document.getElementById('info').textContent='Select a genome dot or use Find a genome to view its details.';document.querySelectorAll('circle').forEach(c=>{const active=c.dataset.id===id;c.setAttribute('stroke',active?'#197a9a':'white');c.setAttribute('stroke-width',active?'3':'1');if(active)document.getElementById('info').textContent=c.dataset.info})}document.getElementById('find').onchange=e=>select(e.target.value);document.querySelectorAll('circle').forEach(c=>{const choose=()=>{document.getElementById('find').value=c.dataset.id;select(c.dataset.id)};c.addEventListener('click',choose);c.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();choose()}})});</script></html>'''
    policy = ("Each sample shows up to %s reference connections ranked by highest accepted reciprocal genome ANI; "
              "only references selected by at least one sample are shown. Reference/reference edges are hidden." % limit
              if limit else "All reference nodes are included; displayed lines use the nearest-neighbor/spanning-forest rule across all nodes.")
    neighbor_policy = ("Sample/sample lines use NEIGHBORS nearest neighbors plus a minimum spanning forest (0 means all)."
                       if limit else "Lines use NEIGHBORS nearest neighbors per node plus a minimum spanning forest (0 means all).")
    page = page.replace("REFERENCE_POLICY", policy).replace("NEIGHBOR_POLICY", neighbor_policy)
    page = page.replace("__OPTIONS__", options).replace("__SVG__", svg).replace("__WIDTH__", str(canvas_width)).replace("__HEIGHT__", str(canvas_height)).replace("MINAF", str(args.min_af)).replace("NEIGHBORS", str(args.neighbors))
    (args.outdir / "ani_network.html").write_text(page)


def refresh_existing_analysis(directory, summary_only=False, network_references_per_sample=10, assembly_info=None):
    """Add nearest-reference output and redraw saved coordinates without any ANI."""
    directory = Path(directory).resolve()
    _, nodes = read_csv(directory / "network_nodes.csv")
    _, pairs = read_csv(directory / "ani_pairs.csv")
    if not nodes or len({n["id"] for n in nodes}) != len(nodes):
        raise ValueError("network_nodes.csv must contain unique genome IDs")
    if any(n["kind"] not in ("sample", "reference") for n in nodes):
        raise ValueError("Unknown node kind in network_nodes.csv")
    valid_ids = {n["id"] for n in nodes}
    if any(p["genome_a"] not in valid_ids or p["genome_b"] not in valid_ids for p in pairs):
        raise ValueError("ani_pairs.csv contains a genome missing from network_nodes.csv")
    enrich_reference_nodes(nodes, assembly_info)
    result = nearest_reference_summary(nodes, pairs, directory)
    if summary_only:
        return result
    positions = {}
    for node in nodes:
        x, y = float(node["x"]), float(node["y"])
        if not math.isfinite(x) or not math.isfinite(y):
            raise ValueError("Nonfinite saved graph coordinate")
        positions[node["id"]] = dict(x=x, y=y, component=int(node["component"]))
    edges = []
    for pair in pairs:
        if pair["status"] == "accepted":
            record = dict(pair)
            for field in ("ani_used", "aligned_fraction_a_to_b", "aligned_fraction_b_to_a", "min_aligned_fraction"):
                record[field] = float(record[field])
            record["distance"] = max(100 - identity(record["ani_used"]), .001)
            edges.append(record)
    if (directory / "references.csv").is_file():
        _, references = read_csv(directory / "references.csv")
        counts = {r["assembly"]: r["total_orf_hits"] for r in references}
        for node in nodes:
            if node["kind"] == "reference":
                node["orf_hits"] = counts.get(node["id"], node.get("orf_hits", ""))
    _, diagnostics = read_csv(directory / "layout_diagnostics.csv")
    manifest_path = directory / "run_manifest.json"
    manifest = read_json(manifest_path)
    old = manifest.get("arguments", {})
    args = argparse.Namespace(outdir=directory, min_af=old.get("min_af", .5), neighbors=old.get("neighbors", 5),
                              network_references_per_sample=network_references_per_sample)
    draw_network(nodes, edges, positions, diagnostics, args)
    manifest["helper_version"] = VERSION
    manifest.setdefault("arguments", {})["network_references_per_sample"] = network_references_per_sample
    manifest["nearest_reference_definition"] = "highest accepted reciprocal genome ANI; ANI distance = 100 - ANI"
    write_json(manifest_path, manifest)
    return result


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--version", action="version", version=VERSION)
    p.add_argument("--results", type=Path, required=True, help="Folder with one sample-named SprayNPray .csv per sample")
    p.add_argument("--genomes", type=Path, help="Folder with matching sample FASTAs")
    p.add_argument("--genome-extension", default=".fa", help="Sample FASTA filename suffix")
    p.add_argument("--assembly-info", type=Path, help="Optional NCBI assembly TSV (.gz supported) for reference strain/isolate and other click details")
    p.add_argument("--reference-dir", type=Path, help="Optional local reference FASTAs, named by exact GCA/GCF accession")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="Output directory; completed ANI can be reused")
    p.add_argument("--summary-only", action="store_true", help="Only produce hit summary tables; no ANI, downloads or plots")
    p.add_argument("--no-download", action="store_true", help="Use only local/cached reference genomes")
    p.add_argument("--allow-missing-references", action="store_true", help="Keep references without FASTAs in tables; visible as unlinked nodes only with --network-references-per-sample 0")
    p.add_argument("--top-references-per-sample", type=int, default=0, help="Limit reference selection by hit count per sample; 0 includes every reference")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--fastani", default="fastANI", help="FastANI executable")
    p.add_argument("--datasets", default="datasets", help="NCBI datasets executable for downloading missing references")
    p.add_argument("--frag-len", type=int, default=3000, help="FastANI fragment length")
    p.add_argument("--min-af", type=float, default=0.5, help="Minimum aligned fraction required in both directions for plotted ANI")
    p.add_argument("--min-ani", type=float, default=80, help="Minimum mean ANI used for the map")
    p.add_argument("--network-references-per-sample", type=int, default=10,
                   help="Display each sample's top N references ranked by accepted genome ANI; 0 includes all reference nodes")
    p.add_argument("--neighbors", type=int, default=5, help="Nearest sample/sample edges plus spanning forest (all node types when reference limit is 0); 0 shows all such pairs")
    p.add_argument("--fastani-results", type=Path, help="Reuse existing five-column FastANI output explicitly")
    p.add_argument("--force-ani", action="store_true", help="Recompute ANI instead of using a verified cache")
    p.add_argument("--ani-scope", choices=("samples", "all"), default="samples",
                   help="samples: all sample pairs and reciprocal sample/reference pairs; all: also compare references to each other")
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    if not args.results.is_dir():
        p.error("--results must be an existing directory")
    if args.assembly_info is not None and not args.assembly_info.is_file():
        p.error("--assembly-info must be an existing TSV file")
    if args.reference_dir is not None and not args.reference_dir.is_dir():
        p.error("--reference-dir must be an existing directory")
    if args.threads < 1 or args.frag_len < 100 or args.top_references_per_sample < 0 or args.neighbors < 0 or args.network_references_per_sample < 0:
        p.error("Invalid threads, fragment length, reference limit or neighbor count")
    if not 0 <= args.min_af <= 1 or not 0 <= args.min_ani <= 100:
        p.error("--min-af must be 0 to 1; --min-ani must be 0 to 100")
    if args.fastani_results and args.force_ani:
        p.error("--fastani-results and --force-ani are mutually exclusive")
    if not args.summary_only and (args.genomes is None or not args.genomes.is_dir()):
        p.error("--genomes must be an existing folder unless --summary-only is used")
    if args.outdir.resolve() == args.results.resolve():
        p.error("Use a separate output folder so summary CSVs are not treated as sample inputs")
    args.outdir = args.outdir.resolve()
    args.outdir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(args.outdir / "run.log")])
    try:
        master, overviews, references = summarize(args.results, args.outdir, args.top_references_per_sample)
        LOG.info("Summarized %s samples; %s reference assemblies selected for ANI", len(overviews), len(references))
        if args.summary_only:
            LOG.info("Finished summary tables")
            return 0
        work = args.outdir / "work"
        work.mkdir(exist_ok=True)
        nodes = prepare_nodes(overviews, references, master, args, work)
        LOG.info("Calculating ANI for %s genomes with sequences (scope: %s)", sum(bool(n["staged_file"]) for n in nodes), args.ani_scope)
        directional, stage = calculate_ani(nodes, args, work)
        edges = pairwise_tables(nodes, directional, args)
        nearest_reference_summary(nodes, edges, args.outdir)
        LOG.info("%s reciprocal pairs passed ANI/aligned-fraction thresholds", len(edges))
        positions, diagnostics = embed_network(nodes, edges)
        draw_network(nodes, edges, positions, diagnostics, args)
        write_json(args.outdir / "run_manifest.json", dict(version=VERSION, ani_stage=stage,
                   arguments={k: str(v.resolve()) if isinstance(v, Path) else v for k, v in vars(args).items()},
                   sample_count=len(overviews), reference_count=len(references), accepted_pairs=len(edges), layout=diagnostics))
        LOG.info("Finished: %s", args.outdir / "ani_network.html")
    except (ValueError, OSError, ImportError, subprocess.CalledProcessError, zipfile.BadZipFile) as exc:
        LOG.error("%s", exc)
        LOG.error("Inspect run.log and work/*.log. Completed ANI is cached for the next run.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
