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
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

VERSION = "0.1.0"
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
                             protein=row["protein_accession"], identity=None if no_hit else identity(row["pident"])))
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
                            min_protein_identity=min(values), max_protein_identity=max(values)))
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
              "hit_percent", "selected_orf_percent", "mean_protein_identity", "min_protein_identity", "max_protein_identity"]
    write_csv(outdir / "master_summary.csv", fields, master)
    write_csv(outdir / "sample_overview.csv", list(overviews[0]), overviews)
    write_csv(outdir / "orf_hits.csv", ["sample", "contig", "orf", "assembly", "organism", "protein", "identity"], hits)
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
                LOG.warning("Reference download batch failed; its missing references remain unlinked: %s", exc)
        missing = [a for a in accessions if resolved.get(a) is None]
    if missing:
        write_csv(args.outdir / "missing_references.csv", ["assembly"], [dict(assembly=a) for a in missing])
        if not args.allow_missing_references:
            raise ValueError("Missing exact reference versions; see missing_references.csv. Supply local FASTAs or use --allow-missing-references to show these as unlinked gray nodes.")
        LOG.warning("%s references have no sequence and will be unlinked", len(missing))
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
                          orf_hits=sum(r["orf_hits"] for r in groups)))
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
    signature = dict(genomes=[dict(id=n["id"], sha256=n["sha256"]) for n in nodes if n["staged_file"]],
                     fragment_length=args.frag_len, min_fraction=0.2, executable=args.fastani)
    raw_path, cache_path = work / "fastani.tsv", work / "ani_cache.json"
    cache = read_json(cache_path)
    if args.fastani_results:
        return parse_fastani(args.fastani_results, nodes, work), "provided"
    if (not args.force_ani and cache.get("signature") == signature and raw_path.is_file()
            and cache.get("output_sha256") == sha256(raw_path)):
        LOG.info("Reusing completed all-vs-all ANI results")
        return parse_fastani(raw_path, nodes, work), "reused"
    if shutil.which(args.fastani) is None:
        raise ValueError("FastANI executable not found: " + args.fastani)
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


def draw_network(nodes, edges, positions, diagnostics, args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    shown = display_edges(nodes, edges, args.neighbors)
    for node in nodes:
        node.update(positions[node["id"]])
    write_csv(args.outdir / "network_nodes.csv", ["id", "kind", "organism", "source_file", "sequence_status",
              "length_bp", "contigs", "n50", "component", "x", "y"], nodes)
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

    xs, ys = [n["x"] for n in nodes], [n["y"] for n in nodes]
    span_x, span_y = max(max(xs) - min(xs), .1), max(max(ys) - min(ys), .1)
    scale = min(1080 / span_x, 650 / span_y)
    xy = {n["id"]: (60 + (n["x"] - min(xs)) * scale, 60 + (max(ys) - n["y"]) * scale) for n in nodes}
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 800" role="img" aria-label="Genome ANI map"><rect width="1200" height="800" fill="white" pointer-events="none"/><g id="drawing">']
    for edge in shown:
        x1, y1 = xy[edge["genome_a"]]
        x2, y2 = xy[edge["genome_b"]]
        title = "%s ↔ %s | ANI %.5f%% | aligned fractions %.3f / %.3f" % (
            edge["genome_a"], edge["genome_b"], edge["ani_used"], edge["aligned_fraction_a_to_b"], edge["aligned_fraction_b_to_a"])
        parts.append('<line x1="%.3f" y1="%.3f" x2="%.3f" y2="%.3f" stroke="#bfc6cc" stroke-width="1" pointer-events="stroke"><title>%s</title></line>' %
                     (x1, y1, x2, y2, html.escape(title)))
    for node in sorted(nodes, key=lambda n: n["kind"] == "sample"):
        x, y = xy[node["id"]]
        color = "#111111" if node["kind"] == "sample" else "#b0b0b0"
        detail = "%s (%s)\n%s\nSequence: %s\nComponent: %s" % (
            node["id"], node["kind"], node["organism"], node["sequence_status"], node["component"])
        if node["kind"] == "reference":
            detail += "\nSprayNPray protein-hit count: %s (representative assembly labels)" % node["orf_hits"]
        parts.append('<circle data-id="%s" data-info="%s" cx="%.3f" cy="%.3f" r="%s" fill="%s" stroke="white" stroke-width="1"><title>%s</title></circle>' % (
            html.escape(node["id"], quote=True), html.escape(detail, quote=True), x, y,
            4 if node["kind"] == "sample" else 6, color, html.escape(detail)))
        parts.append('<text class="%s-label" x="%.3f" y="%.3f" font-size="10" fill="#333" %s>%s</text>' % (
            node["kind"], x + 6, y - 5, 'display="none"' if node["kind"] == "reference" else "", html.escape(node["id"])))
    parts.append("</g></svg>")
    svg = "".join(parts)
    (args.outdir / "ani_network.svg").write_text(svg)
    options = '<option value="">Find a genome…</option>' + "".join(
        '<option value="%s">%s (%s)</option>' % (html.escape(n["id"], quote=True), html.escape(n["id"]), n["kind"]) for n in nodes)
    page = '''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Genome ANI map</title><style>body{font:15px system-ui;margin:28px;color:#20242a;background:#fafafa}h1{font-size:28px;font-weight:600;margin-bottom:8px}p{max-width:1050px;line-height:1.5}.controls{display:flex;gap:18px;align-items:center;flex-wrap:wrap}select,button{padding:8px;border:1px solid #ccd1d6;border-radius:5px;background:white}svg{width:100%;height:72vh;background:white;border:1px solid #e3e6e8;margin-top:16px;touch-action:none}text{pointer-events:none}#info{white-space:pre-line;min-height:55px;font-size:13px;margin-top:12px;color:#445}</style>
<h1>Genome ANI map</h1><p>Black dots: your samples. Gray dots: SprayNPray reference assemblies. Shorter distances approximate higher genome ANI <strong>within each connected component</strong>. Hover lines for exact ANI and aligned fractions. Scroll to zoom and drag to pan.</p>
<div class="controls"><select id="find">OPTIONS</select><label><input id="samples" type="checkbox" checked>Sample labels</label><label><input id="references" type="checkbox">Reference labels</label><button id="reset">Reset view</button></div>SVG<div id="info"></div>
<p>Distances target <code>100 − mean reciprocal ANI</code>; 2D geometry is an approximation, not a phylogenetic tree. Separate components are placed side by side for display: their gaps are arbitrary. Missing or low-coverage ANI comparisons create no distance constraints. Only pairs with both directions and aligned fractions ≥ MINAF are used. Displayed lines use NEIGHBORS nearest neighbors per node plus a minimum spanning forest (0 means all); all accepted pairs inform the layout. Identical genomes can overlap. Protein-hit percentages indicate representative assembly labels, not strain abundance.</p>
<script>const svg=document.querySelector('svg'),drawing=document.getElementById('drawing');let box=[0,0,1200,800],drag=null;
function view(){svg.setAttribute('viewBox',box.join(' '))}svg.addEventListener('wheel',e=>{e.preventDefault();const r=svg.getBoundingClientRect(),u=(e.clientX-r.left)/r.width,v=(e.clientY-r.top)/r.height,k=e.deltaY>0?1.15:1/1.15;box=[box[0]+u*box[2]*(1-k),box[1]+v*box[3]*(1-k),box[2]*k,box[3]*k];view()},{passive:false});
svg.addEventListener('pointerdown',e=>{drag=[e.clientX,e.clientY,box.slice()];svg.setPointerCapture(e.pointerId)});svg.addEventListener('pointermove',e=>{if(!drag)return;const r=svg.getBoundingClientRect();box=[drag[2][0]-(e.clientX-drag[0])*drag[2][2]/r.width,drag[2][1]-(e.clientY-drag[1])*drag[2][3]/r.height,drag[2][2],drag[2][3]];view()});svg.addEventListener('pointerup',()=>drag=null);
document.getElementById('reset').onclick=()=>{box=[0,0,1200,800];view()};for(const [control,kind] of [['samples','sample'],['references','reference']])document.getElementById(control).onchange=e=>document.querySelectorAll('.'+kind+'-label').forEach(t=>t.setAttribute('display',e.target.checked?'':'none'));
function select(id){document.querySelectorAll('circle').forEach(c=>{const active=c.dataset.id===id;c.setAttribute('stroke',active?'#197a9a':'white');c.setAttribute('stroke-width',active?'3':'1');if(active)document.getElementById('info').textContent=c.dataset.info})}document.getElementById('find').onchange=e=>select(e.target.value);document.querySelectorAll('circle').forEach(c=>c.addEventListener('click',()=>{document.getElementById('find').value=c.dataset.id;select(c.dataset.id)}));</script></html>'''
    page = page.replace("OPTIONS", options).replace("SVG", svg).replace("MINAF", str(args.min_af)).replace("NEIGHBORS", str(args.neighbors))
    (args.outdir / "ani_network.html").write_text(page)


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--version", action="version", version=VERSION)
    p.add_argument("--results", type=Path, required=True, help="Folder with one sample-named SprayNPray .csv per sample")
    p.add_argument("--genomes", type=Path, help="Folder with matching sample FASTAs")
    p.add_argument("--genome-extension", default=".fa", help="Sample FASTA filename suffix")
    p.add_argument("--reference-dir", type=Path, help="Optional local reference FASTAs, named by exact GCA/GCF accession")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="Output directory; completed ANI can be reused")
    p.add_argument("--summary-only", action="store_true", help="Only produce hit summary tables; no ANI, downloads or plots")
    p.add_argument("--no-download", action="store_true", help="Use only local/cached reference genomes")
    p.add_argument("--allow-missing-references", action="store_true", help="Show references without FASTAs as unlinked gray nodes")
    p.add_argument("--top-references-per-sample", type=int, default=0, help="Limit reference selection by hit count per sample; 0 includes every reference")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--fastani", default="fastANI", help="FastANI executable")
    p.add_argument("--datasets", default="datasets", help="NCBI datasets executable for downloading missing references")
    p.add_argument("--frag-len", type=int, default=3000, help="FastANI fragment length")
    p.add_argument("--min-af", type=float, default=0.5, help="Minimum aligned fraction required in both directions for plotted ANI")
    p.add_argument("--min-ani", type=float, default=80, help="Minimum mean ANI used for the map")
    p.add_argument("--neighbors", type=int, default=5, help="Nearest displayed edges per node, plus spanning forest; 0 displays all accepted pairs")
    p.add_argument("--fastani-results", type=Path, help="Reuse existing five-column FastANI output explicitly")
    p.add_argument("--force-ani", action="store_true", help="Recompute ANI instead of using a verified cache")
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    if not args.results.is_dir():
        p.error("--results must be an existing directory")
    if args.reference_dir is not None and not args.reference_dir.is_dir():
        p.error("--reference-dir must be an existing directory")
    if args.threads < 1 or args.frag_len < 100 or args.top_references_per_sample < 0 or args.neighbors < 0:
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
        LOG.info("Calculating all-vs-all ANI for %s genomes with sequences", sum(bool(n["staged_file"]) for n in nodes))
        directional, stage = calculate_ani(nodes, args, work)
        edges = pairwise_tables(nodes, directional, args)
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
