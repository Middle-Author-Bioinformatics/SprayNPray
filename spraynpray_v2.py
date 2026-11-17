#!/usr/bin/env python3
"""Assembly-backed SprayNPray prototype: Prodigal -> DIAMOND -> taxonomy.

Python 3.8+; no third-party Python packages. See README_v2.md for methods,
input formats, output definitions and limitations.
"""

import argparse
import csv
import gzip
import hashlib
import json
import logging
import math
import re
import shlex
import shutil
import statistics
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

VERSION = "2.0.0a3"
LOG = logging.getLogger("spraynpray_v2")
RANKS = ("domain", "phylum", "class", "order", "family", "genus", "species")
HIT_FIELDS = ("qseqid", "sseqid", "pident", "length", "qstart", "qend",
              "sstart", "send", "evalue", "bitscore", "qlen", "slen")


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def fasta(path):
    """Yield (entire header without >, sequence); reject malformed records."""
    header, seq = None, []
    with open_text(path) as handle:
        for line_no, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    if not seq:
                        raise ValueError("Empty FASTA record: " + header)
                    yield header, "".join(seq)
                header, seq = line[1:], []
                if not header.strip():
                    raise ValueError("Empty FASTA header at line %s" % line_no)
            else:
                if header is None:
                    raise ValueError("Sequence before FASTA header at line %s" % line_no)
                seq.append("".join(line.split()))
    if header is not None:
        if not seq:
            raise ValueError("Empty FASTA record: " + header)
        yield header, "".join(seq)


def prepare_contigs(source, target):
    contigs = {}
    with open(target, "w") as out:
        for header, seq in fasta(source):
            name = header.split()[0]
            if name in contigs:
                raise ValueError("Duplicate contig ID (first header field): " + name)
            seq = seq.upper()
            if not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", seq):
                raise ValueError("Contig %s contains non-IUPAC DNA characters" % name)
            canonical = sum(seq.count(base) for base in "ACGT")
            contigs[name] = {"length": len(seq), "gc": (
                (seq.count("G") + seq.count("C")) / canonical if canonical else None)}
            out.write(">%s\n%s\n" % (name, seq))
    if not contigs:
        raise ValueError("No contigs in input FASTA")
    return contigs


@dataclass(frozen=True)
class Orf:
    name: str
    contig: str
    start: int
    end: int
    strand: int
    selected: bool


def select_orfs(proteins, contigs, limit, selected_path):
    """Select the first N by leftmost genomic coordinate, on either strand."""
    coords, by_contig = {}, defaultdict(list)
    for header, _ in fasta(proteins):
        fields = header.split(" # ")
        name = fields[0].split()[0]
        contig, sep, ordinal = name.rpartition("_")
        if not sep or not ordinal.isdigit() or contig not in contigs:
            raise ValueError("Cannot map Prodigal protein ID to input contig: " + name)
        if name in coords:
            raise ValueError("Duplicate protein ID: " + name)
        if len(fields) < 4:
            raise ValueError("Expected Prodigal coordinates in protein header: " + name)
        start, end, strand = map(int, fields[1:4])
        if not 1 <= start <= end <= contigs[contig]["length"] or strand not in (-1, 1):
            raise ValueError("Invalid Prodigal coordinates for " + name)
        coords[name] = (contig, start, end, strand)
        by_contig[contig].append(name)
    selected = set()
    for names in by_contig.values():
        names.sort(key=lambda n: (coords[n][1], coords[n][2], n))
        selected.update(names if limit is None else names[:limit])
    with open(selected_path, "w") as out:
        for header, seq in fasta(proteins):
            if header.split()[0] in selected:
                out.write(">%s\n%s\n" % (header, seq))
    return {name: Orf(name, *data, name in selected) for name, data in coords.items()}


def protein_accession(subject):
    """Recognize WP_... and ref|WP_...| IDs without stripping versions."""
    parts = subject.split("|")
    for tag in ("ref", "gb", "emb", "dbj", "sp", "tr"):
        if tag in parts and parts.index(tag) + 1 < len(parts):
            return parts[parts.index(tag) + 1]
    return subject


@dataclass(frozen=True)
class Hit:
    query: str
    subject: str
    identity: float
    length: int
    qcov: float
    scov: float
    evalue: float
    bitscore: float


def read_hits(path, orfs, args):
    """Reapply filters and keep near-top hits; never depend on file ordering."""
    grouped = defaultdict(dict)
    with open_text(path) as handle:
        for line_no, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != len(HIT_FIELDS):
                raise ValueError("DIAMOND line %s: expected %s fields (%s)" % (
                    line_no, len(HIT_FIELDS), " ".join(HIT_FIELDS)))
            q, s, ident, length, qs, qe, ss, se, ev, score, qlen, slen = fields
            if q not in orfs or not orfs[q].selected:
                raise ValueError("DIAMOND query is not a selected Prodigal ORF: " + q)
            ident, ev, score = float(ident), float(ev), float(score)
            length, qs, qe, ss, se, qlen, slen = map(int, (length, qs, qe, ss, se, qlen, slen))
            if not all(math.isfinite(x) for x in (ident, ev, score)) or not (
                    0 <= ident <= 100 and ev >= 0 and score > 0 and length > 0
                    and 1 <= qs <= qe <= qlen and 1 <= ss <= se <= slen):
                raise ValueError("Invalid DIAMOND values at line %s" % line_no)
            qcov, scov = 100 * (qe - qs + 1) / qlen, 100 * (se - ss + 1) / slen
            if (ident < args.min_identity or ev > args.evalue or
                    qcov < args.query_cover or scov < args.subject_cover):
                continue
            hit = Hit(q, protein_accession(s), ident, length, qcov, scov, ev, score)
            old = grouped[q].get(hit.subject)
            if old is None or (hit.bitscore, -hit.evalue, hit.identity) > (
                    old.bitscore, -old.evalue, old.identity):
                grouped[q][hit.subject] = hit
    result = {}
    for q, subjects in grouped.items():
        best = max(h.bitscore for h in subjects.values())
        result[q] = sorted((h for h in subjects.values()
                            if h.bitscore + 1e-9 >= best * (1 - args.top_percent / 100)),
                           key=lambda h: (-h.bitscore, h.evalue, h.subject))
    return result


def load_protein_map(path, wanted):
    """Stream potentially large maps; retain only matched accessions and all sources."""
    mappings = defaultdict(set)
    with open_text(path) as handle:
        for line_no, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 2:
                raise ValueError("Protein map line %s: expected exactly two TAB-separated columns" % line_no)
            protein, genome = (x.strip() for x in fields)
            if protein.lower() in ("protein", "protein_accession", "protein_id") and genome.lower() in (
                    "genome", "genome_accession", "assembly_accession"):
                continue
            if not protein or not re.fullmatch(r"GC[AF]_\d+\.\d+", genome):
                raise ValueError("Invalid protein/genome accession at map line %s" % line_no)
            if protein in wanted:
                mappings[protein].add(genome)
    return mappings


def load_assemblies(path, wanted):
    """Use NCBI's named TAB columns, never positional offsets from pasted text."""
    rows, paired, header = {}, {}, None
    with open_text(path) as handle:
        for line_no, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("##"):
                continue
            if raw.lstrip("# ").startswith("assembly_accession\t"):
                header = raw.rstrip("\r\n").lstrip("# ").split("\t")
                if not {"assembly_accession", "taxid", "organism_name"}.issubset(header):
                    raise ValueError("Assembly metadata header must include assembly_accession, taxid, organism_name")
                if len(set(header)) != len(header):
                    raise ValueError("Duplicate assembly metadata column names")
                continue
            if raw.startswith("#"):
                continue
            if header is None:
                raise ValueError("Missing TAB-separated assembly_accession header in metadata")
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != len(header):
                raise ValueError("Assembly metadata line %s has %s columns; header has %s. "
                                 "Use the original TSV, not text copied from a terminal." % (
                                     line_no, len(fields), len(header)))
            row = dict(zip(header, fields))
            accession = row["assembly_accession"]
            if accession in wanted:
                if accession in rows and rows[accession] != row:
                    raise ValueError("Conflicting assembly metadata rows for " + accession)
                rows[accession] = row
            partner = row.get("gbrs_paired_asm")
            if partner in wanted and row.get("paired_asm_comp") == "identical":
                if partner in paired and paired[partner] != row:
                    raise ValueError("Conflicting identical paired assemblies for " + partner)
                paired[partner] = row
    if header is None:
        raise ValueError("No assembly metadata header found")
    resolved = {a: (rows[a], "exact") if a in rows else (paired[a], "identical_pair")
                for a in wanted if a in rows or a in paired}
    return resolved, header


def dmp_rows(path):
    with open(path) as handle:
        for raw in handle:
            yield [x.strip() for x in raw.split("|")]


class Taxonomy:
    def __init__(self, directory, wanted):
        directory = Path(directory)
        self.parents, self.ranks, self.names, self.merged = {}, {}, {}, {}
        for fields in dmp_rows(directory / "nodes.dmp"):
            taxid, parent = int(fields[0]), int(fields[1])
            self.parents[taxid], self.ranks[taxid] = parent, fields[2]
        merged = directory / "merged.dmp"
        if merged.exists():
            self.merged = {int(f[0]): int(f[1]) for f in dmp_rows(merged)}
        self._lineages = {}
        needed = set()
        for taxid in wanted:
            needed.update(self.lineage(taxid))
        for fields in dmp_rows(directory / "names.dmp"):
            taxid = int(fields[0])
            if taxid in needed and fields[3] == "scientific name":
                self.names[taxid] = fields[1]

    def resolve(self, taxid):
        seen = set()
        while taxid in self.merged:
            if taxid in seen:
                raise ValueError("Cycle in merged.dmp at taxid %s" % taxid)
            seen.add(taxid)
            taxid = self.merged[taxid]
        return taxid if taxid in self.parents else None

    def lineage(self, taxid):
        taxid = self.resolve(taxid)
        if taxid is None:
            return ()
        if taxid in self._lineages:
            return self._lineages[taxid]
        chain, seen, node = [], set(), taxid
        while True:
            if node in seen or node not in self.parents:
                raise ValueError("Broken/cyclic taxonomy lineage at %s" % node)
            seen.add(node)
            chain.append(node)
            parent = self.parents[node]
            if parent == node:
                if node != 1:
                    raise ValueError("Taxonomy terminates at a non-root node: %s" % node)
                break
            node = parent
        result = tuple(reversed(chain))
        self._lineages[taxid] = result
        return result

    def lca(self, taxids):
        chains = [self.lineage(t) for t in set(taxids)]
        if not chains or any(not c for c in chains):
            return None
        last = None
        for group in zip(*chains):
            if len(set(group)) != 1:
                break
            last = group[0]
        return last

    def label(self, taxid):
        return self.names.get(taxid, "taxid:%s" % taxid) if taxid else "unclassified"

    def ranked_lineage(self, taxid):
        values = {}
        for node in self.lineage(taxid):
            rank = self.ranks[node]
            # New and old NCBI dumps use domain and superkingdom respectively.
            if rank in ("domain", "superkingdom"):
                rank = "domain"
            if rank in RANKS:
                values[rank] = self.label(node)
        return values


def taxid_of(row):
    value = row.get("taxid", "")
    return int(value) if value.isdigit() and int(value) > 0 else None


def infer_orf(hits, mappings, assemblies, taxonomy):
    taxids, missing = set(), False
    for hit in hits:
        genomes = mappings.get(hit.subject, ())
        if not genomes:
            missing = True
        for genome in genomes:
            row = assemblies.get(genome)
            taxid = taxonomy.resolve(taxid_of(row[0])) if row else None
            if taxid is None:
                missing = True
            else:
                taxids.add(taxid)
    # Do not manufacture specificity by silently removing unknown candidates.
    lca = taxonomy.lca(taxids) if taxids and not missing else None
    status = ("no_hit" if not hits else "incomplete_metadata" if missing else
              "root_only" if lca == 1 else "classified" if lca else "unclassified")
    return lca if status == "classified" else None, status


def infer_contig(votes, n_selected, taxonomy, args):
    valid = [v for v in votes if v is not None and v != 1]
    if n_selected == 0:
        return None, 0, "no_orfs"
    if not valid:
        return None, 0, "no_taxonomic_evidence"
    if len(valid) < args.min_tax_orfs:
        return None, 0, "insufficient_orfs"
    if len(valid) / n_selected < args.min_classified_fraction:
        return None, 0, "low_classified_fraction"
    counts = Counter(n for v in valid for n in taxonomy.lineage(v))
    candidates = []
    for taxid, count in counts.items():
        if taxid == 1 or count / len(valid) + 1e-12 < args.consensus:
            continue
        chain = taxonomy.lineage(taxid)
        if args.max_resolution == "species":
            # Exclude strain/subspecies nodes and unranked nodes below species.
            species_nodes = [n for n in chain if taxonomy.ranks[n] == "species"]
            if species_nodes and taxid != species_nodes[-1]:
                continue
        candidates.append((len(chain), taxid))
    if not candidates:
        return None, 0, "discordant_evidence"
    _, prediction = max(candidates)
    return prediction, counts[prediction], "classified"


def csv_writer(path, fields):
    handle = open(path, "w", newline="")
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    return handle, writer


def top_hit_record(orf, hits, mappings, assemblies):
    """One actual protein hit and one explicitly representative source assembly."""
    record = dict(orf=orf.name, contig=orf.contig, start=orf.start, end=orf.end,
                  strand=orf.strand, taxa="no hit", top_hit="no hit",
                  protein_accession="", source_genome_accession="",
                  metadata_assembly_accession="", metadata_match="", pident="",
                  evalue="", bitscore="", candidate_genomes=0,
                  source_genome_is_representative=False, status="no_hit")
    if not hits:
        return record, {}
    hit = min(hits, key=lambda h: (-h.bitscore, h.evalue, -h.identity, h.subject))
    genomes = sorted(mappings.get(hit.subject, ()))
    # A shared WP accession cannot identify its original source strain. Prefer
    # a source with metadata, then accession order, solely for the display.
    genome = min(genomes, key=lambda g: (g not in assemblies, g)) if genomes else ""
    row, match = assemblies.get(genome, ({}, "missing"))
    status = ("missing_protein_mapping" if not genome else
              "missing_assembly_metadata" if not row else "matched")
    record.update(taxa=row.get("organism_name", "unclassified"), top_hit=hit.subject,
                  protein_accession=hit.subject, source_genome_accession=genome,
                  metadata_assembly_accession=row.get("assembly_accession", ""),
                  metadata_match=match, pident=hit.identity, evalue=hit.evalue,
                  bitscore=hit.bitscore, candidate_genomes=len(genomes),
                  source_genome_is_representative=len(genomes) > 1, status=status)
    return record, row


def compact_summary(entries):
    """Keep whole top-hit entries below spreadsheet cell limits in all-ORF mode."""
    reserve = "; ... (%s more ORFs; see top_hits.csv)" % len(entries)
    budget = 30000 - len(reserve.encode("utf-16-le")) // 2
    kept, used = [], 0
    for entry in entries:
        size = len(entry.encode("utf-16-le")) // 2 + (2 if kept else 0)
        if used + size > budget:
            break
        kept.append(entry)
        used += size
    omitted = len(entries) - len(kept)
    text = "; ".join(kept)
    if omitted:
        text += ("; " if text else "") + "... (%s more ORFs; see top_hits.csv)" % omitted
    return text, len(kept), omitted


def write_reports(outdir, contigs, orfs, hits, mappings, assemblies, metadata_fields, taxonomy, args):
    selected = [o for o in orfs.values() if o.selected]
    evidence, votes, by_contig = {}, defaultdict(list), defaultdict(list)
    diagnostics = Counter()
    for orf in selected:
        taxid, status = infer_orf(hits.get(orf.name, []), mappings, assemblies, taxonomy)
        evidence[orf.name] = taxid, status
        votes[orf.contig].append(taxid)
        by_contig[orf.contig].append(orf.name)
        diagnostics[status] += 1

    top_hits = {}
    top_fields = ["orf", "taxa", "top_hit", "contig", "start", "end", "strand",
                  "protein_accession", "source_genome_accession", "metadata_assembly_accession",
                  "metadata_match", "pident", "evalue", "bitscore", "candidate_genomes",
                  "source_genome_is_representative", "status"]
    handle, writer = csv_writer(outdir / "top_hits.csv", top_fields + ["assembly_" + f for f in metadata_fields])
    with handle:
        for contig in contigs:
            for name in sorted(by_contig[contig], key=lambda n: (orfs[n].start, orfs[n].end, n)):
                record, row = top_hit_record(orfs[name], hits.get(name, []), mappings, assemblies)
                top_hits[name] = record
                writer.writerow(dict(record, **{"assembly_" + f: row.get(f, "") for f in metadata_fields}))

    fields = ["orf", "taxa", "top_hit", "contig", "protein_accession", "source_genome_accession",
              "metadata_assembly_accession", "metadata_match", "mapping_status", "pident",
              "query_coverage", "subject_coverage", "evalue", "bitscore", "resolved_taxid"]
    fields += ["assembly_" + f for f in metadata_fields]
    handle, writer = csv_writer(outdir / "matches.csv", fields)
    with handle:
        for orf in selected:
            for hit in hits.get(orf.name, []):
                for genome in sorted(mappings.get(hit.subject, ())) or [""]:
                    row, match_type = assemblies.get(genome, ({}, "missing"))
                    taxid = taxonomy.resolve(taxid_of(row))
                    status = ("missing_protein_mapping" if not genome else
                              "missing_assembly_metadata" if not row else
                              "unknown_taxid" if taxid is None else "mapped")
                    name = row.get("organism_name", "unclassified")
                    record = dict(orf=orf.name, taxa=name,
                                  top_hit="%s [%s] %s" % (hit.subject, genome or "unmapped", name),
                                  contig=orf.contig, protein_accession=hit.subject,
                                  source_genome_accession=genome,
                                  metadata_assembly_accession=row.get("assembly_accession", ""),
                                  metadata_match=match_type, mapping_status=status,
                                  pident=hit.identity, query_coverage=hit.qcov,
                                  subject_coverage=hit.scov, evalue=hit.evalue, bitscore=hit.bitscore,
                                  resolved_taxid=taxid or "")
                    record.update({"assembly_" + f: row.get(f, "") for f in metadata_fields})
                    writer.writerow(record)
                    diagnostics[status + "_match_rows"] += 1

    handle, writer = csv_writer(outdir / "orf_taxonomy.csv", [
        "orf", "contig", "start", "end", "strand", "retained_protein_hits",
        "candidate_genomes", "taxid", "taxon", "rank", "status"])
    with handle:
        for orf in selected:
            taxid, status = evidence[orf.name]
            genomes = set(g for h in hits.get(orf.name, []) for g in mappings.get(h.subject, ()))
            writer.writerow(dict(orf=orf.name, contig=orf.contig, start=orf.start, end=orf.end,
                                 strand=orf.strand, retained_protein_hits=len(hits.get(orf.name, [])),
                                 candidate_genomes=len(genomes), taxid=taxid or "",
                                 taxon=taxonomy.label(taxid), rank=taxonomy.ranks.get(taxid, ""), status=status))

    fields = ["contig", "contig_length", "hits_per_kb", "cov", "GC-content", "Average_AAI",
              "closest_blast_hits", "total_orfs", "selected_orfs", "orfs_with_hits",
              "classified_orfs", "classified_fraction", "predicted_taxid", "predicted_taxon",
              "predicted_rank", "supporting_orfs", "support_fraction", "status",
              "summary_orfs_listed", "summary_orfs_omitted"] + list(RANKS)
    handle, writer = csv_writer(outdir / "contigs.csv", fields)
    totals = Counter(o.contig for o in orfs.values())
    with handle:
        for contig, info in contigs.items():
            names = by_contig[contig]
            names.sort(key=lambda n: (orfs[n].start, orfs[n].end, n))
            matched = [n for n in names if hits.get(n)]
            taxid, support, status = infer_contig(votes[contig], len(names), taxonomy, args)
            identities = []
            descriptions = []
            for n in matched:
                identities.append(top_hits[n]["pident"])
            for n in names:
                hit = top_hits[n]
                if hit["status"] == "no_hit":
                    descriptions.append("%s: no hit" % n)
                else:
                    descriptions.append("%s: %s [%s; %s; %.2f%% identity]" % (
                        n, hit["taxa"], hit["protein_accession"],
                        hit["source_genome_accession"] or "unmapped", hit["pident"]))
            summary, listed, omitted = compact_summary(descriptions)
            classified = sum(v is not None for v in votes[contig])
            record = dict(contig=contig, contig_length=info["length"],
                          hits_per_kb=1000 * len(matched) / info["length"], cov="Unknown",
                          **{"GC-content": info["gc"] if info["gc"] is not None else "NA"},
                          Average_AAI=statistics.mean(identities) if identities else "NA",
                          closest_blast_hits=summary, total_orfs=totals[contig],
                          selected_orfs=len(names), orfs_with_hits=len(matched),
                          classified_orfs=classified, classified_fraction=classified / len(names) if names else 0,
                          predicted_taxid=taxid or "", predicted_taxon=taxonomy.label(taxid),
                          predicted_rank=taxonomy.ranks.get(taxid, ""), supporting_orfs=support,
                          support_fraction=support / classified if classified else 0, status=status,
                          summary_orfs_listed=listed, summary_orfs_omitted=omitted)
            record.update(taxonomy.ranked_lineage(taxid))
            writer.writerow(record)
            diagnostics["contigs_" + status] += 1
    return dict(diagnostics)


def run_command(command, log_path, commands):
    LOG.info("Running: %s", shlex.join(command))
    commands.append(command)
    with open(log_path, "w") as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    try:
        with open(path) as handle:
            value = json.load(handle)
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError):
        return {}


def write_json_atomic(path, value):
    temporary = path.with_name(path.name + ".tmp")
    with open(temporary, "w") as handle:
        json.dump(value, handle, indent=2)
        handle.write("\n")
    temporary.replace(path)


def database_fingerprint(path):
    path = Path(path)
    if not path.is_file():
        path = Path(str(path) + ".dmnd")
    stat = path.stat()
    return dict(path=str(path.resolve()), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def cached_stage(cache, stage, signature, output):
    record = cache.get(stage, {})
    return (isinstance(record, dict) and record.get("complete") is True
            and record.get("signature") == signature and output.is_file()
            and record.get("output_sha256") == file_sha256(output))


def save_stage(cache_path, cache, stage, signature, output):
    cache[stage] = dict(complete=True, signature=signature, output_sha256=file_sha256(output))
    write_json_atomic(cache_path, cache)


def legacy_settings_match(previous, args, fields):
    old = previous.get("arguments", {})
    return isinstance(old, dict) and all(k in old and old[k] == getattr(args, k) for k in fields)


def legacy_file_valid(path, manifest_path):
    return (path.is_file() and manifest_path.is_file()
            and path.stat().st_mtime_ns <= manifest_path.stat().st_mtime_ns)


def require_executable(executable):
    if shutil.which(executable) is None:
        raise ValueError("Executable not found on PATH: " + executable)


ALIGNMENT_SETTINGS = ("evalue", "min_identity", "query_cover", "subject_cover", "top_percent", "sensitivity", "diamond")


def run_predictions_and_alignment(args, contig_path, contig_digest, previous, commands):
    """Resume verified stages, including completed a1/a2 runs, and return parsed hits."""
    outdir = args.outdir
    cache_path, manifest_path = outdir / "stage_cache.json", outdir / "run_manifest.json"
    allow_legacy = not cache_path.exists() and previous.get("version") in ("2.0.0a1", "2.0.0a2")
    cache, stages = read_json(cache_path), {}
    old_contigs = outdir / "input_contigs.fasta"
    legacy_contigs = (allow_legacy and legacy_file_valid(old_contigs, manifest_path)
                      and file_sha256(old_contigs) == contig_digest)
    contig_path.replace(old_contigs)
    proteins = args.proteins or outdir / "prodigal.faa"
    signature = dict(contigs_sha256=contig_digest, meta=args.meta,
                     genetic_code=args.genetic_code, prodigal=args.prodigal)
    if args.proteins is not None:
        stages["prodigal"] = "provided"
    else:
        reuse = cached_stage(cache, "prodigal", signature, proteins)
        legacy_reuse = (legacy_contigs and legacy_file_valid(proteins, manifest_path)
                        and legacy_settings_match(previous, args, ("meta", "genetic_code", "prodigal")))
        if reuse or legacy_reuse:
            LOG.info("Reusing completed Prodigal predictions: %s", proteins)
            stages["prodigal"] = "reused" if reuse else "legacy_reused"
        else:
            require_executable(args.prodigal)
            partial = outdir / "prodigal.faa.partial"
            partial_gff = outdir / "prodigal.gff.partial"
            run_command([args.prodigal, "-i", str(old_contigs), "-a", str(partial),
                         "-o", str(partial_gff), "-f", "gff", "-p",
                         "meta" if args.meta else "single", "-g", str(args.genetic_code)],
                        outdir / "prodigal.log", commands)
            # Validate complete records before promoting the tool's output.
            for _ in fasta(partial):
                pass
            partial.replace(proteins)
            partial_gff.replace(outdir / "prodigal.gff")
            stages["prodigal"] = "ran"
        save_stage(cache_path, cache, "prodigal", signature, proteins)

    selected_path = outdir / "selected_orfs.faa"
    previous_selected = (file_sha256(selected_path) if legacy_contigs
                         and legacy_file_valid(selected_path, manifest_path) else None)
    # Use a distinct target even if --proteins points to selected_orfs.faa.
    selected_temporary = outdir / "selected_orfs.faa.partial"
    # The caller parsed contig lengths while normalizing the input.
    contigs = {h.split()[0]: {"length": len(s)} for h, s in fasta(old_contigs)}
    orfs = select_orfs(proteins, contigs, args.orfs_per_contig, selected_temporary)
    selected_digest = file_sha256(selected_temporary)
    selected_temporary.replace(selected_path)
    n_selected = sum(o.selected for o in orfs.values())
    LOG.info("Selected %s of %s ORFs across %s contigs", n_selected, len(orfs), len(contigs))
    blast = args.diamond_results or outdir / "diamond.tsv"
    if args.diamond_results is not None:
        LOG.info("Using supplied DIAMOND hits: %s", blast)
        stages["diamond"] = "provided"
        return orfs, read_hits(blast, orfs, args), stages

    signature = dict(selected_orfs_sha256=selected_digest,
                     database=database_fingerprint(args.diamond_db),
                     settings={k: getattr(args, k) for k in ALIGNMENT_SETTINGS},
                     fields=list(HIT_FIELDS), max_hsps=1)
    reuse = not args.force_reblast and cached_stage(cache, "diamond", signature, blast)
    legacy_reuse = False
    if (not args.force_reblast and legacy_contigs and previous_selected == selected_digest
            and legacy_file_valid(blast, manifest_path)
            and legacy_settings_match(previous, args, ALIGNMENT_SETTINGS)):
        old_db = previous.get("arguments", {}).get("diamond_db")
        if old_db:
            try:
                legacy_reuse = (database_fingerprint(old_db) == signature["database"]
                                and signature["database"]["mtime_ns"] <= blast.stat().st_mtime_ns)
            except OSError:
                pass
    if reuse or legacy_reuse:
        try:
            hits = read_hits(blast, orfs, args)
        except ValueError as exc:
            LOG.warning("Cached DIAMOND output is invalid (%s); rerunning alignment", exc)
        else:
            LOG.info("Reusing completed DIAMOND hits: %s", blast)
            stages["diamond"] = "reused" if reuse else "legacy_reused"
            save_stage(cache_path, cache, "diamond", signature, blast)
            return orfs, hits, stages

    partial = outdir / "diamond.tsv.partial"
    if n_selected:
        require_executable(args.diamond)
        command = [args.diamond, "blastp", "--db", str(args.diamond_db),
                   "--query", str(selected_path), "--out", str(partial),
                   "--outfmt", "6"] + list(HIT_FIELDS) + [
                       "--threads", str(args.threads), "--evalue", str(args.evalue),
                       "--id", str(args.min_identity), "--query-cover", str(args.query_cover),
                       "--subject-cover", str(args.subject_cover), "--top", str(args.top_percent),
                       "--max-hsps", "1"]
        if args.sensitivity != "default":
            command.append("--" + args.sensitivity)
        run_command(command, outdir / "diamond.log", commands)
    else:
        partial.write_text("")
    hits = read_hits(partial, orfs, args)
    partial.replace(blast)
    save_stage(cache_path, cache, "diamond", signature, blast)
    stages["diamond"] = "ran" if n_selected else "no_orfs"
    return orfs, hits, stages


def orf_limit(value):
    if value.lower() == "all":
        return None
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Use a positive integer or 'all'")
    if number < 1:
        raise argparse.ArgumentTypeError("ORF limit must be positive or 'all'")
    return number


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--version", action="version", version=VERSION)
    p.add_argument("-i", "--contigs", type=Path, required=True, help="Contigs FASTA (optionally .gz)")
    p.add_argument("-d", "--diamond-db", type=Path, help="DIAMOND database (.dmnd or prefix)")
    p.add_argument("--protein-to-genome", type=Path, required=True, help="Two-column protein/genome TSV (optionally .gz)")
    p.add_argument("--assembly-info", type=Path, required=True, help="Original NCBI assembly summary TSV (optionally .gz)")
    p.add_argument("--taxdump", type=Path, required=True, help="Directory containing nodes.dmp and names.dmp")
    p.add_argument("-o", "--outdir", type=Path, required=True, help="Output directory; automatically resumes completed stages")
    p.add_argument("--meta", action="store_true", help="Run Prodigal -p meta instead of single")
    p.add_argument("--orfs-per-contig", type=orf_limit, default=25, metavar="N|all", help="First N ORFs by genomic coordinate")
    p.add_argument("-t", "--threads", type=int, default=8, help="DIAMOND threads (Prodigal is single-threaded)")
    p.add_argument("--genetic-code", type=int, default=11, help="Prodigal translation table")
    p.add_argument("--evalue", type=float, default=1e-15)
    p.add_argument("--min-identity", type=float, default=0, help="Minimum amino-acid identity (percent)")
    p.add_argument("--query-cover", type=float, default=50, help="Minimum query coverage (percent)")
    p.add_argument("--subject-cover", type=float, default=50, help="Minimum subject coverage (percent)")
    p.add_argument("--top-percent", type=float, default=5, help="Retain protein hits within this percent of the best bitscore")
    p.add_argument("--sensitivity", choices=("default", "fast", "mid-sensitive", "sensitive", "very-sensitive", "ultra-sensitive"), default="sensitive")
    p.add_argument("--consensus", type=float, default=0.8, help="Required ORF vote fraction (>0.5 to 1)")
    p.add_argument("--min-tax-orfs", type=int, default=3, help="Minimum ORFs with usable taxonomy per contig")
    p.add_argument("--min-classified-fraction", type=float, default=0.5, help="Minimum usable fraction of selected ORFs")
    p.add_argument("--max-resolution", choices=("species", "taxid"), default="species", help="Allow or exclude predictions below species")
    p.add_argument("--prodigal", default="prodigal", help="Prodigal executable")
    p.add_argument("--diamond", default="diamond", help="DIAMOND executable")
    p.add_argument("--proteins", type=Path, help="Reuse Prodigal .faa with original IDs and coordinate headers")
    p.add_argument("--diamond-results", type=Path, help="Reuse 12-column DIAMOND TSV in the documented field order")
    p.add_argument("--force-reblast", action="store_true", help="Rerun DIAMOND even when a completed compatible alignment exists")
    return p


def validate(args, p):
    if args.force_reblast and args.diamond_results is not None:
        p.error("--force-reblast cannot be combined with --diamond-results")
    if args.threads < 1 or args.min_tax_orfs < 1 or args.genetic_code < 1:
        p.error("threads, min-tax-orfs, and genetic-code must be positive")
    for name in ("min_identity", "query_cover", "subject_cover", "top_percent"):
        if not 0 <= getattr(args, name) <= 100:
            p.error(name.replace("_", "-") + " must be between 0 and 100")
    if not math.isfinite(args.evalue) or args.evalue <= 0:
        p.error("evalue must be finite and positive")
    if not 0.5 < args.consensus <= 1:
        p.error("consensus must be greater than 0.5 and at most 1")
    if not 0 <= args.min_classified_fraction <= 1:
        p.error("min-classified-fraction must be between 0 and 1")
    for path in (args.contigs, args.protein_to_genome, args.assembly_info,
                 args.taxdump / "nodes.dmp", args.taxdump / "names.dmp", args.proteins, args.diamond_results):
        if path is not None and not path.is_file():
            p.error("File not found: " + str(path))
    if args.diamond_results is None:
        if args.diamond_db is None:
            p.error("--diamond-db is required unless --diamond-results is supplied")
        if not args.diamond_db.is_file() and not Path(str(args.diamond_db) + ".dmnd").is_file():
            p.error("DIAMOND database not found: " + str(args.diamond_db))
    if args.outdir.exists():
        if not args.outdir.is_dir():
            p.error("Output path is not a directory")
        recognized = ("stage_cache.json", "run_manifest.json", "input_contigs.fasta", "prodigal.faa", "diamond.tsv")
        if any(args.outdir.iterdir()) and not any((args.outdir / n).is_file() for n in recognized):
            p.error("Nonempty output directory does not contain recognizable SprayNPray V2 run files")


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    validate(args, p)
    args.outdir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(args.outdir / "run.log")])
    commands = []
    try:
        LOG.info("SprayNPray V2 prototype %s", VERSION)
        previous = read_json(args.outdir / "run_manifest.json")
        with tempfile.TemporaryDirectory(prefix=".spraynpray-", dir=args.outdir) as temp:
            contig_path = Path(temp) / "normalized.fasta"
            contigs = prepare_contigs(args.contigs, contig_path)
            orfs, hits, stages = run_predictions_and_alignment(
                args, contig_path, file_sha256(contig_path), previous, commands)
        n_selected = sum(o.selected for o in orfs.values())
        wanted = {h.subject for hs in hits.values() for h in hs}
        LOG.info("Loading source genomes for %s matched protein accessions", len(wanted))
        mappings = load_protein_map(args.protein_to_genome, wanted)
        genomes = {g for gs in mappings.values() for g in gs}
        assemblies, metadata_fields = load_assemblies(args.assembly_info, genomes)
        LOG.info("Resolved assembly metadata for %s/%s source genomes", len(assemblies), len(genomes))
        taxonomy = Taxonomy(args.taxdump, {taxid_of(row) for row, _ in assemblies.values()} - {None})
        diagnostics = write_reports(args.outdir, contigs, orfs, hits, mappings, assemblies,
                                    metadata_fields, taxonomy, args)
        manifest = dict(version=VERSION, arguments={k: str(v.resolve()) if isinstance(v, Path) else v
                                                  for k, v in vars(args).items()},
                        commands=commands, stages=stages, diagnostics=diagnostics, contigs=len(contigs),
                        total_orfs=len(orfs), selected_orfs=n_selected,
                        protein_accessions_with_hits=len(wanted),
                        unmapped_protein_accessions=len(wanted - set(mappings)),
                        missing_assembly_accessions=len(genomes - set(assemblies)))
        write_json_atomic(args.outdir / "run_manifest.json", manifest)
        if wanted - set(mappings) or genomes - set(assemblies):
            LOG.warning("Some hits lack mappings/metadata. Inspect matches.csv and run_manifest.json")
        if not hits:
            LOG.warning("No qualifying protein hits; all contigs are still reported")
        LOG.info("Finished: %s", args.outdir / "contigs.csv")
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        LOG.error("%s", exc)
        LOG.error("Run failed. Inspect run.log and tool logs; rerun the command to resume completed stages.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
