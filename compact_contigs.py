#!/usr/bin/env python3
"""Compact an existing V2 contigs.csv without rerunning Prodigal or DIAMOND.

Preserves all contig statistics and predictions. If the run's orf_taxonomy.csv
is available beside the source CSV (or passed with --orf-taxonomy), uses those
existing ORF assignments. Otherwise summarizes original organism labels and
candidate counts, without inventing an ORF LCA from organism names.
"""

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

from spraynpray_v2 import compact_summary


def read_csv(path):
    # The legacy summary can have fields many megabytes long.
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            break
        except OverflowError:
            limit //= 10
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        if fields is None:
            raise ValueError("CSV has no header: " + str(path))
        rows = list(reader)
        if any(None in r or any(v is None for v in r.values()) for r in rows):
            raise ValueError("CSV rows do not match the header: " + str(path))
    return fields, rows


def organism_label_summary(text):
    """Summarize old name/accession lists, retaining first-occurrence order."""
    result = []
    for entry in text.split("; ") if text else []:
        orf, separator, candidates = entry.partition(": ")
        if not separator:
            raise ValueError("Unrecognized legacy match entry; use the original V2 CSV")
        labels, genomes = [], set()
        for candidate in candidates.split(" | "):
            match = re.fullmatch(r"(.*?) \[(GC[AF]_\d+\.\d+)\]", candidate)
            if match is None:
                if candidates == "unmapped protein":
                    break
                raise ValueError("Unrecognized legacy candidate for " + orf)
            label, genome = match.groups()
            if label not in labels:
                labels.append(label)
            genomes.add(genome)
        if not labels:
            if candidates == "unmapped protein":
                result.append("%s: unmapped protein [0 candidate genomes]" % orf)
                continue
            raise ValueError("Unrecognized legacy candidate list for " + orf)
        displayed = " | ".join(labels[:3])
        if len(labels) > 3:
            displayed += " | ..."
        result.append("%s: %s [%s organism labels; %s candidate genomes]" % (
            orf, displayed, len(labels), len(genomes)))
    return result


def compact_file(source, output, orf_taxonomy=None):
    if source.resolve() == output.resolve():
        raise ValueError("Use a separate output path to preserve the original detailed CSV")
    if output.exists():
        raise ValueError("Output already exists: " + str(output))
    fields, rows = read_csv(source)
    if "closest_blast_hits" not in fields or "contig" not in fields:
        raise ValueError("Expected contig and closest_blast_hits columns")
    by_contig = defaultdict(list)
    if orf_taxonomy is not None:
        orf_fields, orf_rows = read_csv(orf_taxonomy)
        required = {"contig", "orf", "taxid", "taxon", "rank", "status", "candidate_genomes", "start", "end"}
        if not required.issubset(orf_fields):
            raise ValueError("ORF taxonomy CSV is missing required columns")
        seen = set()
        for row in orf_rows:
            if row["orf"] in seen:
                raise ValueError("Duplicate ORF in taxonomy CSV: " + row["orf"])
            seen.add(row["orf"])
            by_contig[row["contig"]].append(row)
        for row in rows:
            expected = int(row["selected_orfs"])
            if len(by_contig[row["contig"]]) != expected:
                raise ValueError("ORF taxonomy does not match the selected ORF count for contig " + row["contig"])
    new_fields = ["summary_orfs_listed", "summary_orfs_omitted", "summary_basis"]
    fields += [f for f in new_fields if f not in fields]
    for row in rows:
        if orf_taxonomy is not None:
            entries = []
            for orf in sorted(by_contig[row["contig"]], key=lambda o: (int(o["start"]), int(o["end"]), o["orf"])):
                resolution = "taxid %s; %s" % (orf["taxid"], orf["rank"]) if orf["taxid"] else orf["status"]
                entries.append("%s: %s [%s; %s candidate genomes]" % (
                    orf["orf"], orf["taxon"], resolution, orf["candidate_genomes"]))
            basis = "existing_orf_taxonomy"
        else:
            entries = organism_label_summary(row["closest_blast_hits"])
            basis = "original_organism_labels_not_LCA"
        summary, listed, omitted = compact_summary(entries, detail_file="original run files")
        row.update(closest_blast_hits=summary, summary_orfs_listed=listed,
                   summary_orfs_omitted=omitted, summary_basis=basis)
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-i", "--input", type=Path, required=True, help="Original V2 contigs.csv")
    p.add_argument("-o", "--output", type=Path, required=True, help="New compact CSV path")
    p.add_argument("--orf-taxonomy", type=Path, help="Existing orf_taxonomy.csv (auto-detects beside input)")
    args = p.parse_args(argv)
    orfs = args.orf_taxonomy
    if orfs is None and (args.input.parent / "orf_taxonomy.csv").is_file():
        orfs = args.input.parent / "orf_taxonomy.csv"
    try:
        count = compact_file(args.input, args.output, orfs)
    except (ValueError, OSError) as exc:
        p.exit(1, str(exc) + "\n")
    print("Wrote %s contigs to %s" % (count, args.output))
    if orfs is None:
        print("Summarized original organism labels; existing contig predictions were preserved.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
