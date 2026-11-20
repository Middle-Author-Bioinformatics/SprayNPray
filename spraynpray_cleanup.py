#!/usr/bin/env python3
"""Plan or remove known intermediates from an archived hybrid_results folder.

Dry-run by default. Does not target final assemblies/graphs, QC reports, ANI
tables, native trees/alignments, sample hit summaries, logs or manifests.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import stat
import time
from collections import Counter
from pathlib import Path

VERSION = "1.0.0"


def regular(path):
    return path.is_file() and not path.is_symlink() and path.stat().st_size > 0


def fasta_ready(path):
    if not regular(path):
        return False
    with path.open('rb') as f:
        return f.read(1) == b'>'


def clean_path(path, root):
    """Refuse every symlink component, including symlinked directories."""
    try:
        path.relative_to(root)
    except ValueError:
        return False
    return path != root and '..' not in path.relative_to(root).parts and not any(p.is_symlink() for p in [path, *path.parents] if p != root.parent)


def completed_tree(folder):
    for tree in list(folder.glob('GToTree_output/*.tre')) + list(folder.glob('runs/attempt_*/*.tre')):
        native = tree.parent
        log = native / 'gtotree-runlog.txt'
        alignment = list(native.glob('Aligned_SCGs*.faa')) + list(native.glob('Aligned_SCGs*.fa'))
        if regular(tree) and regular(log) and any(regular(p) for p in alignment):
            text = log.read_text(errors='replace')
            if 'Done!!' in text and 'Total process runtime:' in text:
                return True
    return False


def candidates(root, profile='safe', remove_reference_cache=False):
    """A narrow allowlist; unknown files and incomplete samples are retained."""
    buckets = []
    for sample in sorted(root.glob('sample_*')):
        if not sample.is_dir() or not clean_path(sample, root):
            continue
        success = False
        for assembler, graph, done in [('unicycler','assembly.gfa','unicycler.done'), ('flye','assembly_graph.gfa','flye.done')]:
            folder = sample / assembler
            if not (fasta_ready(folder / 'assembly.fasta') and ((sample / done).is_file() and not (sample / done).is_symlink())):
                continue
            success = True
            if not regular(folder / graph):
                continue
            if assembler == 'unicycler':
                for path in folder.iterdir():
                    if re.fullmatch(r'\d{3}_.+\.gfa', path.name):
                        buckets.append((path, 'assembler_intermediate_graph'))
                for name in ['spades', 'miniasm_assembly', 'long_read_alignments']:
                    if (folder / name).is_dir():
                        buckets.append((folder / name, 'assembler_work'))
            else:
                for path in folder.iterdir():
                    if path.name in ['00-assembly','10-consensus','20-repeat','30-contigger','40-polishing'] and path.is_dir():
                        buckets.append((path, 'flye_work'))
        if profile == 'compact' and success and regular(sample / 'quast/report.tsv'):
            for name in ['R1.fastq.gz','R2.fastq.gz','ONT.fastq.gz']:
                buckets.append((sample / 'trimmed' / name, 'trimmed_reads'))
        hits = sample / 'spraynpray'
        if profile == 'compact' and all(regular(hits / name) for name in ['top_hits.csv','contigs.csv','run_manifest.json']):
            for name in ['diamond.tsv','input_contigs.fasta','prodigal.faa','selected_orfs.faa','prodigal.gff']:
                buckets.append((hits / name, 'search_intermediate'))
    if profile == 'compact' and all(regular(root / 'cohort_ani' / name) for name in ['ani_pairs.csv','ani_matrix.csv','network_nodes.csv']):
        buckets.append((root / 'cohort_ani/work/genomes', 'ani_staged_genomes'))
        buckets.append((root / 'cohort_ani/work/reference_download.zip', 'download_archive'))
        if remove_reference_cache:
            buckets.append((root / 'cohort_ani/work/references', 'downloaded_reference_cache'))
    for folder in sorted(root.glob('cohort_phylogen*')):
        if folder.is_dir() and clean_path(folder, root) and completed_tree(folder):
            if profile == 'compact' and all(regular(folder / name) for name in ['selected_genomes.csv','labels.tsv']):
                buckets.append((folder / 'genomes', 'phylogeny_staged_genomes'))
    seen = set()
    for path, category in buckets:
        if not clean_path(path, root):
            continue
        if path.is_dir():
            iterator = (Path(base) / name for base, dirs, files in os.walk(path, followlinks=False)
                        for name in files)
        else:
            iterator = [path]
        for file in iterator:
            if file in seen or not clean_path(file, root):
                continue
            try:
                info = file.lstat()
            except FileNotFoundError:
                continue
            if not stat.S_ISREG(info.st_mode):
                continue
            seen.add(file)
            yield dict(path=str(file.relative_to(root)), category=category, bytes=info.st_size,
                       allocated_bytes=getattr(info,'st_blocks',0)*512, inode=info.st_ino, device=info.st_dev,
                       links=info.st_nlink, mtime_ns=info.st_mtime_ns)


def plan(root, profile='safe', remove_reference_cache=False):
    rows = list(candidates(root, profile, remove_reference_cache))
    links = Counter((r['device'],r['inode']) for r in rows)
    counted = set()
    for row in rows:
        key = row['device'],row['inode']
        row['estimated_reclaim_bytes'] = (row['allocated_bytes'] if key not in counted and links[key] == row['links'] else 0)
        counted.add(key)
    return rows


def apply_plan(root, rows):
    deleted, skipped = [], []
    for row in rows:
        path = root / row['path']
        if not clean_path(path, root):
            skipped.append(dict(row, reason='symlink or outside root'))
            continue
        try:
            info = path.lstat()
            if not stat.S_ISREG(info.st_mode) or (info.st_dev,info.st_ino,info.st_size,info.st_mtime_ns) != (
                    row['device'],row['inode'],row['bytes'],row['mtime_ns']):
                skipped.append(dict(row, reason='file changed since planning'))
                continue
            path.unlink()
            deleted.append(row)
        except FileNotFoundError:
            skipped.append(dict(row, reason='already absent'))
    # Only remove empty parent directories of files just deleted.
    parents = sorted({p for r in deleted for p in (root/r['path']).parents if p != root and root in p.parents},
                     key=lambda p:len(p.parts), reverse=True)
    for folder in parents:
        if clean_path(folder, root):
            try:
                folder.rmdir()
            except OSError:
                pass
    return deleted, skipped


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True,help='hybrid_results folder')
    p.add_argument('--profile',choices=['safe','compact'],default='safe',help='compact also removes trimmed reads and disposable search/staging files')
    p.add_argument('--remove-reference-cache',action='store_true',help='With compact: remove downloaded reference FASTAs; future analyses may need downloads')
    p.add_argument('--apply',action='store_true',help='Actually delete the listed files; omitted means dry-run')
    p.add_argument('--report',type=Path,help='Plan TSV (default ROOT/cleanup_plan.tsv)')
    args = p.parse_args(argv)
    root=args.root.expanduser().resolve()
    if not root.is_dir() or root==Path(root.anchor) or not (root/'assemblies_combined').is_dir():
        p.error('Expected a hybrid_results directory containing assemblies_combined')
    if args.remove_reference_cache and args.profile!='compact':
        p.error('--remove-reference-cache requires --profile compact')
    rows=plan(root,args.profile,args.remove_reference_cache)
    report=args.report or root/'cleanup_plan.tsv'
    report.parent.mkdir(parents=True,exist_ok=True)
    with report.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['path','category','bytes','estimated_reclaim_bytes'],extrasaction='ignore',delimiter='\t')
        writer.writeheader();writer.writerows(rows)
    for category in sorted({r['category'] for r in rows}):
        group=[r for r in rows if r['category']==category]
        print('%-30s %6d files  %8.2f GiB logical' % (category,len(group),sum(r['bytes'] for r in group)/2**30))
    print('Estimated disk blocks reclaimable: %.2f GiB (hardlinks accounted for; filesystem compression/snapshots may differ)' % (sum(r['estimated_reclaim_bytes'] for r in rows)/2**30))
    print('Plan:',report)
    if args.profile=='compact':
        print('Compact removes derived trimmed FASTQs from completed/QC-reported samples and analysis staging files. Keep original reads separately; removed stages will need regeneration for reruns.')
    if not args.apply:
        print('DRY RUN: no analysis files deleted. Add --apply to execute this profile.')
        return 0
    deleted,skipped=apply_plan(root,rows)
    audit=root/('cleanup_audit_%s_%s.json' % (time.strftime('%Y%m%d_%H%M%S'),os.getpid()))
    audit.write_text(json.dumps(dict(version=VERSION,root=str(root),profile=args.profile,deleted=deleted,skipped=skipped),indent=2)+'\n')
    print('Deleted %s files; skipped %s. Audit: %s' % (len(deleted),len(skipped),audit))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
