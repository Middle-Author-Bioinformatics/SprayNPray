"""Archive safety and offline reporting of incomplete cohort evidence."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import spraynpray_cleanup as cleanup
import spraynpray_dashboard as dashboard


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def table(path, rows, fields=None, delimiter=','):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=fields or list(rows[0]), delimiter=delimiter)
        writer.writeheader(); writer.writerows(rows)


def cohort(root):
    for sample in ['sample_a', 'sample_b', 'sample_c']:
        write(root/'assemblies_combined'/f'{sample}.fa', '>contig circular=true\nACGTACGTNN\n>small\nGCGC\n')
        write(root/sample/'unicycler/assembly.fasta', '>contig\nACGT\n')
        write(root/sample/'unicycler/assembly.gfa', 'final graph')
        write(root/sample/'unicycler.done', '')
        write(root/sample/'quast/report.tsv', 'Assembly\tassembly\nN50\t1234\n# misassemblies\t0\n')
    write(root/'sample_flye/flye/assembly.fasta', '>flye\nACGT\n')
    write(root/'sample_flye/flye/assembly_graph.gfa', 'final graph')
    write(root/'sample_flye/flye.done', '')
    (root/'sample_missing').mkdir()
    rows=[dict(genome_a='sample_a',genome_b='sample_b',ani_used='99.9',ani_mean='99.9',min_aligned_fraction='0.9',status='accepted'),
          dict(genome_a='sample_a',genome_b='sample_c',ani_used='97.4',ani_mean='97.4',min_aligned_fraction='0.88',status='accepted'),
          dict(genome_a='sample_b',genome_b='sample_c',ani_used='',ani_mean='99.5',min_aligned_fraction='0.1',status='low_aligned_fraction'),
          dict(genome_a='sample_a',genome_b='GCF_000000001.1',ani_used='99.95',ani_mean='99.95',min_aligned_fraction='0.92',status='accepted')]
    table(root/'cohort_ani/ani_pairs.csv', rows)
    write(root/'cohort_ani/ani_matrix.csv', 'genome,sample_a\nsample_a,100\n')
    table(root/'cohort_ani/network_nodes.csv', [dict(id='GCF_000000001.1',kind='reference',organism='Test organism',strain='x</script><script>alert(1)</script>',isolate='')])
    phylo=root/'cohort_phylogeny_fresh/GToTree_output'
    tree=write(phylo/'GToTree_output.tre', '((sample_sample_a:0.01,sample_sample_b:0.01)95:0.02,(g000002:0.03,reference_GCF_000000001.1:0.02):0.02);')
    write(phylo/'gtotree-runlog.txt', 'Done!!\nTotal process runtime: 1 hours\n')
    write(phylo/'Aligned_SCGs_mod_names.faa', '>sample_sample_a\nACDEFG\n')
    table(phylo.parent/'selected_genomes.csv', [dict(id=s,kind='sample',tree_label='sample__'+s,input_fasta=f'genomes/g{i:06d}.fa') for i,s in enumerate(['sample_a','sample_b','sample_c'])])
    write(phylo.parent/'labels.tsv', 'g000000.fa\tsample_sample_a\n')
    table(phylo/'Genomes_summary_info.tsv', [dict(assembly_id='g000000',label='sample_sample_a',perc_comp='98.5',perc_redund='0')],delimiter='\t')
    return tree


class CleanupTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)/'results'; cohort(self.root)

    def test_safe_allowlist_preserves_final_and_incomplete_outputs(self):
        intermediate=write(self.root/'sample_a/unicycler/001_spades_graph_k021.gfa','intermediate')
        unknown=write(self.root/'sample_a/unicycler/notes.txt','keep me')
        incomplete=write(self.root/'sample_missing/unicycler/001_spades_graph_k021.gfa','unfinished')
        trimmed=write(self.root/'sample_a/trimmed/R1.fastq.gz','reads')
        rows=cleanup.plan(self.root)
        self.assertEqual([r['path'] for r in rows],[str(intermediate.relative_to(self.root))])
        self.assertTrue(intermediate.exists())  # planning is read-only
        deleted,skipped=cleanup.apply_plan(self.root,rows)
        self.assertEqual(len(deleted),1); self.assertFalse(skipped)
        for p in [unknown,incomplete,trimmed,self.root/'sample_a/unicycler/assembly.fasta',self.root/'sample_a/unicycler/assembly.gfa']:
            self.assertTrue(p.exists(),p)

    def test_compact_and_reference_cache_are_separate(self):
        for path in ['sample_a/trimmed/R1.fastq.gz','cohort_ani/work/genomes/g000000.fa',
                     'cohort_ani/work/references/GCF_000000001.1.fna',
                     'cohort_phylogeny_fresh/genomes/g000000.fa','cohort_phylogeny_abandoned/genomes/g000000.fa',
                     'sample_flye/flye/40-polishing/reads.bam']:
            write(self.root/path,'some data')
        paths={r['path'] for r in cleanup.plan(self.root,'compact')}
        self.assertIn('sample_a/trimmed/R1.fastq.gz',paths)
        self.assertIn('cohort_ani/work/genomes/g000000.fa',paths)
        self.assertIn('cohort_phylogeny_fresh/genomes/g000000.fa',paths)
        self.assertNotIn('cohort_phylogeny_abandoned/genomes/g000000.fa',paths)
        self.assertNotIn('cohort_ani/work/references/GCF_000000001.1.fna',paths)
        self.assertIn('cohort_ani/work/references/GCF_000000001.1.fna',{r['path'] for r in cleanup.plan(self.root,'compact',True)})

    def test_symlinks_and_changed_files_are_not_deleted(self):
        target=write(self.root.parent/'external.gfa','external')
        (self.root/'sample_a/unicycler/001_link.gfa').symlink_to(target)
        external_dir=self.root.parent/'external'; write(external_dir/'file','external directory')
        (self.root/'sample_a/unicycler/spades').symlink_to(external_dir,target_is_directory=True)
        self.assertFalse(cleanup.plan(self.root))
        intermediate=write(self.root/'sample_a/unicycler/002_graph.gfa','old')
        rows=cleanup.plan(self.root); intermediate.write_text('changed content')
        deleted,skipped=cleanup.apply_plan(self.root,rows)
        self.assertFalse(deleted);self.assertEqual(len(skipped),1);self.assertTrue(target.exists())
        self.assertFalse(cleanup.clean_path(self.root/'../external.gfa',self.root))

    def test_cli_dry_run_then_apply_has_tsv_plan_and_audit(self):
        intermediate=write(self.root/'sample_a/unicycler/001_intermediate.gfa','x')
        cleanup.main(['--root',str(self.root)])
        self.assertTrue(intermediate.exists())
        self.assertTrue((self.root/'cleanup_plan.tsv').read_text().startswith('path\tcategory\t'))
        cleanup.main(['--root',str(self.root),'--apply'])
        self.assertFalse(intermediate.exists())
        audits=list(self.root.glob('cleanup_audit_*.json'))
        self.assertEqual(len(json.loads(audits[0].read_text())['deleted']),1)

    def test_retained_hardlink_has_no_reclaimable_blocks(self):
        source=write(self.root/'sample_a/unicycler/001_intermediate.gfa','x'*8192)
        os.link(source,self.root/'sample_a/unicycler/assembly_copy.gfa')
        self.assertEqual(cleanup.plan(self.root)[0]['estimated_reclaim_bytes'],0)


class DashboardTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)/'results'; self.tree=cohort(self.root)
        self.out=Path(self.temp.name)/'dashboard'

    def build(self):
        return dashboard.build(argparse.Namespace(root=self.root,outdir=self.out,tree=None,title='Test cohort',group_ani=99))

    @unittest.skipUnless(shutil.which('node'), 'Node.js unavailable for dashboard logic check')
    def test_offline_javascript_controls(self):
        self.build()
        result=subprocess.run(['node',str(Path(__file__).with_name('check_dashboard_ui.js')),str(self.out/'dashboard.html')],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)

    def test_bundle_recovers_flye_and_reports_missing_evidence(self):
        data=self.build(); s={r['sample']:r for r in data['samples']}
        self.assertEqual(len(s),5);self.assertEqual(s['sample_flye']['assembly_status'],'completed')
        self.assertFalse(s['sample_flye']['in_combined']);self.assertFalse(s['sample_flye']['in_tree'])
        self.assertEqual(s['sample_missing']['assembly_status'],'missing')
        self.assertEqual(s['sample_a']['length_bp'],14);self.assertEqual(s['sample_a']['n50'],10)
        self.assertEqual(s['sample_a']['marker_recovery_percent'],98.5)
        self.assertEqual(s['sample_a']['marker_redundancy_percent'],0)
        self.assertEqual(data['accepted_pair_count'],2);self.assertEqual(data['expected_pair_count'],10)
        self.assertIsNone(next(p for p in data['pairs'] if p['a']=='sample_b')['ani'])
        self.assertTrue(s['sample_c']['in_tree']);self.assertEqual(s['sample_a']['nearest_reference'],'GCF_000000001.1')
        with zipfile.ZipFile(self.out.with_suffix('.zip')) as z:
            names=z.namelist();self.assertIn('dashboard/data/phylogeny.tre',names)
            self.assertEqual(gzip.decompress(z.read('dashboard/assemblies/sample_a.fa.gz')),(self.root/'assemblies_combined/sample_a.fa').read_bytes())
        content=(self.out/'dashboard.html').read_text()
        self.assertNotIn('<script>alert(1)</script>',content)
        extracted=re.search(r'<script id="payload" type="application/json">(.*?)</script>',content,re.S).group(1)
        self.assertEqual(json.loads(extracted)['references'][0]['strain'],'x</script><script>alert(1)</script>')
        self.assertNotIn('https://cdn',content);self.assertNotIn('fetch(',content)

    def test_unknown_or_low_coverage_ani_never_becomes_zero(self):
        pairs,nearest,refs=dashboard.read_ani(self.root/'cohort_ani',dashboard.sample_files(self.root))
        self.assertIsNone(pairs[2]['ani']);self.assertEqual(pairs[2]['observed_ani'],99.5)
        groups=dashboard.similarity_groups(['sample_a','sample_b','sample_c','sample_flye'],pairs,99)
        self.assertEqual(groups['sample_a'],groups['sample_b'])
        self.assertNotEqual(groups['sample_a'],groups['sample_c'])
        self.assertEqual(groups['sample_flye'],'Unresolved')
        chain=[dict(a='a',b='b',ani=99.9),dict(a='b',b='c',ani=99.8),dict(a='a',b='c',ani=97)]
        self.assertNotEqual(dashboard.similarity_groups(['a','b','c'],chain,99)['a'],dashboard.similarity_groups(['a','b','c'],chain,99)['c'])

    def test_tree_ambiguity_requires_choice(self):
        write(self.root/'cohort_phylogeny_other/GToTree_output/other.tre','(a:1,b:1);')
        with self.assertRaisesRegex(ValueError,'Multiple trees'):self.build()
        self.assertEqual(dashboard.choose_tree(self.root,self.tree),self.tree)

    def test_quoted_newick_and_invalid_branch_lengths(self):
        tree=dashboard.NewickParser("('sample a':1e-4[ok],'O''Brien':0)95;").parse()
        self.assertEqual([n['label'] for n in dashboard.leaves(tree)],['sample a',"O'Brien"])
        for bad in ['(a:-1,b:1);','(a:nan,b:1);','(a:1,b:1)','(,b:1);']:
            with self.assertRaises(ValueError):dashboard.NewickParser(bad).parse()

    def test_missing_tree_and_ani_still_produces_status_dashboard(self):
        self.tree.unlink();(self.root/'cohort_ani/ani_pairs.csv').unlink()
        data=self.build();self.assertIsNone(data['tree']);self.assertIsNone(data['ani_min'])
        self.assertEqual(data['accepted_pair_count'],0)

    def test_refuses_unowned_output_and_filename_collision(self):
        write(self.out/'notes.txt','unrelated user content')
        with self.assertRaisesRegex(ValueError,'choose a new folder'):self.build()
        self.assertEqual((self.out/'notes.txt').read_text(),'unrelated user content')


if __name__=='__main__':unittest.main()
