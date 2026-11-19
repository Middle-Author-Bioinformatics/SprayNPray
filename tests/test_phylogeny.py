"""ANI-based selection, local genome staging and GToTree orchestration checks."""

import csv
import gzip
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import spraynpray_compare as compare
import spraynpray_phylogeny as phylo


class PhylogenyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.analysis = self.base / "cohort_ani"
        self.analysis.mkdir()
        self.out = self.base / "tree"
        self.nodes = [dict(id=s, kind="sample") for s in ("sample_a", "sample_b", "unlinked")]
        self.nodes += [dict(id="GCF_%09d.1" % i, kind="reference", organism="Pseudomonas aeruginosa", strain="strain_%s" % i) for i in range(1, 13)]
        self.pairs = []
        for i, ref in enumerate(self.nodes[3:]):
            self.pairs += [dict(genome_a="sample_a", genome_b=ref["id"], ani_used=99-i/10, min_aligned_fraction=.9, status="accepted"),
                           dict(genome_a=ref["id"], genome_b="sample_b", ani_used=97+i/10, min_aligned_fraction=.9, status="accepted")]

    def write_inputs(self):
        for index, node in enumerate(self.nodes):
            path = self.base / ("input_%s.fa.gz" % index)
            with gzip.open(path, "wt") as f:
                f.write(">c\n" + "ACGT" * (index+2) + "\n")
            node["source_file"] = str(path)
        compare.write_csv(self.analysis / "network_nodes.csv", ["id", "kind", "organism", "strain", "source_file"], self.nodes)
        compare.write_csv(self.analysis / "ani_pairs.csv", ["genome_a", "genome_b", "ani_used", "min_aligned_fraction", "status"], self.pairs)

    def command(self, *extra):
        return ["--analysis-dir", str(self.analysis), "--outdir", str(self.out)] + list(extra)

    def test_five_reference_selection_is_per_sample_shared_refs_and_ties(self):
        selected, ranking = phylo.select_genomes(self.nodes, self.pairs, 5)
        self.assertEqual(len(ranking), 10)
        self.assertEqual(len(selected), 13)
        self.assertIn("unlinked", [n["id"] for n in selected])
        self.assertEqual([r["reference"] for r in ranking[:5]], [n["id"] for n in self.nodes[3:8]])
        self.assertEqual([r["reference"] for r in ranking[5:]], [n["id"] for n in reversed(self.nodes[-5:])])
        expected_pairs = {(r["sample"], r["reference"]) for r in ranking}
        shared = [dict(p, ani_used=99, min_aligned_fraction=.8) for p in self.pairs]
        shared[2]["min_aligned_fraction"] = .9
        selected, ranking = phylo.select_genomes(self.nodes, shared, 1)
        self.assertEqual(ranking[0]["reference"], self.nodes[4]["id"])
        self.assertEqual(len(phylo.select_genomes(self.nodes, self.pairs, 0)[0]), len(self.nodes))
        _, selected_edges = compare.select_network_references(self.nodes, self.pairs, limit=5, neighbors=0)
        graph_pairs = {(e["genome_a"], e["genome_b"]) if e["genome_a"].startswith("sample") else
                       (e["genome_b"], e["genome_a"]) for e in selected_edges}
        self.assertEqual(expected_pairs, graph_pairs)
        self.assertEqual(len(selected_edges), 10)

    def test_preparation_stages_compressed_fastas_and_safe_strain_labels(self):
        self.write_inputs()
        self.assertEqual(phylo.main(self.command("--prepare-only", "--hmm-set", "Gammaproteobacteria")), 0)
        _, selected = compare.read_csv(self.out / "selected_genomes.csv")
        self.assertEqual(len(selected), 13)
        self.assertEqual(len({r["tree_label"] for r in selected}), 13)
        self.assertTrue(any("Pseudomonas_aeruginosa" in r["tree_label"] and "strain_" in r["tree_label"] for r in selected))
        for row in selected:
            self.assertTrue((self.out / row["input_fasta"]).read_text().startswith(">c\nACGT"))
        with (self.out / "labels.tsv").open() as f:
            labels = list(csv.reader(f, delimiter="\t"))
        self.assertEqual(len(labels), 13)
        self.assertEqual(labels[0][0], "g000000.fa")
        script = (self.out / "run_gtotree.sh").read_text()
        self.assertIn("-H Gammaproteobacteria", script)
        self.assertIn("-n 1 -M 1", script)
        self.assertFalse((self.out / "run_manifest.json").exists())

    def test_gtotree_command_resume_changed_settings_and_empty_tree_failure(self):
        self.write_inputs()
        calls = []
        def fake(command, log, cwd=None):
            calls.append(command)
            output = cwd / command[command.index("-o") + 1]
            output.mkdir(parents=True)
            with (cwd / "labels.tsv").open() as handle:
                labels = [row[1] for row in csv.reader(handle, delimiter="\t")]
            (output / "Aligned_SCGs_mod_names.tre").write_text("(" + ",".join(label + ":0.1" for label in labels) + ");\n")
            log.write_text("mock GToTree success")
        with patch.object(phylo.shutil, "which", return_value="GToTree"), patch.object(compare, "run_command", side_effect=fake):
            self.assertEqual(phylo.main(self.command("--threads", "4")), 0)
            self.assertEqual(len(calls), 1)
            self.assertEqual(phylo.main(self.command("--threads", "2")), 0)
            self.assertEqual(len(calls), 1)
            self.assertEqual(phylo.main(self.command("--nucleotide", "--tree-program", "IQ-TREE")), 0)
            self.assertEqual(len(calls), 2)
            self.assertIn("-z", calls[-1])
            self.assertIn("IQ-TREE", calls[-1])
            self.assertTrue((self.out / "runs/attempt_001/Aligned_SCGs_mod_names.tre").is_file())
            def empty(command, log, cwd=None):
                output = cwd / command[command.index("-o") + 1]
                output.mkdir(parents=True)
                (output / "Aligned_SCGs_mod_names.tre").write_text("")
            with patch.object(compare, "run_command", side_effect=empty):
                self.assertEqual(phylo.main(self.command("--force")), 1)
            self.assertEqual(compare.read_json(self.out / "run_manifest.json")["status"], "failed")
            self.assertEqual(phylo.main(self.command()), 0)
            self.assertTrue((self.out / "runs/attempt_004/Aligned_SCGs_mod_names.tre").is_file())

    def test_missing_genomes_duplicate_pairs_and_changed_sequences_fail(self):
        self.write_inputs()
        Path(self.nodes[0]["source_file"]).unlink()
        with patch.object(compare, "run_command", side_effect=AssertionError("must not run")):
            self.assertEqual(phylo.main(self.command()), 1)
        with self.assertRaisesRegex(ValueError, "Duplicate ANI"):
            phylo.select_genomes(self.nodes, self.pairs + self.pairs[:1], 5)
        self.write_inputs()
        work = self.analysis / "work"
        work.mkdir()
        compare.write_json(work / "ani_cache.json", dict(signature=dict(genomes=[dict(id="sample_a", sha256="wrong")])))
        self.assertEqual(phylo.main(self.command("--prepare-only")), 1)

    def test_tree_membership_reports_genomes_filtered_by_gtotree(self):
        self.out.mkdir()
        tree = self.out / "example.tre"
        tree.write_text("((sample__a:0.1,reference__b:0.2)0.9:0.1,sample__c:0.3);")
        records = [dict(id=x, kind="sample", tree_label="sample__" + x) for x in ("a", "c", "excluded")]
        self.assertEqual(phylo.report_tree_membership(tree, records, self.out), ["excluded"])
        _, rows = compare.read_csv(self.out / "tree_membership.csv")
        self.assertEqual([r["present_in_tree"] for r in rows], ["True", "True", "False"])

    def test_verified_cached_fasta_is_reused_after_original_path_moves(self):
        cached = self.analysis / "work/genomes/g000000.fa"
        cached.parent.mkdir(parents=True)
        cached.write_text(">a\nACGT\n")
        node = dict(id="a", kind="sample", source_file="/missing/a.fa", cached_index=0, ani_sha256=compare.sha256(cached))
        self.assertEqual(phylo.locate_fasta(node, self.analysis), cached.resolve())
        cached.write_text(">a\nAAAA\n")
        with self.assertRaisesRegex(ValueError, "FASTA unavailable"):
            phylo.locate_fasta(node, self.analysis)

    def test_replacement_genome_folder_and_unrelated_output_protection(self):
        self.write_inputs()
        folder = self.base / "replacement"
        folder.mkdir()
        (folder / "sample_a.fa").write_text(">a\nACGT\n")
        self.assertEqual(phylo.locate_fasta(self.nodes[0], self.analysis, genomes=folder), (folder / "sample_a.fa").resolve())
        self.out.mkdir()
        important = self.out / "user-file.txt"
        important.write_text("keep")
        with self.assertRaises(SystemExit):
            phylo.main(self.command("--prepare-only"))
        self.assertEqual(important.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
