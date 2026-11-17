"""Synthetic evidence tests; optional real-tool integration when tools are installed."""

import csv
import gzip
import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("spraynpray_v2", ROOT / "spraynpray_v2.py")
v2 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = v2
SPEC.loader.exec_module(v2)


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.dump = self.base / "taxdump"
        self.dump.mkdir()
        # Two species in one genus, one strain under species A, another genus.
        nodes = [(1, 1, "no rank"), (2, 1, "domain"), (10, 2, "genus"),
                 (11, 10, "species"), (12, 10, "species"), (13, 11, "no rank"),
                 (20, 2, "genus"), (21, 20, "species")]
        (self.dump / "nodes.dmp").write_text("".join(
            "%s\t|\t%s\t|\t%s\t|\n" % n for n in nodes))
        labels = {1: "root", 2: "Bacteria", 10: "Pseudomonas", 11: "Pseudomonas alpha",
                  12: "Pseudomonas beta", 13: "Pseudomonas alpha strain 1",
                  20: "Othergenus", 21: "Othergenus gamma"}
        (self.dump / "names.dmp").write_text("".join(
            "%s\t|\t%s\t|\t\t|\tscientific name\t|\n" % n for n in labels.items()))
        (self.dump / "merged.dmp").write_text("99\t|\t11\t|\n")
        self.tax = v2.Taxonomy(self.dump, {11, 12, 13, 21, 99, 999})
        self.args = SimpleNamespace(consensus=.8, min_tax_orfs=3,
                                    min_classified_fraction=.5, max_resolution="species",
                                    min_identity=0, evalue=1e-15, query_cover=50,
                                    subject_cover=50, top_percent=5)

    def hit(self, subject="WP_1.1", score=100, query="a_1", identity=95):
        return v2.Hit(query, subject, identity, 100, 100, 100, 1e-30, score)

    def assembly_file(self, rows):
        path = self.base / "assemblies.tsv"
        fields = ["assembly_accession", "taxid", "organism_name", "gbrs_paired_asm",
                  "paired_asm_comp", "group", "bioproject", "biosample"]
        with path.open("w", newline="") as out:
            out.write("## A real comment\n#" + "\t".join(fields) + "\n")
            for row in rows:
                out.write("\t".join(row) + "\n")
        return path

    def test_shared_protein_preserves_all_genomes(self):
        mapping = self.base / "map.tsv.gz"
        with gzip.open(mapping, "wt") as out:
            out.write("protein_accession\tgenome_accession\n"
                      "WP_1.1\tGCF_000001.1\nWP_1.1\tGCF_000002.1\n"
                      "WP_1.1\tGCF_000001.1\nWP_2.1\tGCF_000003.1\n")
        sources = v2.load_protein_map(mapping, {"WP_1.1"})
        self.assertEqual(sources, {"WP_1.1": {"GCF_000001.1", "GCF_000002.1"}})
        assemblies = {"GCF_000001.1": ({"taxid": "11"}, "exact"),
                      "GCF_000002.1": ({"taxid": "12"}, "exact")}
        self.assertEqual(v2.infer_orf([self.hit()], sources, assemblies, self.tax), (10, "classified"))

    def test_missing_candidate_does_not_create_specificity(self):
        self.assertEqual(v2.infer_orf([self.hit()], {"WP_1.1": {"known", "missing"}},
                                     {"known": ({"taxid": "11"}, "exact")}, self.tax),
                         (None, "incomplete_metadata"))
        self.assertEqual(v2.infer_orf([self.hit()], {}, {}, self.tax), (None, "incomplete_metadata"))
        self.assertEqual(v2.infer_orf([], {}, {}, self.tax), (None, "no_hit"))

    def test_consensus_and_evidence_coverage(self):
        self.assertEqual(v2.infer_contig([11, 11, 11, 11, 10], 5, self.tax, self.args), (11, 4, "classified"))
        self.assertEqual(v2.infer_contig([11, 11, 11, 10, 10], 5, self.tax, self.args), (10, 5, "classified"))
        self.assertEqual(v2.infer_contig([11, 11], 2, self.tax, self.args)[2], "insufficient_orfs")
        self.assertEqual(v2.infer_contig([11] * 3 + [None] * 7, 10, self.tax, self.args)[2], "low_classified_fraction")
        self.assertEqual(v2.infer_contig([], 0, self.tax, self.args)[2], "no_orfs")
        self.assertEqual(v2.infer_contig([None], 1, self.tax, self.args)[2], "no_taxonomic_evidence")

    def test_species_cap_and_broad_votes_do_not_support_descendants(self):
        self.assertEqual(v2.infer_contig([13] * 3, 3, self.tax, self.args), (11, 3, "classified"))
        self.args.max_resolution = "taxid"
        self.assertEqual(v2.infer_contig([13] * 3, 3, self.tax, self.args), (13, 3, "classified"))
        self.assertEqual(v2.infer_contig([10] * 3, 3, self.tax, self.args), (10, 3, "classified"))

    def test_taxdump_merged_missing_and_domain(self):
        self.assertEqual(self.tax.resolve(99), 11)
        self.assertIsNone(self.tax.resolve(999))
        self.assertEqual(self.tax.lca([99, 12]), 10)
        self.assertEqual(self.tax.ranked_lineage(13)["domain"], "Bacteria")
        self.tax.ranks[2] = "superkingdom"
        self.assertEqual(self.tax.ranked_lineage(13)["domain"], "Bacteria")
        self.assertIsNone(self.tax.lca([11, 999]))

    def test_metadata_pairs_exact_precedence_and_csv_commas(self):
        path = self.assembly_file([
            ("GCA_000001.1", "11", "Species, strain X", "GCF_000001.1", "identical", "bacteria", "P1", "S1"),
            ("GCA_000002.1", "12", "Other", "GCF_000002.1", "different", "bacteria", "P2", "S2"),
            ("GCF_000003.1", "12", "Exact", "na", "na", "bacteria", "P3", "S3")])
        rows, fields = v2.load_assemblies(path, {"GCF_000001.1", "GCF_000002.1", "GCF_000003.1"})
        self.assertEqual(rows["GCF_000001.1"][1], "identical_pair")
        self.assertEqual(rows["GCF_000001.1"][0]["organism_name"], "Species, strain X")
        self.assertNotIn("GCF_000002.1", rows)
        self.assertEqual(rows["GCF_000003.1"][1], "exact")
        self.assertIn("bioproject", fields)

    def test_malformed_metadata_rejected(self):
        path = self.assembly_file([])
        with path.open("a") as out:
            out.write("GCF_000001.1\t11\tToo few fields\n")
        with self.assertRaisesRegex(ValueError, "columns"):
            v2.load_assemblies(path, set())

    def test_first_25_selection_uses_coordinate_and_keeps_underscores(self):
        path = self.base / "genes.faa"
        # Deliberately reversed input record order, alternating strands.
        path.write_text("".join(
            ">contig_with_under_%s # %s # %s # %s # ID=1_%s\nMKK\n" %
            (i, i * 10, i * 10 + 8, -1 if i % 2 else 1, i) for i in range(30, 0, -1)))
        contigs = {"contig_with_under": {"length": 400}}
        dest = self.base / "selected.faa"
        orfs = v2.select_orfs(path, contigs, 25, dest)
        self.assertEqual(sum(o.selected for o in orfs.values()), 25)
        self.assertTrue(orfs["contig_with_under_25"].selected)
        self.assertFalse(orfs["contig_with_under_26"].selected)
        self.assertEqual(len(list(v2.fasta(dest))), 25)
        self.assertEqual(sum(o.selected for o in v2.select_orfs(path, contigs, None, dest).values()), 30)

    def test_alignment_filters_ties_and_order_independence(self):
        blast = self.base / "diamond.tsv"
        rows = [
            "a_1\tref|WP_1.1|\t90\t100\t1\t100\t1\t100\t1e-30\t100\t100\t100",
            "a_1\tWP_2.1\t85\t100\t1\t100\t1\t100\t1e-30\t95\t100\t100",
            "a_1\tWP_3.1\t85\t100\t1\t100\t1\t100\t1e-30\t94\t100\t100",
            "a_1\tWP_4.1\t99\t20\t1\t20\t1\t20\t1e-30\t150\t100\t100",
            "a_1\tWP_5.1\t99\t100\t1\t100\t1\t100\t1e-3\t150\t100\t100"]
        orfs = {"a_1": v2.Orf("a_1", "a", 1, 300, 1, True)}
        blast.write_text("\n".join(reversed(rows + rows[:1])) + "\n")
        hits = v2.read_hits(blast, orfs, self.args)
        self.assertEqual([h.subject for h in hits["a_1"]], ["WP_1.1", "WP_2.1"])
        blast.write_text(rows[0].replace("a_1", "not_selected") + "\n")
        with self.assertRaisesRegex(ValueError, "not a selected"):
            v2.read_hits(blast, orfs, self.args)

    def fixture_inputs(self):
        contigs = self.base / "contigs.fna"
        contigs.write_text(">a description\n" + "ATG" * 400 + "\n>b\n" + "N" * 100 + "\n")
        proteins = self.base / "proteins.faa"
        proteins.write_text("".join(
            ">a_%s # %s # %s # 1 # ID=1_%s\n%s\n" %
            (i, 1 + (i - 1) * 300, i * 300, i, "M" * 100) for i in range(1, 4)))
        blast = self.base / "diamond.tsv"
        blast.write_text("".join(
            "a_%s\tWP_1.1\t95\t100\t1\t100\t1\t100\t1e-30\t100\t100\t100\n" % i
            for i in range(1, 4)))
        mapping = self.base / "protein_to_genome.tsv"
        mapping.write_text("WP_1.1\tGCF_000001.1\nWP_1.1\tGCF_000002.1\n")
        metadata = self.assembly_file([
            ("GCF_000001.1", "11", "Species, alpha", "na", "na", "bacteria", "P1", "S1"),
            ("GCF_000002.1", "12", "Species beta", "na", "na", "bacteria", "P2", "S2")])
        return contigs, proteins, blast, mapping, metadata

    def test_cli_reuse_reports_and_empty_alignment(self):
        contigs, proteins, blast, mapping, metadata = self.fixture_inputs()
        args = ["-i", str(contigs), "--proteins", str(proteins), "--diamond-results", str(blast),
                "--protein-to-genome", str(mapping), "--assembly-info", str(metadata), "--taxdump", str(self.dump)]
        out = self.base / "output with spaces"
        result = subprocess.run([sys.executable, str(ROOT / "spraynpray_v2.py")] + args + ["-o", str(out)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (out / "contigs.csv").open() as f:
            rows = list(csv.DictReader(f))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["predicted_taxid"], "10")
        self.assertEqual(rows[1]["status"], "no_orfs")
        with (out / "matches.csv").open() as f:
            matches = list(csv.DictReader(f))
        self.assertEqual(len(matches), 6)
        self.assertEqual(matches[0]["taxa"], "Species, alpha")
        self.assertEqual(matches[0]["assembly_bioproject"], "P1")
        manifest = json.loads((out / "run_manifest.json").read_text())
        self.assertEqual(manifest["selected_orfs"], 3)
        # Empty alignment is a valid no-hit run, with every contig retained.
        blast.write_text("")
        empty_out = self.base / "empty"
        result = subprocess.run([sys.executable, str(ROOT / "spraynpray_v2.py")] + args + ["-o", str(empty_out)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (empty_out / "contigs.csv").open() as f:
            rows = list(csv.DictReader(f))
        self.assertEqual(rows[0]["status"], "no_taxonomic_evidence")
        # A recognized output directory can regenerate reports on a repeat run.
        again = subprocess.run([sys.executable, str(ROOT / "spraynpray_v2.py")] + args + ["-o", str(out)], capture_output=True, text=True)
        self.assertEqual(again.returncode, 0, again.stderr)

    def test_one_actual_top_hit_per_selected_orf(self):
        contigs, proteins, blast, _, _ = self.fixture_inputs()
        out = self.base / "reports"
        out.mkdir()
        info = v2.prepare_contigs(contigs, out / "input.fasta")
        orfs = v2.select_orfs(proteins, info, 25, out / "selected.faa")
        hits = {"a_1": [self.hit("WP_2.1", score=95), self.hit("WP_1.1", score=100)],
                "a_2": [self.hit("WP_2.1", score=100, query="a_2", identity=90),
                        self.hit("WP_1.1", score=100, query="a_2", identity=95)]}
        genomes = {"GCF_%09d.1" % i for i in range(1000)}
        assemblies = {g: ({"assembly_accession": g, "taxid": "11",
                           "organism_name": "Pseudomonas alpha"}, "exact") for g in genomes}
        v2.write_reports(out, info, orfs, hits, {"WP_1.1": genomes, "WP_2.1": genomes}, assemblies,
                         ["assembly_accession", "taxid", "organism_name"], self.tax, self.args)
        with (out / "top_hits.csv").open() as f:
            top = list(csv.DictReader(f))
        self.assertEqual(len(top), 3)
        self.assertEqual([r["orf"] for r in top], ["a_1", "a_2", "a_3"])
        self.assertEqual([r["protein_accession"] for r in top], ["WP_1.1", "WP_1.1", ""])
        self.assertEqual(top[0]["source_genome_accession"], "GCF_000000000.1")
        self.assertEqual(top[0]["candidate_genomes"], "1000")
        self.assertEqual(top[0]["source_genome_is_representative"], "True")
        self.assertEqual(top[2]["status"], "no_hit")
        with (out / "contigs.csv").open() as f:
            row = next(csv.DictReader(f))
        summary = row["closest_blast_hits"]
        self.assertEqual(summary.count("% identity]"), 2)
        self.assertNotIn("WP_2.1", summary)
        self.assertEqual(summary.count("GCF_"), 2)
        self.assertIn("a_3: no hit", summary)
        self.assertLess(len(summary), 500)
        with (out / "matches.csv").open() as f:
            self.assertEqual(sum(1 for _ in csv.DictReader(f)), 4000)

    def test_top_hit_ties_and_missing_mapping(self):
        orf = v2.Orf("a_1", "a", 1, 300, 1, True)
        a = self.hit("WP_1.1", identity=90)
        b = self.hit("WP_2.1", identity=95)
        record, _ = v2.top_hit_record(orf, [a, b], {}, {})
        self.assertEqual(record["protein_accession"], "WP_2.1")
        self.assertEqual(record["status"], "missing_protein_mapping")
        self.assertEqual(v2.top_hit_record(orf, [b, a], {}, {})[0], record)
        record, _ = v2.top_hit_record(orf, [self.hit("WP_2.1"), self.hit("WP_1.1")], {}, {})
        self.assertEqual(record["protein_accession"], "WP_1.1")
        sources = {"WP_1.1": {"GCF_000001.1", "GCF_000002.1"}}
        record, _ = v2.top_hit_record(orf, [a], sources, {"GCF_000002.1": ({"organism_name": "known"}, "exact")})
        self.assertEqual(record["source_genome_accession"], "GCF_000002.1")
        self.assertEqual(record["taxa"], "known")

    def test_all_orf_summary_is_bounded(self):
        entries = ["a_%s: Pseudomonas alpha [WP_1.1; GCF_000001.1; 99.00%% identity]" % i for i in range(6000)]
        summary, listed, omitted = v2.compact_summary(entries)
        self.assertLessEqual(len(summary.encode("utf-16-le")) // 2, 30000)
        self.assertEqual(listed + omitted, 6000)
        self.assertGreater(omitted, 0)
        self.assertIn("top_hits.csv", summary)

    def stage_fixture(self):
        contigs, proteins, _, mapping, metadata = self.fixture_inputs()
        db = self.base / "database.dmnd"
        db.write_bytes(b"test database")
        out = self.base / "stage_run"
        out.mkdir()
        args = v2.parser().parse_args([
            "-i", str(contigs), "-d", str(db), "--protein-to-genome", str(mapping),
            "--assembly-info", str(metadata), "--taxdump", str(self.dump), "-o", str(out), "--meta"])
        calls = []

        def fake_command(command, log_path, commands):
            calls.append(command)
            commands.append(command)
            if "-a" in command:
                Path(command[command.index("-a") + 1]).write_text(proteins.read_text())
                Path(command[command.index("-o") + 1]).write_text("##gff-version 3\n")
            else:
                selected = Path(command[command.index("--query") + 1])
                target = Path(command[command.index("--out") + 1])
                target.write_text("".join(
                    "%s\tWP_1.1\t95\t100\t1\t100\t1\t100\t1e-30\t100\t100\t100\n" % h.split()[0]
                    for h, _ in v2.fasta(selected)))

        def run(previous=None):
            normalized = out / "request.fasta"
            v2.prepare_contigs(contigs, normalized)
            with patch.object(v2, "run_command", side_effect=fake_command), patch.object(v2, "require_executable"):
                return v2.run_predictions_and_alignment(args, normalized, v2.file_sha256(normalized), previous or {}, [])
        return args, calls, run

    def test_resume_validates_inputs_parameters_database_and_tampering(self):
        args, calls, run = self.stage_fixture()
        self.assertEqual(run()[2], {"prodigal": "ran", "diamond": "ran"})
        self.assertEqual(len(calls), 2)
        args.consensus = 1.0
        args.threads = 16
        self.assertEqual(run()[2], {"prodigal": "reused", "diamond": "reused"})
        self.assertEqual(len(calls), 2)
        args.orfs_per_contig = 2
        self.assertEqual(run()[2]["diamond"], "ran")
        self.assertEqual(len(calls), 3)
        args.query_cover = 60
        self.assertEqual(run()[2]["diamond"], "ran")
        args.diamond_db.write_bytes(b"changed database")
        self.assertEqual(run()[2]["diamond"], "ran")
        (args.outdir / "diamond.tsv").write_text("truncated\n")
        self.assertEqual(run()[2]["diamond"], "ran")
        args.force_reblast = True
        self.assertEqual(run()[2]["diamond"], "ran")
        args.force_reblast = False
        with args.contigs.open("a") as f:
            f.write(">new_contig\n" + "ATG" * 100 + "\n")
        self.assertEqual(run()[2]["prodigal"], "ran")

    def test_legacy_completed_run_is_adopted_without_external_tools(self):
        args, calls, run = self.stage_fixture()
        run()
        (args.outdir / "stage_cache.json").unlink()
        previous = dict(version="2.0.0a1", arguments={
            k: str(v.resolve()) if isinstance(v, Path) else v for k, v in vars(args).items()})
        v2.write_json_atomic(args.outdir / "run_manifest.json", previous)
        self.assertEqual(run(previous)[2], {"prodigal": "legacy_reused", "diamond": "legacy_reused"})
        self.assertEqual(len(calls), 2)
        self.assertTrue((args.outdir / "stage_cache.json").is_file())
        self.assertEqual(run()[2], {"prodigal": "reused", "diamond": "reused"})

    def test_partial_unmarked_and_empty_completed_alignments(self):
        args, calls, run = self.stage_fixture()
        (args.outdir / "diamond.tsv.partial").write_text("unfinished")
        (args.outdir / "diamond.tsv").write_text("unmarked")
        self.assertEqual(run()[2]["diamond"], "ran")
        self.assertEqual(len(calls), 2)
        blast = args.outdir / "diamond.tsv"
        cache_path = args.outdir / "stage_cache.json"
        cache = v2.read_json(cache_path)
        blast.write_text("")
        v2.save_stage(cache_path, cache, "diamond", cache["diamond"]["signature"], blast)
        self.assertEqual(run()[1], {})
        self.assertEqual(len(calls), 2)

    def test_failed_alignment_does_not_become_a_completed_cache(self):
        args, calls, run = self.stage_fixture()
        normalized = args.outdir / "request.fasta"
        v2.prepare_contigs(args.contigs, normalized)
        _, proteins, _, _, _ = self.fixture_inputs()

        def fail_alignment(command, log_path, commands):
            if "-a" in command:
                Path(command[command.index("-a") + 1]).write_text(proteins.read_text())
                Path(command[command.index("-o") + 1]).write_text("##gff-version 3\n")
            else:
                Path(command[command.index("--out") + 1]).write_text("incomplete")
                raise subprocess.CalledProcessError(1, command)

        with patch.object(v2, "run_command", side_effect=fail_alignment), patch.object(v2, "require_executable"):
            with self.assertRaises(subprocess.CalledProcessError):
                v2.run_predictions_and_alignment(args, normalized, v2.file_sha256(normalized), {}, [])
        self.assertNotIn("diamond", v2.read_json(args.outdir / "stage_cache.json"))
        self.assertFalse((args.outdir / "diamond.tsv").exists())
        self.assertEqual(run()[2], {"prodigal": "reused", "diamond": "ran"})
        self.assertEqual(len(calls), 1)

    def test_duplicate_contig_ids_rejected(self):
        path = self.base / "duplicate.fasta"
        path.write_text(">a first\nACGT\n>a second\nACGT\n")
        with self.assertRaisesRegex(ValueError, "Duplicate contig"):
            v2.prepare_contigs(path, self.base / "normalized.fasta")

    @unittest.skipUnless(shutil.which("prodigal") and shutil.which("diamond"), "requires real Prodigal and DIAMOND")
    def test_real_prodigal_diamond_pipeline(self):
        contigs = ROOT / "toy_dataset" / "test_dataset.fasta"
        protein_path = self.base / "reference.faa"
        subprocess.run(["prodigal", "-i", str(contigs), "-a", str(protein_path),
                        "-o", str(self.base / "reference.gff"), "-f", "gff", "-p", "meta"],
                       check=True, capture_output=True)
        reference = self.base / "database.faa"
        protein_records = list(v2.fasta(protein_path))
        self.assertGreater(len(protein_records), 25)
        reference.write_text("".join(
            ">WP_%09d.1\n%s\n" % (i, seq) for i, (_, seq) in enumerate(protein_records, 1)))
        mapping = self.base / "map.tsv"
        mapping.write_text("".join("WP_%09d.1\tGCF_000001.1\n" % i
                                   for i in range(1, len(protein_records) + 1)))
        metadata = self.assembly_file([
            ("GCF_000001.1", "11", "Pseudomonas alpha", "na", "na", "bacteria", "P1", "S1")])
        db = self.base / "database"
        subprocess.run(["diamond", "makedb", "--in", str(reference), "--db", str(db)],
                       check=True, capture_output=True)
        for limit in ("25", "all"):
            out = self.base / ("real_" + limit)
            command = [
                sys.executable, str(ROOT / "spraynpray_v2.py"), "-i", str(contigs),
                "-d", str(db), "--protein-to-genome", str(mapping), "--assembly-info", str(metadata),
                "--taxdump", str(self.dump), "--meta", "--threads", "2", "--orfs-per-contig", limit,
                "-o", str(out)]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            with (out / "contigs.csv").open() as f:
                rows = list(csv.DictReader(f))
            classified = [r for r in rows if r["status"] == "classified"]
            self.assertTrue(classified)
            self.assertTrue(all(r["predicted_taxid"] == "11" for r in classified))
            if limit == "25":
                self.assertTrue(all(int(r["selected_orfs"]) <= 25 for r in rows))
            else:
                self.assertTrue(all(r["total_orfs"] == r["selected_orfs"] for r in rows))
            with (out / "top_hits.csv").open() as f:
                top = list(csv.DictReader(f))
            self.assertEqual(len(top), sum(int(r["selected_orfs"]) for r in rows))
            self.assertTrue(all(r["source_genome_accession"] == "GCF_000001.1"
                                for r in top if r["status"] != "no_hit"))
            before = {name: (out / name).stat().st_mtime_ns for name in ("prodigal.faa", "diamond.tsv")}
            repeated = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(repeated.returncode, 0, repeated.stderr)
            self.assertEqual(before, {name: (out / name).stat().st_mtime_ns for name in before})
            manifest = json.loads((out / "run_manifest.json").read_text())
            self.assertEqual(manifest["stages"], {"prodigal": "reused", "diamond": "reused"})
            self.assertEqual(manifest["commands"], [])


if __name__ == "__main__":
    unittest.main()
