"""Cohort summary, ANI evidence, missing-data and distance-map regression tests."""

import csv
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("spraynpray_compare", ROOT / "spraynpray_compare.py")
compare = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(compare)


class CompareTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.results = self.base / "results"
        self.results.mkdir()
        self.out = self.base / "out"
        self.out.mkdir()
        self.work = self.out / "work"
        self.work.mkdir()
        self.args = SimpleNamespace(outdir=self.out, min_af=.5, min_ani=80, neighbors=2)

    def summary(self, sample="sample_a", text=None, count=4):
        path = self.results / (sample + ".csv")
        if text is None:
            text = ("c_1: Species, strain X [WP_1.1; GCF_000001.1; 100.00% identity]; "
                    "c_2: Species, strain X [WP_2.1; GCF_000001.1; 98.00% identity]; "
                    "c_3: Species Y [WP_3.1; GCA_000002.2; 95.00% identity]; c_4: no hit")
        with path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["contig", "selected_orfs", "summary_orfs_omitted", "closest_blast_hits"])
            writer.writeheader()
            writer.writerow(dict(contig="c", selected_orfs=count, summary_orfs_omitted=0, closest_blast_hits=text))
        return path

    def nodes(self):
        return [dict(id="a", kind="sample", staged_file="genomes/a.fa", source_file="/input/a.fa"),
                dict(id="b", kind="sample", staged_file="genomes/b.fa", source_file="/input/b.fa"),
                dict(id="GCF_000001.1", kind="reference", staged_file="genomes/r.fa", source_file="/refs/r.fa")]

    def measurement(self, ani, af=.9):
        return dict(ani=ani, mapped_fragments=int(af * 100), query_fragments=100, aligned_fraction=af)

    def test_summary_counts_percentages_and_mean_are_per_orf(self):
        master, overview, _ = compare.load_sample(self.summary())
        self.assertEqual(overview["selected_orfs"], 4)
        self.assertEqual(overview["matched_orfs"], 3)
        self.assertEqual(overview["no_hit_orfs"], 1)
        self.assertEqual(master[0]["orf_hits"], 2)
        self.assertAlmostEqual(master[0]["hit_percent"], 200 / 3)
        self.assertEqual(master[0]["selected_orf_percent"], 50)
        self.assertEqual(master[0]["mean_protein_identity"], 99)
        self.assertAlmostEqual(sum(r["hit_percent"] for r in master), 100)
        self.assertEqual(master[0]["organism_names"], "Species, strain X")

    def test_summary_folder_and_reference_subset_do_not_change_tallies(self):
        self.summary()
        self.summary("sample_b")
        master, overview, refs = compare.summarize(self.results, self.out, top_references=1)
        self.assertEqual(len(master), 4)
        self.assertEqual(len(overview), 2)
        self.assertEqual(refs, ["GCF_000001.1"])
        with (self.out / "master_summary.csv").open() as f:
            self.assertEqual(len(list(csv.DictReader(f))), 4)

    def test_top_hits_csv_unmapped_and_no_hit(self):
        path = self.results / "sample_a.csv"
        with path.open("w", newline="") as f:
            fields = ["orf", "contig", "protein_accession", "source_genome_accession", "pident", "taxa", "status"]
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerow(dict(orf="a_1", contig="a", protein_accession="WP_1.1", source_genome_accession="", pident=99, taxa="", status="missing_protein_mapping"))
            writer.writerow(dict(orf="a_2", contig="a", protein_accession="", source_genome_accession="", pident="", taxa="no hit", status="no_hit"))
        master, overview, _ = compare.load_sample(path)
        self.assertEqual(master[0]["reference_assembly"], "unmapped")
        self.assertEqual(master[0]["hit_percent"], 100)
        self.assertEqual(overview["unmapped_hit_orfs"], 1)

    def test_truncated_malformed_duplicate_and_invalid_identity_rejected(self):
        path = self.summary()
        path.write_text(path.read_text().replace(",4,0,", ",4,2,"))
        with self.assertRaisesRegex(ValueError, "Truncated"):
            compare.load_sample(path)
        with self.assertRaisesRegex(ValueError, "Unsupported"):
            list(compare.parse_summary_entries("a_1: Species [GCF_000001.1] | extra"))
        with self.assertRaisesRegex(ValueError, "between 0 and 100"):
            list(compare.parse_summary_entries("a_1: Species [WP_1.1; GCF_000001.1; 101.00% identity]"))
        self.summary(text="c_1: no hit; c_1: no hit", count=2)
        with self.assertRaisesRegex(ValueError, "Duplicate ORF"):
            compare.load_sample(path)

    def test_fastani_both_directions_coverage_and_missing_pairs(self):
        nodes = self.nodes()
        directional = {("a", "b"): self.measurement(99), ("b", "a"): self.measurement(98),
                       ("a", "GCF_000001.1"): self.measurement(100, .1),
                       ("GCF_000001.1", "a"): self.measurement(100, .9)}
        edges = compare.pairwise_tables(nodes, directional, self.args)
        self.assertEqual(len(edges), 1)
        self.assertEqual(edges[0]["ani_used"], 98.5)
        self.assertEqual(edges[0]["distance"], 1.5)
        with (self.out / "ani_pairs.csv").open() as f:
            pairs = list(csv.DictReader(f))
        self.assertEqual([r["status"] for r in pairs], ["accepted", "low_aligned_fraction", "not_reported"])
        with (self.out / "ani_matrix.csv").open() as f:
            matrix = list(csv.DictReader(f))
        self.assertEqual(matrix[0]["GCF_000001.1"], "")
        self.assertEqual(matrix[1]["a"], "98.5")
        directional[("b", "GCF_000001.1")] = self.measurement(99)
        compare.pairwise_tables(nodes, directional, self.args)
        with (self.out / "ani_pairs.csv").open() as f:
            self.assertEqual(list(csv.DictReader(f))[2]["status"], "missing_direction")

    def test_fastani_parser_aliases_duplicate_validation(self):
        raw = self.base / "ani.tsv"
        raw.write_text("genomes/a.fa\tgenomes/b.fa\t99.5\t90\t100\n")
        values = compare.parse_fastani(raw, self.nodes(), self.work)
        self.assertEqual(values[("a", "b")]["aligned_fraction"], .9)
        raw.write_text("genomes/a.fa\tgenomes/b.fa\t99.5\t101\t100\n")
        with self.assertRaisesRegex(ValueError, "fragment"):
            compare.parse_fastani(raw, self.nodes(), self.work)
        raw.write_text("unknown\tb\t99\t90\t100\n")
        with self.assertRaisesRegex(ValueError, "Unknown"):
            compare.parse_fastani(raw, self.nodes(), self.work)

    def test_fasta_staging_preserves_sequences_and_n50(self):
        source = self.base / "a.fa"
        source.write_text(">a\n" + "ACGT" * 100 + "\n>b\n" + "N" * 100 + "\n")
        target = self.work / "a.fa"
        stats = compare.stage_fasta(source, target)
        self.assertEqual((stats["length_bp"], stats["contigs"], stats["n50"]), (500, 2, 400))
        before = target.stat().st_mtime_ns
        compare.stage_fasta(source, target)
        self.assertEqual(before, target.stat().st_mtime_ns)
        source.write_text(">a\nMPEPTIDE\n")
        with self.assertRaisesRegex(ValueError, "Non-DNA"):
            compare.stage_fasta(source, target)

    def test_reference_zip_versions_and_path_safety(self):
        archive = self.base / "reference.zip"
        reference = "GCF_000001.1"
        with zipfile.ZipFile(archive, "w") as z:
            z.writestr("ncbi_dataset/data/%s/%s_genomic.fna" % (reference, reference), ">a\nACGT\n")
            z.writestr("ncbi_dataset/data/GCF_000001.2/GCF_000001.2_genomic.fna", ">wrong_version\nACGT\n")
            z.writestr("../../escape.fa", ">bad\nACGT\n")
        result = compare.extract_reference_package(archive, self.work, {reference})
        self.assertEqual(list(result), [reference])
        self.assertFalse((self.base.parent / "escape.fa").exists())

    def test_ani_cache_reuse_and_invalidation(self):
        nodes = self.nodes()
        for i, node in enumerate(nodes):
            node["sha256"] = str(i)
        args = SimpleNamespace(frag_len=3000, fastani="fastANI", force_ani=False, fastani_results=None, threads=2)
        count = []

        def fake_command(command, log, cwd=None):
            count.append(command)
            (cwd / "fastani.tsv.partial").write_text("genomes/a.fa\tgenomes/b.fa\t99\t90\t100\n")

        with patch.object(compare, "run_command", side_effect=fake_command), patch.object(shutil, "which", return_value="fastANI"):
            self.assertEqual(compare.calculate_ani(nodes, args, self.work)[1], "ran")
            args.threads = 4
            self.assertEqual(compare.calculate_ani(nodes, args, self.work)[1], "reused")
            self.assertEqual(len(count), 1)
            nodes[0]["sha256"] = "changed"
            self.assertEqual(compare.calculate_ani(nodes, args, self.work)[1], "ran")
            (self.work / "fastani.tsv").write_text("damaged")
            self.assertEqual(compare.calculate_ani(nodes, args, self.work)[1], "ran")
            self.assertEqual(len(count), 3)

    def test_reference_download_command_and_reuse(self):
        args = SimpleNamespace(reference_dir=None, no_download=False, datasets="datasets",
                               allow_missing_references=False, outdir=self.out)
        calls = []

        def fake_download(command, log, cwd=None):
            calls.append(command)
            requested = Path(command[command.index("--inputfile") + 1]).read_text().splitlines()
            archive = Path(command[command.index("--filename") + 1])
            with zipfile.ZipFile(archive, "w") as z:
                for accession in requested:
                    z.writestr("ncbi_dataset/data/%s/%s_genomic.fna" % (accession, accession), ">reference\nACGT\n")

        with patch.object(compare, "run_command", side_effect=fake_download), patch.object(shutil, "which", return_value="datasets"):
            resolved = compare.resolve_references(["GCF_000001.1", "GCA_000002.2"], args, self.work)
            self.assertTrue(all(p.is_file() for p in resolved.values()))
            compare.resolve_references(["GCF_000001.1", "GCA_000002.2"], args, self.work)
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0][:4], ["datasets", "download", "genome", "accession"])
        self.assertEqual(calls[0][calls[0].index("--include") + 1], "genome")
        args.no_download = True
        args.allow_missing_references = True
        self.assertIsNone(compare.resolve_references(["GCF_000001.2"], args, self.work)["GCF_000001.2"])
        self.assertTrue((self.out / "missing_references.csv").is_file())

    def test_failed_fastani_does_not_create_a_completed_cache(self):
        nodes = self.nodes()
        for node in nodes:
            node["sha256"] = node["id"]
        args = SimpleNamespace(frag_len=3000, fastani="fastANI", force_ani=False, fastani_results=None, threads=2)

        def fail(command, log, cwd=None):
            (cwd / "fastani.tsv.partial").write_text("unfinished")
            raise subprocess.CalledProcessError(1, command)

        with patch.object(compare, "run_command", side_effect=fail), patch.object(shutil, "which", return_value="fastANI"):
            with self.assertRaises(subprocess.CalledProcessError):
                compare.calculate_ani(nodes, args, self.work)
        self.assertFalse((self.work / "ani_cache.json").exists())
        self.assertFalse((self.work / "fastani.tsv").exists())

    def test_embedding_uses_true_distances_and_keeps_components_separate(self):
        nodes = [dict(id=x, kind="sample") for x in ("a", "b", "c", "unlinked")]
        # Exact triangle with side lengths 1, 2, 2.5 in ANI percentage points.
        edges = [dict(genome_a="a", genome_b="b", distance=1), dict(genome_a="b", genome_b="c", distance=2),
                 dict(genome_a="a", genome_b="c", distance=2.5)]
        coords, diagnostics = compare.embed_network(nodes, edges)
        for edge in edges:
            a, b = coords[edge["genome_a"]], coords[edge["genome_b"]]
            self.assertAlmostEqual(math.hypot(a["x"] - b["x"], a["y"] - b["y"]), edge["distance"], places=4)
        self.assertNotEqual(coords["a"]["component"], coords["unlinked"]["component"])
        self.assertLess(diagnostics[0]["normalized_stress"], 1e-4)
        self.assertEqual(len(compare.display_edges(nodes, edges, 1)), 2)

    def test_nearest_reference_uses_ani_handles_ties_and_missing_matches(self):
        nodes = [dict(id=x, kind="sample") for x in ("a", "b")] + [dict(id=x, kind="reference") for x in ("r1", "r2", "r3")]
        pairs = [dict(genome_a="a", genome_b="r1", ani_used=99, min_aligned_fraction=.9, status="accepted"),
                 dict(genome_a="r2", genome_b="a", ani_used=100, min_aligned_fraction=.8, status="accepted"),
                 dict(genome_a="a", genome_b="r3", ani_used=100, min_aligned_fraction=.9, status="accepted"),
                 dict(genome_a="b", genome_b="r1", ani_used="", min_aligned_fraction=.1, status="low_aligned_fraction")]
        nearest = compare.nearest_reference_summary(nodes, pairs, self.out)
        self.assertEqual(nearest[0]["nearest_reference"], "r3")
        self.assertEqual(nearest[0]["estimated_ani"], 100)
        self.assertEqual(nearest[0]["ani_distance"], 0)
        self.assertEqual(nearest[0]["tied_reference_count"], 2)
        self.assertEqual(nearest[1]["estimated_ani"], "")
        self.assertEqual(nearest[1]["status"], "no_accepted_reference")

    def test_top_hit_metadata_survives_summary_and_reference_nodes(self):
        path = self.results / "sample_a.csv"
        fields = ["orf", "protein_accession", "source_genome_accession", "pident", "taxa", "status",
                  "assembly_infraspecific_name", "assembly_isolate", "assembly_bioproject", "assembly_taxid",
                  "assembly_assembly_accession", "metadata_match"]
        compare.write_csv(path, fields, [dict(orf="a_1", protein_accession="WP_1.1", source_genome_accession="GCF_000001.1",
            pident=99, taxa="Pseudomonas aeruginosa", status="matched", assembly_infraspecific_name="strain=PAO1",
            assembly_isolate="isolate-X", assembly_bioproject="PRJNA1", assembly_taxid="287",
            assembly_assembly_accession="GCA_000001.1", metadata_match="identical_pair")])
        master, overview, refs = compare.summarize(self.results, self.out)
        self.assertEqual(master[0]["strain"], "PAO1")
        self.assertEqual(master[0]["isolate"], "isolate-X")
        fasta = self.base / "sample_a.fa"
        fasta.write_text(">a\nACGT\n")
        args = SimpleNamespace(genomes=self.base, genome_extension=".fa")
        with patch.object(compare, "resolve_references", return_value={"GCF_000001.1": fasta}):
            nodes = compare.prepare_nodes(overview, refs, master, args, self.work)
        self.assertEqual(nodes[1]["strain"], "PAO1")
        self.assertEqual(nodes[1]["metadata_match"], "identical_pair")
        self.assertNotIn("strain", nodes[0])

    def test_metadata_enrichment_uses_exact_versions_and_missing_values(self):
        path = self.base / "metadata.tsv.gz"
        import gzip
        with gzip.open(path, "wt") as f:
            f.write("## assembly summary\n#assembly_accession\torganism_name\tinfraspecific_name\tisolate\ttaxid\n"
                    "GCF_000001.1\tPseudomonas aeruginosa\tstrain=PAO1\tna\t287\n"
                    "GCA_000002.1\tOther species\tstrain=wrong\twrong\t1\n")
        nodes = [dict(id="GCF_000001.1", kind="reference"), dict(id="GCF_000002.1", kind="reference")]
        compare.enrich_reference_nodes(nodes, path)
        self.assertEqual(nodes[0]["strain"], "PAO1")
        self.assertEqual(nodes[0].get("isolate", ""), "")
        self.assertEqual(nodes[0]["taxid"], "287")
        self.assertNotIn("strain", nodes[1])
        self.assertEqual(compare.annotation_values({"infraspecific_name": "serovar=Typhi"})["strain"], "")

    def test_click_details_separate_sample_and_reference_strain(self):
        nodes = [dict(id="sample", kind="sample", length_bp=1000, contigs=2, n50=600),
                 dict(id="ref", kind="reference", organism="Pseudomonas aeruginosa", strain="PAO1", isolate="X")]
        edges = [dict(genome_a="sample", genome_b="ref", ani_used=99.5, min_aligned_fraction=.9, status="accepted")]
        sample = compare.genome_details(nodes[0], nodes, edges)
        self.assertIn("Strain: Not provided", sample)
        self.assertIn("Reference strain: PAO1", sample)
        self.assertIn("Mean reciprocal genome ANI: 99.50000%", sample)
        self.assertIn("Length (bp): 1000", sample)
        reference = compare.genome_details(nodes[1], nodes, edges)
        self.assertIn("Strain: PAO1", reference)
        self.assertIn("Isolate: X", reference)
        self.assertIn("Closest input sample: sample", reference)

    def test_network_caps_reference_edges_per_sample_by_ani_and_keeps_samples(self):
        samples = [dict(id=x, kind="sample") for x in ("a", "b", "unlinked")]
        refs = [dict(id="r%02d" % i, kind="reference", orf_hits=1000-i) for i in range(15)]
        nodes = samples + refs
        edges = []
        for i, ref in enumerate(refs[:12]):
            for sample, ani in (("a", 99 - i / 10), ("b", 96 + i / 10)):
                edges.append(dict(genome_a=ref["id"], genome_b=sample, ani_used=ani,
                                  min_aligned_fraction=.9, distance=100-ani, status="accepted"))
        sample_pair = dict(genome_a="a", genome_b="b", distance=1, status="accepted")
        edges.append(sample_pair)
        # A close reference/reference edge must not pull in an unrelated node.
        edges.append(dict(genome_a="r00", genome_b="r14", distance=.001, status="accepted"))
        visible, shown = compare.select_network_references(nodes, list(reversed(edges)))
        links = {sample: {e["genome_a"] for e in shown if e["genome_b"] == sample and e["genome_a"].startswith("r")} for sample in ("a", "b")}
        self.assertEqual(links["a"], {"r%02d" % i for i in range(5)})
        self.assertEqual(links["b"], {"r%02d" % i for i in range(7, 12)})
        self.assertEqual(len(shown), 11)
        self.assertIn(sample_pair, shown)
        self.assertEqual({n["id"] for n in visible}, {n["id"] for n in samples + refs[:5] + refs[7:12]})
        self.assertEqual(len({n["id"] for n in visible}), len(visible))

    def test_network_reference_ties_zero_and_unmatched_samples(self):
        nodes = [dict(id="a", kind="sample")] + [dict(id=x, kind="reference") for x in ("r1", "r2", "r3", "unused")]
        edges = [dict(genome_a="a", genome_b=ref, ani_used=99, min_aligned_fraction=af,
                      distance=1, status="accepted") for ref, af in (("r3", .9), ("r1", .8), ("r2", .9))]
        excluded = dict(genome_a="a", genome_b="unused", ani_used=100,
                        min_aligned_fraction=.1, distance=.001, status="low_aligned_fraction")
        visible, shown = compare.select_network_references(nodes, edges + [excluded], limit=1)
        self.assertEqual([n["id"] for n in visible], ["a", "r2"])
        self.assertEqual(shown[0]["genome_b"], "r2")
        self.assertEqual(compare.select_network_references(nodes, edges, limit=0, neighbors=0), (nodes, edges))
        self.assertEqual(compare.select_network_references(nodes, []), (nodes[:1], []))
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            compare.select_network_references(nodes, edges, limit=-1)
        args = compare.parser().parse_args(["--results", "results", "--outdir", "out"])
        self.assertEqual(args.network_references_per_sample, 5)
        self.assertEqual(args.top_references_per_sample, 0)

    def test_helper_filters_existing_network_and_can_restore_all_without_ani(self):
        nodes = [dict(id="sample", kind="sample")] + [dict(id="r%02d" % i, kind="reference") for i in range(13)]
        for node in nodes:
            node.update(organism="Test species", source_file="local.fa", staged_file="local.fa", sequence_status="available")
        directional = {}
        for i, ref in enumerate(nodes[1:]):
            for pair in (("sample", ref["id"]), (ref["id"], "sample")):
                directional[pair] = self.measurement(99-i/10)
        edges = compare.pairwise_tables(nodes, directional, self.args)
        positions = {node["id"]: dict(x=i/10, y=i%3, component=1) for i, node in enumerate(nodes)}
        diagnostics = [dict(component=1, nodes=len(nodes), measured_pairs=len(edges), normalized_stress=0, max_distance_error=0)]
        self.args.network_references_per_sample = 0
        compare.draw_network(nodes, edges, positions, diagnostics, self.args)
        pairs_before = compare.sha256(self.out / "ani_pairs.csv")
        metadata = self.base / "assembly.tsv"
        metadata.write_text("#assembly_accession\torganism_name\tinfraspecific_name\tisolate\n"
                            "r00\tTest species\tstrain=PAO1 <special>\tiso-X\n")
        # Default helper selection changes plots only, not ANI evidence or saved coordinates.
        with patch.object(compare, "run_command", side_effect=AssertionError("no external tools")), \
             patch.object(compare, "embed_network", side_effect=AssertionError("no layout recomputation")):
            result = compare.refresh_existing_analysis(self.out, assembly_info=metadata)
            self.assertEqual(result[0]["nearest_reference"], "r00")
            _, saved = compare.read_csv(self.out / "network_nodes.csv")
            self.assertEqual(len(saved), 14)
            self.assertEqual(sum(n["displayed"] == "True" for n in saved), 6)
            for n in saved:
                self.assertEqual(float(n["x"]), positions[n["id"]]["x"])
                self.assertEqual(float(n["y"]), positions[n["id"]]["y"])
            _, shown = compare.read_csv(self.out / "network_edges.csv")
            self.assertEqual(len(shown), 5)
            page = (self.out / "ani_network.html").read_text()
            self.assertEqual(page.count('<circle tabindex='), 6)
            self.assertNotIn('data-id="r12"', page)
            self.assertIn('up to 5 reference connections', page)
            self.assertIn('Strain: PAO1 &lt;special&gt;', page)
            self.assertIn('aria-label="Genome details"', page)
            self.assertIn("e.target.closest('circle')", page)
            self.assertIn("c.addEventListener('click',choose)", page)
            self.assertIn("c.addEventListener('keydown'", page)
            compare.refresh_existing_analysis(self.out, network_references_per_sample=0)
            self.assertEqual((self.out / "ani_network.html").read_text().count('<circle tabindex='), 14)
        self.assertEqual(pairs_before, compare.sha256(self.out / "ani_pairs.csv"))

    @unittest.skipUnless(shutil.which("node"), "requires Node.js for generated interaction checks")
    def test_generated_html_click_keyboard_selection_and_circle_pointerdown(self):
        import re
        import xml.etree.ElementTree as ET
        nodes = self.nodes()
        for node in nodes:
            node.update(organism="Pseudomonas aeruginosa", sequence_status="available")
        nodes[2]["strain"] = "PAO1"
        directional = {(a["id"], b["id"]): self.measurement(99.5) for a in nodes for b in nodes if a != b}
        edges = compare.pairwise_tables(nodes, directional, self.args)
        positions = {n["id"]: dict(x=i, y=i%2, component=1) for i, n in enumerate(nodes)}
        compare.draw_network(nodes, edges, positions, [dict(component=1, nodes=3)], self.args)
        page = (self.out / "ani_network.html").read_text()
        svg = ET.fromstring(re.search(r'(<svg.*?</svg>)', page, re.S).group(1))
        circles = [dict(id=e.attrib["data-id"], info=e.attrib["data-info"]) for e in svg.iter()
                   if e.tag.endswith("}circle")]
        script = re.search(r'<script>(.*?)</script>', page, re.S).group(1)
        harness = r"""
const vm=require('node:vm'),assert=require('node:assert/strict');
const element=(dataset={})=>({dataset,events:{},attrs:{},textContent:'',value:'',
  addEventListener(type,callback){this.events[type]=callback},setAttribute(k,v){this.attrs[k]=v}});
const circles=INPUT_CIRCLES.map(c=>element(c)), svg=element(),elements={};
for(const id of ['drawing','find','info','reset','samples','references'])elements[id]=element();
let captures=0;svg.setPointerCapture=()=>captures++;
const document={getElementById:id=>elements[id],querySelector:()=>svg,
  querySelectorAll:selector=>selector==='circle'?circles:[]};
vm.runInNewContext(INPUT_SCRIPT,{document,DOMPoint:class{}});
svg.events.pointerdown({target:{closest:()=>circles[0]}});
assert.equal(captures,0,'circle selection must not be captured by background panning');
circles[0].events.click();assert.equal(elements.info.textContent,circles[0].dataset.info);
assert.equal(elements.find.value,circles[0].dataset.id);assert.equal(circles[0].attrs.stroke,'#197a9a');
let prevented=false;const reference=circles.find(c=>c.dataset.id.startsWith('GCF_'));
reference.events.keydown({key:'Enter',preventDefault:()=>prevented=true});
assert.equal(prevented,true);assert.ok(elements.info.textContent.includes('Strain: PAO1'));
assert.ok(elements.info.textContent.includes('Mean reciprocal genome ANI: 99.50000%'));
elements.find.onchange({target:{value:circles[1].dataset.id}});
assert.equal(elements.info.textContent,circles[1].dataset.info);
elements.find.onchange({target:{value:''}});assert.ok(elements.info.textContent.startsWith('Select a genome'));
"""
        # Substitute code after data so metadata cannot alter the JavaScript harness.
        harness = harness.replace("INPUT_CIRCLES", json.dumps(circles)).replace("INPUT_SCRIPT", json.dumps(script))
        subprocess.run([shutil.which("node"), "-"], input=harness, text=True, capture_output=True, check=True)

    def test_reference_label_boxes_do_not_overlap_or_move_genome_points(self):
        nodes = [dict(id="r%03d" % i, kind="reference", organism="Pseudomonas aeruginosa very long strain name " * 2) for i in range(60)]
        nodes += [dict(id="sample", kind="sample")]
        points = {n["id"]: dict(x=1, y=1) for n in nodes}
        points["sample"] = dict(x=2, y=1)
        original = json.dumps(points, sort_keys=True)
        xy, labels, width, height = compare.expanded_reference_layout(nodes, points)
        self.assertEqual(original, json.dumps(points, sort_keys=True))
        self.assertEqual(len(labels), 60)
        self.assertGreater(width, 1200)
        for side in (0, 1):
            column = sorted((l for l in labels if l["side"] == side), key=lambda l: l["y"])
            for a, b in zip(column, column[1:]):
                self.assertGreaterEqual(b["y"], a["y"] + a["height"] + 7.99)
            self.assertGreaterEqual(column[0]["y"], 39.99)
            self.assertLessEqual(column[-1]["y"] + column[-1]["height"], height - 39.99)
        self.assertEqual(xy["r000"], xy["r059"])
        self.assertGreater(xy["sample"][0], xy["r000"][0])

    def test_sample_scope_skips_reference_only_pairs_and_reuses_full_cache(self):
        nodes = self.nodes() + [dict(id="r2", kind="reference", staged_file="genomes/r2.fa", source_file="/refs/r2.fa")]
        for node in nodes:
            node["sha256"] = node["id"]
        args = SimpleNamespace(frag_len=3000, fastani="fastANI", force_ani=False, fastani_results=None, threads=2, ani_scope="samples")
        calls = []

        def fake(command, log, cwd=None):
            queries = (cwd / command[command.index("--ql") + 1]).read_text().splitlines()
            targets = (cwd / command[command.index("--rl") + 1]).read_text().splitlines()
            calls.append((queries, targets))
            (cwd / command[command.index("-o") + 1]).write_text("".join(
                "%s\t%s\t99\t90\t100\n" % (q, r) for q in queries for r in targets))

        with patch.object(compare, "run_command", side_effect=fake), patch.object(shutil, "which", return_value="fastANI"):
            values, stage = compare.calculate_ani(nodes, args, self.work)
            self.assertEqual(stage, "ran")
            self.assertEqual(len(calls), 2)
            self.assertEqual(sum(len(q) * len(r) for q, r in calls), 12)  # 2^2 + 2*2*2; full union = 16
            self.assertNotIn(("GCF_000001.1", "r2"), values)
            self.assertIn(("a", "r2"), values)
            self.assertIn(("r2", "a"), values)
            self.assertEqual(compare.calculate_ani(nodes, args, self.work)[1], "reused")
            self.assertEqual(len(calls), 2)
            (self.work / "fastani.tsv").unlink()
            compare.calculate_ani(nodes, args, self.work)  # completed directional jobs reassemble without tool calls
            self.assertEqual(len(calls), 2)
            args.ani_scope = "all"
            compare.calculate_ani(nodes, args, self.work)
            self.assertEqual(len(calls), 3)
            args.ani_scope = "samples"
            self.assertEqual(compare.calculate_ani(nodes, args, self.work)[1], "reused")
            self.assertEqual(len(calls), 3)
        graph_args = SimpleNamespace(outdir=self.out, min_af=.5, min_ani=80, ani_scope="samples")
        compare.pairwise_tables(nodes, values, graph_args)
        with (self.out / "ani_pairs.csv").open() as f:
            rows = list(csv.DictReader(f))
        self.assertFalse(any(r["kind_a"] == r["kind_b"] == "reference" for r in rows))

    def test_old_all_scope_cache_is_reused_by_new_default(self):
        nodes = self.nodes()
        for node in nodes:
            node["sha256"] = node["id"]
        signature = dict(genomes=[dict(id=n["id"], sha256=n["sha256"]) for n in nodes],
                         fragment_length=3000, min_fraction=.2, executable="fastANI")
        raw = self.work / "fastani.tsv"
        raw.write_text("genomes/a.fa\tgenomes/b.fa\t99\t90\t100\n")
        compare.write_json(self.work / "ani_cache.json", dict(signature=signature, output_sha256=compare.sha256(raw)))
        args = SimpleNamespace(frag_len=3000, fastani="fastANI", force_ani=False, fastani_results=None, threads=2, ani_scope="samples")
        with patch.object(compare, "run_command", side_effect=AssertionError("must reuse")):
            self.assertEqual(compare.calculate_ani(nodes, args, self.work)[1], "reused")

    @unittest.skipUnless(shutil.which("fastANI"), "requires real FastANI")
    def test_real_fastani_and_offline_network(self):
        # Large random bacterial-like DNA and a known closely related derivative.
        import random
        rng = random.Random(9)
        sequence = "".join(rng.choices("ACGT", k=250000))
        mutant = list(sequence)
        for index in range(0, len(mutant), 100):
            mutant[index] = "A" if mutant[index] != "A" else "C"
        genomes = self.base / "genomes"
        references = self.base / "references"
        genomes.mkdir()
        references.mkdir()
        (genomes / "sample_a.fa").write_text(">a\n" + sequence + "\n")
        (genomes / "sample_b.fa").write_text(">b\n" + "".join(mutant) + "\n")
        (references / "GCF_000001.1.fna").write_text(">reference\n" + sequence + "\n")
        text = "c_1: Test organism [WP_1.1; GCF_000001.1; 100.00% identity]"
        self.summary("sample_a", text, 1)
        self.summary("sample_b", text, 1)
        command = [sys.executable, str(ROOT / "spraynpray_compare.py"), "--results", str(self.results),
                   "--genomes", str(genomes), "--reference-dir", str(references), "--no-download",
                   "--outdir", str(self.out), "--threads", "2"]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.out / "ani_pairs.csv").open() as f:
            pairs = list(csv.DictReader(f))
        self.assertEqual(len(pairs), 3)
        self.assertTrue(all(r["status"] == "accepted" for r in pairs))
        self.assertTrue(all(98 < float(r["ani_used"]) <= 100 for r in pairs))
        for name in ("ani_network.png", "ani_network.svg", "ani_network.html"):
            self.assertGreater((self.out / name).stat().st_size, 1000)
        self.assertNotIn("https://cdn", (self.out / "ani_network.html").read_text())
        before = (self.work / "fastani.tsv").stat().st_mtime_ns
        again = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(again.returncode, 0, again.stderr)
        self.assertEqual(before, (self.work / "fastani.tsv").stat().st_mtime_ns)
        self.assertEqual(json.loads((self.out / "run_manifest.json").read_text())["ani_stage"], "reused")
        pairs_before = compare.sha256(self.out / "ani_pairs.csv")
        with (self.out / "network_nodes.csv").open() as f:
            positions_before = [(r["id"], float(r["x"]), float(r["y"])) for r in csv.DictReader(f)]
        helper = subprocess.run([sys.executable, str(ROOT / "spraynpray_helper.py"), "--analysis-dir", str(self.out)], capture_output=True, text=True)
        self.assertEqual(helper.returncode, 0, helper.stderr)
        self.assertEqual(pairs_before, compare.sha256(self.out / "ani_pairs.csv"))
        self.assertEqual(before, (self.work / "fastani.tsv").stat().st_mtime_ns)
        with (self.out / "network_nodes.csv").open() as f:
            positions_after = [(r["id"], float(r["x"]), float(r["y"])) for r in csv.DictReader(f)]
        self.assertEqual(positions_before, positions_after)
        with (self.out / "nearest_references.csv").open() as f:
            nearest = list(csv.DictReader(f))
        self.assertEqual(len(nearest), 2)
        self.assertTrue(all(r["nearest_reference"] == "GCF_000001.1" for r in nearest))
        page = (self.out / "ani_network.html").read_text()
        self.assertIn('id="references" type="checkbox" checked', page)
        self.assertIn('stroke-dasharray="3 3"', page)


if __name__ == "__main__":
    unittest.main()
