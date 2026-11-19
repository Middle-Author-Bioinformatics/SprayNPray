"""Wrapper discovery, cohort-input selection, and real end-to-end/resume tests."""

import csv
import gzip
import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import spraynpray_pipeline as pipeline
import spraynpray_compare as compare


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)

    def test_discovery_handles_compression_and_rejects_colliding_ids(self):
        folder = self.base / "genomes"
        folder.mkdir()
        (folder / "sample_a.fa").write_text(">a\nACGT\n")
        with gzip.open(folder / "sample_b.fna.gz", "wt") as f:
            f.write(">b\nACGT\n")
        (folder / "readme.txt").write_text("ignored")
        self.assertEqual(list(pipeline.discover_genomes(folder)), ["sample_a", "sample_b"])
        (folder / "sample_a.fasta").write_text(">a\nACGT\n")
        with self.assertRaisesRegex(ValueError, "duplicate"):
            pipeline.discover_genomes(folder)

    def test_wrapper_passes_compare_scope_and_reference_controls(self):
        args = pipeline.parser().parse_args([
            "--genomes", str(self.base), "--diamond-db", "db", "--protein-to-genome", "map.tsv",
            "--assembly-info", "metadata.tsv", "--taxdump", "taxdump", "--outdir", str(self.base / "output"),
            "--threads", "16", "--no-download", "--reference-dir", str(self.base), "--top-references-per-sample", "3"])
        command = pipeline.compare_command(args)
        self.assertEqual(command[command.index("--ani-scope") + 1], "samples")
        self.assertEqual(command[command.index("--threads") + 1], "16")
        self.assertEqual(command[command.index("--top-references-per-sample") + 1], "3")
        self.assertIn("--no-download", command)
        self.assertEqual(command[command.index("--assembly-info") + 1], str(Path("metadata.tsv").resolve()))
        self.assertEqual(command[command.index("--network-references-per-sample") + 1], "5")
        args.network_references_per_sample = 0
        command = pipeline.compare_command(args)
        self.assertEqual(command[command.index("--network-references-per-sample") + 1], "0")

    @unittest.skipUnless(all(shutil.which(t) for t in ("prodigal", "diamond", "fastANI")), "requires all three real tools")
    def test_real_pipeline_two_samples_complete_and_resume(self):
        genomes = self.base / "input_genomes"
        references = self.base / "references"
        genomes.mkdir()
        references.mkdir()
        toy = ROOT / "toy_dataset" / "test_dataset.fasta"
        shutil.copyfile(toy, genomes / "sample_a.fa")
        with open(toy, "rb") as source, gzip.open(genomes / "sample_b.fna.gz", "wb") as out:
            shutil.copyfileobj(source, out)
        reference = "GCF_999999999.1"
        shutil.copyfile(toy, references / (reference + ".fna"))
        proteins = self.base / "prodigal.faa"
        subprocess.run(["prodigal", "-i", str(toy), "-a", str(proteins), "-o", str(self.base / "genes.gff"),
                        "-f", "gff", "-p", "meta"], check=True, capture_output=True)
        fasta = self.base / "database.faa"
        # Protein data are not DNA: use a separate simple FASTA reader here.
        sequences, chunks = [], []
        for line in proteins.read_text().splitlines():
            if line.startswith(">"):
                if chunks:
                    sequences.append("".join(chunks))
                chunks = []
            else:
                chunks.append(line)
        if chunks:
            sequences.append("".join(chunks))
        fasta.write_text("".join(">WP_%09d.1\n%s\n" % (i, seq) for i, seq in enumerate(sequences, 1)))
        mapping = self.base / "map.tsv"
        mapping.write_text("".join("WP_%09d.1\t%s\n" % (i, reference) for i in range(1, len(sequences) + 1)))
        metadata = self.base / "assembly.tsv"
        metadata.write_text("#assembly_accession\ttaxid\torganism_name\n%s\t11\tTest species\n" % reference)
        dump = self.base / "taxdump"
        dump.mkdir()
        (dump / "nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tdomain\t|\n11\t|\t2\t|\tspecies\t|\n")
        (dump / "names.dmp").write_text("1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tBacteria\t|\t\t|\tscientific name\t|\n11\t|\tTest species\t|\t\t|\tscientific name\t|\n")
        db = self.base / "database"
        subprocess.run(["diamond", "makedb", "--in", str(fasta), "--db", str(db)], check=True, capture_output=True)
        out = self.base / "pipeline_output"
        command = [sys.executable, str(ROOT / "spraynpray_pipeline.py"), "--genomes", str(genomes), "--diamond-db", str(db),
                   "--protein-to-genome", str(mapping), "--assembly-info", str(metadata), "--taxdump", str(dump),
                   "--outdir", str(out), "--meta", "--threads", "2", "--sample-jobs", "2", "--reference-dir", str(references), "--no-download"]
        first = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(first.returncode, 0, first.stderr + (out / "logs" / "compare.log").read_text() if (out / "logs" / "compare.log").exists() else first.stderr)
        self.assertEqual(json.loads((out / "pipeline_manifest.json").read_text())["status"], "complete")
        with (out / "cohort_ani" / "nearest_references.csv").open() as f:
            nearest = list(csv.DictReader(f))
        self.assertEqual(len(nearest), 2)
        self.assertTrue(all(r["nearest_reference"] == reference for r in nearest))
        self.assertTrue(all(float(r["estimated_ani"]) > 99 for r in nearest))
        times = {}
        for sample in ("sample_a", "sample_b"):
            manifest = json.loads((out / "samples" / sample / "run_manifest.json").read_text())
            self.assertEqual(manifest["arguments"]["threads"], 1)
            for name in ("prodigal.faa", "diamond.tsv"):
                path = out / "samples" / sample / name
                times[path] = path.stat().st_mtime_ns
            fields, rows = compare.read_csv(out / "spraynpray_combined" / (sample + ".csv"))
            self.assertIn("protein_accession", fields)
            self.assertGreater(len(rows), 0)
        ani = out / "cohort_ani" / "work" / "fastani.tsv"
        times[ani] = ani.stat().st_mtime_ns
        again = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(again.returncode, 0, again.stderr)
        self.assertTrue(all(path.stat().st_mtime_ns == timestamp for path, timestamp in times.items()))
        # Removing a sample must not leave its old combined CSV in the next cohort.
        (genomes / "sample_b.fna.gz").unlink()
        final = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(final.returncode, 0, final.stderr)
        self.assertFalse((out / "spraynpray_combined" / "sample_b.csv").exists())
        self.assertTrue((out / "samples" / "sample_b").is_dir())
        with (out / "cohort_ani" / "nearest_references.csv").open() as f:
            self.assertEqual([r["sample"] for r in csv.DictReader(f)], ["sample_a"])


if __name__ == "__main__":
    unittest.main()
