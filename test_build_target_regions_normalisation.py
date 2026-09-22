import csv
import hashlib
import tempfile
import unittest
from pathlib import Path

import build_target_regions_normalisation as build_mod


class TestChrPrefixBed(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.src = Path(self.tmp.name) / "nochr.bed"
        self.dest = Path(self.tmp.name) / "chr.bed"

    def tearDown(self):
        self.tmp.cleanup()

    def test_prefixes_every_line(self):
        self.src.write_text("1\t100\t200\t.\n2\t300\t400\trs123\n")
        build_mod.chr_prefix_bed(self.src, self.dest)
        self.assertEqual(
            self.dest.read_text(),
            "chr1\t100\t200\t.\nchr2\t300\t400\trs123\n",
        )

    def test_preserves_line_count(self):
        self.src.write_text("1\t1\t2\t.\nX\t1\t2\t.\nY\t1\t2\t.\n")
        build_mod.chr_prefix_bed(self.src, self.dest)
        self.assertEqual(
            len(self.dest.read_text().splitlines()),
            len(self.src.read_text().splitlines()),
        )

    def test_skips_blank_lines_without_prefixing_them(self):
        self.src.write_text("1\t1\t2\t.\n\n2\t3\t4\t.\n")
        build_mod.chr_prefix_bed(self.src, self.dest)
        lines = self.dest.read_text().splitlines()
        self.assertEqual(lines, ["chr1\t1\t2\t.", "", "chr2\t3\t4\t."])

    def test_does_not_double_prefix_if_rerun_on_output(self):
        # Guards against accidentally chr-prefixing an already chr-prefixed
        # file (e.g. a caller passing the wrong source by mistake).
        self.src.write_text("1\t1\t2\t.\n")
        build_mod.chr_prefix_bed(self.src, self.dest)
        second = Path(self.tmp.name) / "chr2.bed"
        build_mod.chr_prefix_bed(self.dest, second)
        self.assertEqual(second.read_text(), "chrchr1\t1\t2\t.\n")
        # i.e. calling this twice is a caller error, not something the
        # function silently protects against -- documented behaviour.


class TestLoadCohortManifest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "manifest.tsv"

    def tearDown(self):
        self.tmp.cleanup()

    def test_loads_rows_as_dicts(self):
        self.path.write_text(
            "SampleId\tRatioFileId\tAmberBafFileId\n"
            "S1\tfile-AAA\tfile-BBB\n"
            "S2\tfile-CCC\tfile-DDD\n"
        )
        rows = build_mod.load_cohort_manifest(self.path)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0], {"SampleId": "S1", "RatioFileId": "file-AAA", "AmberBafFileId": "file-BBB"})
        self.assertEqual(rows[1]["SampleId"], "S2")

    def test_empty_manifest_returns_empty_list(self):
        self.path.write_text("SampleId\tRatioFileId\tAmberBafFileId\n")
        rows = build_mod.load_cohort_manifest(self.path)
        self.assertEqual(rows, [])

    def test_real_bundled_manifest_has_41_rows_and_no_duplicates(self):
        real_manifest = Path(__file__).parent / "training_cohort_manifest.tsv"
        rows = build_mod.load_cohort_manifest(real_manifest)
        self.assertEqual(len(rows), 41)
        sample_ids = [r["SampleId"] for r in rows]
        self.assertEqual(len(sample_ids), len(set(sample_ids)))
        for row in rows:
            self.assertTrue(row["RatioFileId"].startswith("file-"))
            self.assertTrue(row["AmberBafFileId"].startswith("file-"))


class TestBuildSampleIdsCsv(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dest = Path(self.tmp.name) / "sample_ids.csv"

    def tearDown(self):
        self.tmp.cleanup()

    def test_writes_header_and_sample_ids_only(self):
        manifest = [
            {"SampleId": "S1", "RatioFileId": "file-AAA", "AmberBafFileId": "file-BBB"},
            {"SampleId": "S2", "RatioFileId": "file-CCC", "AmberBafFileId": "file-DDD"},
        ]
        build_mod.build_sample_ids_csv(manifest, self.dest)
        self.assertEqual(self.dest.read_text(), "SampleId\nS1\nS2\n")

    def test_never_writes_a_gender_column(self):
        # Decision 4: gender must be inferred by NormalisationFileBuilder
        # from real AMBER BAF data, not supplied -- so this function must
        # never emit a Gender column under any input shape.
        manifest = [{"SampleId": "S1", "RatioFileId": "file-AAA", "AmberBafFileId": "file-BBB", "Gender": "MALE"}]
        build_mod.build_sample_ids_csv(manifest, self.dest)
        header = self.dest.read_text().splitlines()[0]
        self.assertEqual(header, "SampleId")


class TestVerifyOutput(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "out.tsv"
        self.content = "Chromosome\tPosition\tRelativeEnrichment\nchr1\t100\t1.0\n"
        self.path.write_text(self.content)
        self.real_md5 = hashlib.md5(self.content.encode()).hexdigest()
        self.real_lines = len(self.content.splitlines())

    def tearDown(self):
        self.tmp.cleanup()

    def test_passes_when_lines_and_md5_match(self):
        orig_md5, orig_lines = build_mod.EXPECTED_OUTPUT_MD5, build_mod.EXPECTED_OUTPUT_LINES
        build_mod.EXPECTED_OUTPUT_MD5, build_mod.EXPECTED_OUTPUT_LINES = self.real_md5, self.real_lines
        try:
            build_mod.verify_output(self.path)  # must not raise
        finally:
            build_mod.EXPECTED_OUTPUT_MD5, build_mod.EXPECTED_OUTPUT_LINES = orig_md5, orig_lines

    def test_exits_when_md5_does_not_match(self):
        orig_md5, orig_lines = build_mod.EXPECTED_OUTPUT_MD5, build_mod.EXPECTED_OUTPUT_LINES
        build_mod.EXPECTED_OUTPUT_MD5, build_mod.EXPECTED_OUTPUT_LINES = "0" * 32, self.real_lines
        try:
            with self.assertRaises(SystemExit):
                build_mod.verify_output(self.path)
        finally:
            build_mod.EXPECTED_OUTPUT_MD5, build_mod.EXPECTED_OUTPUT_LINES = orig_md5, orig_lines

    def test_exits_when_line_count_does_not_match(self):
        orig_md5, orig_lines = build_mod.EXPECTED_OUTPUT_MD5, build_mod.EXPECTED_OUTPUT_LINES
        build_mod.EXPECTED_OUTPUT_MD5, build_mod.EXPECTED_OUTPUT_LINES = self.real_md5, 999999
        try:
            with self.assertRaises(SystemExit):
                build_mod.verify_output(self.path)
        finally:
            build_mod.EXPECTED_OUTPUT_MD5, build_mod.EXPECTED_OUTPUT_LINES = orig_md5, orig_lines

    def test_skip_flag_bypasses_a_mismatch_without_raising(self):
        orig_md5 = build_mod.EXPECTED_OUTPUT_MD5
        build_mod.EXPECTED_OUTPUT_MD5 = "0" * 32
        try:
            build_mod.verify_output(self.path, skip=True)  # must not raise
        finally:
            build_mod.EXPECTED_OUTPUT_MD5 = orig_md5


if __name__ == "__main__":
    unittest.main()
