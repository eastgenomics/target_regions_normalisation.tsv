import csv
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import build_target_regions_normalisation as build_mod


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
        with mock.patch.object(build_mod, "EXPECTED_OUTPUT_MD5", self.real_md5), \
             mock.patch.object(build_mod, "EXPECTED_OUTPUT_LINES", self.real_lines):
            build_mod.verify_output(self.path)  # must not raise

    def test_exits_when_md5_does_not_match(self):
        with mock.patch.object(build_mod, "EXPECTED_OUTPUT_MD5", "0" * 32), \
             mock.patch.object(build_mod, "EXPECTED_OUTPUT_LINES", self.real_lines):
            with self.assertRaises(SystemExit):
                build_mod.verify_output(self.path)

    def test_exits_when_line_count_does_not_match(self):
        with mock.patch.object(build_mod, "EXPECTED_OUTPUT_MD5", self.real_md5), \
             mock.patch.object(build_mod, "EXPECTED_OUTPUT_LINES", 999999):
            with self.assertRaises(SystemExit):
                build_mod.verify_output(self.path)

    def test_skip_flag_bypasses_a_mismatch_without_raising(self):
        with mock.patch.object(build_mod, "EXPECTED_OUTPUT_MD5", "0" * 32):
            build_mod.verify_output(self.path, skip=True)  # must not raise


class TestRunNormalisationFileBuilderCommandConstruction(unittest.TestCase):
    """Verifies the actual argv list passed to sh() -- no real java/dx call."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.workdir = Path(self.tmp.name)
        self.captured = []
        self.orig_sh = build_mod.sh

        def fake_sh(cmd, **kw):
            self.captured.append((cmd, kw))

        build_mod.sh = fake_sh

    def tearDown(self):
        build_mod.sh = self.orig_sh
        self.tmp.cleanup()

    def test_uses_shell_false_style_argv_list_not_a_string(self):
        build_mod.run_normalisation_file_builder(
            self.workdir / "cobalt.jar", self.workdir / "bed.bed", self.workdir / "gc.cnp",
            "38", self.workdir / "cobalt_dir", self.workdir / "amber_dir",
            self.workdir / "sample_ids.csv", self.workdir / "out.tsv", self.workdir,
        )
        cmd, kw = self.captured[0]
        self.assertIsInstance(cmd, list)
        self.assertTrue(all(isinstance(c, str) for c in cmd))
        self.assertEqual(kw.get("cwd"), self.workdir)

    def test_output_file_arg_is_the_resolved_absolute_path_even_with_a_subdirectory(self):
        # Regression test: -output_file must not be reduced to a bare
        # filename, or an --out with a subdirectory silently writes to the
        # wrong place (verify_output would then look in the subdirectory
        # and raise FileNotFoundError).
        nested_out = self.workdir / "results" / "out.tsv"
        build_mod.run_normalisation_file_builder(
            self.workdir / "cobalt.jar", self.workdir / "bed.bed", self.workdir / "gc.cnp",
            "38", self.workdir / "cobalt_dir", self.workdir / "amber_dir",
            self.workdir / "sample_ids.csv", nested_out, self.workdir,
        )
        cmd, _ = self.captured[0]
        idx = cmd.index("-output_file")
        self.assertEqual(cmd[idx + 1], str(nested_out.resolve()))

    def test_other_inputs_are_passed_as_bare_names_relative_to_cwd(self):
        build_mod.run_normalisation_file_builder(
            self.workdir / "cobalt.jar", self.workdir / "bed.bed", self.workdir / "gc.cnp",
            "38", self.workdir / "cobalt_dir", self.workdir / "amber_dir",
            self.workdir / "sample_ids.csv", self.workdir / "out.tsv", self.workdir,
        )
        cmd, _ = self.captured[0]
        self.assertIn("cobalt.jar", cmd)
        self.assertIn("bed.bed", cmd)
        self.assertIn("cobalt_dir/", cmd)


class TestMainWorkdirHandling(unittest.TestCase):
    def test_skip_download_without_explicit_workdir_exits(self):
        orig_argv = sys.argv
        sys.argv = ["prog", "--skip-download"]
        try:
            with self.assertRaises(SystemExit):
                build_mod.main()
        finally:
            sys.argv = orig_argv

    def test_empty_workdir_string_is_rejected_not_treated_as_omitted(self):
        orig_argv = sys.argv
        sys.argv = ["prog", "--workdir", ""]
        try:
            with self.assertRaises(SystemExit) as ctx:
                build_mod.main()
            self.assertEqual(ctx.exception.code, "--workdir must not be empty")
        finally:
            sys.argv = orig_argv

    def test_empty_workdir_string_with_skip_download_is_rejected(self):
        orig_argv = sys.argv
        sys.argv = ["prog", "--skip-download", "--workdir", ""]
        try:
            with self.assertRaises(SystemExit) as ctx:
                build_mod.main()
            self.assertEqual(ctx.exception.code, "--workdir must not be empty")
        finally:
            sys.argv = orig_argv


if __name__ == "__main__":
    unittest.main()
