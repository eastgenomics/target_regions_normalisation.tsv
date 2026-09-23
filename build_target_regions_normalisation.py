#!/usr/bin/env python3
"""
build_target_regions_normalisation.py

Builds target_regions_normalisation.tsv — the COBALT panel-normalisation file
used by eggd_cgp-cobalt to correct per-window read-depth ratios for CGP
backbone target-region biases. Built once from a fixed 41-sample training
cohort ("EF v1") via Hartwig's own NormalisationFileBuilder
(com.hartwig.hmftools.cobalt.norm.NormalisationFileBuilder, ships inside
cobalt.jar) — this is Hartwig's documented "Panel resources training
procedure" (see pipeline/README_TARGETED.md in hartwigmedical/hmftools),
not a bespoke process, executed against our own cohort data.

This is a one-off, panel-setup-time build, not a per-sample or per-run
pipeline step. Re-run only if the panel design, training cohort, or
reference-genome convention changes.

See DECISIONS below for every manual/judgement-based choice this script
makes and why — each is backed by direct evidence gathered by tracing
NormalisationFileBuilder's actual source, and by controlled A/B reruns
against the historically-built, currently-live file.

Inputs (DNAnexus file IDs, immutable/content-addressed):
  --cobalt-jar-source     project-Fkb6Gkj433GVVvj73J7x8KbV:file-J893p9Q470j4zY3zzpVBjP11
                            (cobalt.jar, COBALT 3.0-beta.5 — the APPROVED,
                            documented production jar; see DECISIONS)
  --backbone-bed-source   project-Fkb6Gkj433GVVvj73J7x8KbV:file-JBp3vx8433GQbP2bJY9J9X7V
                            (backbone_padded_150bp.bed, chr-prefixed — padded
                            directly from the chr-prefixed raw source, see
                            DECISIONS)
  --gc-profile-source     file-J88xxvQ4QyVPb8K6VFqX1FKB
                            (GC_profile.1000bp.38.cnp, already chr-prefixed,
                            no swap needed)
  --cohort-manifest       training_cohort_manifest.tsv (bundled with this
                            script) — the frozen 41-sample EF v1 cohort:
                            SampleId, RatioFileId, AmberBafFileId, all in
                            project-J88p7V0470jPZ2Vv1VqB3BFz

Output:
  target_regions_normalisation.tsv

Requires: dxpy (DNAnexus Python API, authenticated via the same
credentials as the dx CLI), java (JRE, no Docker needed —
NormalisationFileBuilder is a plain JAR invocation), python3 stdlib only
otherwise.
"""
import argparse
import csv
import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

import dxpy

DEFAULT_COBALT_JAR_SOURCE = "project-Fkb6Gkj433GVVvj73J7x8KbV:file-J893p9Q470j4zY3zzpVBjP11"
DEFAULT_BACKBONE_BED_SOURCE = "project-Fkb6Gkj433GVVvj73J7x8KbV:file-JBp3vx8433GQbP2bJY9J9X7V"
DEFAULT_GC_PROFILE_SOURCE = "file-J88xxvQ4QyVPb8K6VFqX1FKB"
DEFAULT_REF_GENOME_VERSION = "38"

# Cohort project — all ratio/AMBER-BAF files in the bundled manifest live here.
COHORT_PROJECT = "project-J88p7V0470jPZ2Vv1VqB3BFz"

EXPECTED_OUTPUT_MD5 = "ab8a12cf4fa129e7c227c6fdb3d0612f"  # pinned 2026-09-22, confirmed
# reproducible via two independent from-scratch runs (a manual/agent-driven
# rebuild, then this script). 74,502 data rows + 1 header = 74,503 lines.
EXPECTED_OUTPUT_LINES = 74503

# ---------------------------------------------------------------------------
# DECISIONS — every manual, judgement-based resolution from the investigation
# that produced this script, baked in as explicit data rather than left as
# one-off interactive commands. Full rationale and evidence for each is in
# README.md's "Decisions" section, not repeated here.
# ---------------------------------------------------------------------------

# Decision 1 — Use the approved COBALT 3.0-beta.5 jar
# (file-J893p9Q470j4zY3zzpVBjP11), not the older, unapproved jar the
# historical build actually used (Confluence DV 4807622662, DI-3865).

# Decision 2 — Use -ref_genome_version 38, never 37 -- required by the
# approved jar, and the BED's chr-prefix convention must match it anyway.

# Decision 3 — use a permanently-stored chr-prefixed BED
# (backbone_padded_150bp.bed, padded directly from the chr-prefixed raw
# source) rather than chr-prefixing a no-chr sibling at build time or
# reusing a separate no-chr GC profile -- GC-profile matching isn't
# chr-prefix agnostic, so target_regions_bed must match its convention.

# Decision 4 — no Gender column in the sample manifest; gender is inferred
# by NormalisationFileBuilder from real AMBER BAF data. A 41-sample A/B
# rerun confirmed the V37/V38 PAR-boundary difference changes zero gender
# calls for this cohort.

# Decision 5 — this 41-sample cohort manifest is frozen, versioned data
# (matches the currently-live production file), not re-derived on the fly.


def sh(cmd, **kw):
    print(f"+ {' '.join(str(c) for c in cmd)}", file=sys.stderr)
    return subprocess.run(cmd, shell=False, check=True, **kw)


def dx_download(project_file_id, out_path):
    # Accepts either "project-XXXX:file-YYYY" or a bare "file-YYYY" -- the
    # same two forms `dx download` itself accepts.
    project_id, _, file_id = str(project_file_id).rpartition(":")
    dxpy.download_dxfile(file_id, str(out_path), project=project_id or None)


def load_cohort_manifest(path):
    """Returns a list of dicts: SampleId, RatioFileId, AmberBafFileId."""
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def build_sample_ids_csv(manifest, dest_path):
    """SampleId-only CSV, deliberately no Gender column -- see Decision 4."""
    with open(dest_path, "w") as f:
        f.write("SampleId\n")
        for row in manifest:
            f.write(row["SampleId"] + "\n")


def download_cohort_files(manifest, cobalt_dir, amber_dir):
    cobalt_dir, amber_dir = Path(cobalt_dir), Path(amber_dir)
    cobalt_dir.mkdir(parents=True, exist_ok=True)
    amber_dir.mkdir(parents=True, exist_ok=True)
    for row in manifest:
        sample = row["SampleId"]
        dx_download(
            f"{COHORT_PROJECT}:{row['RatioFileId']}",
            cobalt_dir / f"{sample}.cobalt.ratio.tsv.gz",
        )
        dx_download(
            f"{COHORT_PROJECT}:{row['AmberBafFileId']}",
            amber_dir / f"{sample}.amber.baf.tsv.gz",
        )


def run_normalisation_file_builder(
    cobalt_jar, backbone_bed, gc_profile, ref_genome_version,
    cobalt_dir, amber_dir, sample_ids_csv, out_path, workdir,
):
    # -output_file is the one argument that may legitimately point outside
    # workdir (a caller-supplied --out with a subdirectory) -- resolve it to
    # an absolute path rather than reducing it to a bare filename, or the
    # builder would silently write to the wrong place (see the test that
    # covers this: test_out_path_with_subdirectory_is_preserved).
    sh(
        [
            "java", "-cp", Path(cobalt_jar).name,
            "com.hartwig.hmftools.cobalt.norm.NormalisationFileBuilder",
            "-cobalt_dir", f"{Path(cobalt_dir).name}/",
            "-target_regions_bed", Path(backbone_bed).name,
            "-gc_profile", Path(gc_profile).name,
            "-ref_genome_version", str(ref_genome_version),
            "-sample_id_file", Path(sample_ids_csv).name,
            "-amber_dir", f"{Path(amber_dir).name}/",
            "-output_file", str(Path(out_path).resolve()),
        ],
        cwd=workdir,
    )


def verify_output(out_path, skip=False):
    out_path = Path(out_path)
    content = out_path.read_bytes()
    lines = content.count(b"\n") + (1 if content and not content.endswith(b"\n") else 0)
    actual_md5 = hashlib.md5(content).hexdigest()
    if skip:
        print(f"Output checksum check skipped (--skip-checksum-assert): "
              f"{lines} lines, md5 {actual_md5}", file=sys.stderr)
        return
    if lines != EXPECTED_OUTPUT_LINES or actual_md5 != EXPECTED_OUTPUT_MD5:
        sys.exit(
            f"Output mismatch: expected {EXPECTED_OUTPUT_LINES} lines / md5 "
            f"{EXPECTED_OUTPUT_MD5}, got {lines} lines / md5 {actual_md5}.\n"
            f"This would silently change a resource this pipeline depends "
            f"on, so the build stops here instead. If this is an intentional "
            f"rebuild (new cohort, panel change, jar update), review the "
            f"resulting diff against the previous file, then update "
            f"EXPECTED_OUTPUT_MD5/EXPECTED_OUTPUT_LINES. To proceed once "
            f"without updating the pin, pass --skip-checksum-assert."
        )
    print(f"Output verified: {lines} lines, md5 {actual_md5}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cobalt-jar-source", default=DEFAULT_COBALT_JAR_SOURCE)
    ap.add_argument("--backbone-bed-source", default=DEFAULT_BACKBONE_BED_SOURCE)
    ap.add_argument("--gc-profile-source", default=DEFAULT_GC_PROFILE_SOURCE)
    ap.add_argument("--ref-genome-version", default=DEFAULT_REF_GENOME_VERSION,
                     help="Do not change to 37 -- see Decision 2 in the module docstring")
    ap.add_argument(
        "--cohort-manifest",
        default=str(Path(__file__).parent / "training_cohort_manifest.tsv"),
        help="Frozen 41-sample EF v1 cohort manifest, see Decision 5",
    )
    ap.add_argument("--workdir", default=None,
                     help="Defaults to a fresh, private temp dir (tempfile.mkdtemp) -- "
                     "pass an explicit path to reuse one, e.g. with --skip-download")
    ap.add_argument("--out", default="target_regions_normalisation.tsv")
    ap.add_argument("--skip-download", action="store_true", help="Reuse existing files in --workdir")
    ap.add_argument("--skip-checksum-assert", action="store_true")
    args = ap.parse_args()

    if args.workdir is not None:
        if args.workdir == "":
            sys.exit("--workdir must not be empty")
        workdir = Path(args.workdir)
        workdir.mkdir(parents=True, exist_ok=True)
    else:
        if args.skip_download:
            sys.exit("--skip-download requires an explicit --workdir (nothing to reuse otherwise)")
        # A fixed, world-writable /tmp path is an insecure default -- another
        # local user could pre-create it or race to replace cobalt.jar
        # between download and the java invocation. A fresh private tempdir
        # closes that off; pass --workdir explicitly to reuse a location.
        workdir = Path(tempfile.mkdtemp(prefix="build_target_regions_normalisation-"))
        print(f"Using generated workdir: {workdir}", file=sys.stderr)

    manifest = load_cohort_manifest(args.cohort_manifest)

    cobalt_jar = workdir / "cobalt.jar"
    backbone_bed = workdir / "backbone_padded_150bp.bed"
    gc_profile = workdir / "GC_profile.1000bp.38.cnp"
    cobalt_dir = workdir / "cobalt_bootstrap_all"
    amber_dir = workdir / "amber_flat"
    sample_ids_csv = workdir / "sample_ids.csv"
    out_path = (workdir / args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not args.skip_download:
        dx_download(args.cobalt_jar_source, cobalt_jar)
        dx_download(args.backbone_bed_source, backbone_bed)
        dx_download(args.gc_profile_source, gc_profile)
        download_cohort_files(manifest, cobalt_dir, amber_dir)
    else:
        for p in (cobalt_jar, backbone_bed, gc_profile, cobalt_dir, amber_dir):
            if not p.exists():
                sys.exit(f"--skip-download given but {p} is missing")

    build_sample_ids_csv(manifest, sample_ids_csv)

    run_normalisation_file_builder(
        cobalt_jar, backbone_bed, gc_profile, args.ref_genome_version,
        cobalt_dir, amber_dir, sample_ids_csv, out_path, workdir,
    )

    verify_output(out_path, skip=args.skip_checksum_assert)
    print(f"Built {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
