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
  --backbone-bed-source   project-Fkb6Gkj433GVVvj73J7x8KbV:file-J88gVF84Y8X123K6JX8jB8Z5
                            (backbone_padded_150bp_nochr.bed — no-chr; this
                            script chr-prefixes it itself, see DECISIONS)
  --gc-profile-source     file-J88xxvQ4QyVPb8K6VFqX1FKB
                            (GC_profile.1000bp.38.cnp, already chr-prefixed,
                            no swap needed)
  --cohort-manifest       training_cohort_manifest.tsv (bundled with this
                            script) — the frozen 41-sample EF v1 cohort:
                            SampleId, RatioFileId, AmberBafFileId, all in
                            project-J88p7V0470jPZ2Vv1VqB3BFz

Output:
  target_regions_normalisation.tsv

Requires: dx (DNAnexus CLI, authenticated), java (JRE, no Docker needed —
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

DEFAULT_COBALT_JAR_SOURCE = "project-Fkb6Gkj433GVVvj73J7x8KbV:file-J893p9Q470j4zY3zzpVBjP11"
DEFAULT_BACKBONE_BED_SOURCE = "project-Fkb6Gkj433GVVvj73J7x8KbV:file-J88gVF84Y8X123K6JX8jB8Z5"
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
# that produced this script, baked in as explicit data/comments rather than
# left as one-off interactive commands. Each entry records WHY and the
# evidence, so a reviewer doesn't have to re-run the original investigation.
# ---------------------------------------------------------------------------

# Decision 1 — Use COBALT 3.0-beta.5 (file-J893p9Q470j4zY3zzpVBjP11), not the
# older jar (file-J88y0Jj4VP405z36Jq1zyQQj, hmf-common-cobalt-2.2) that the
# historical build actually used. 3.0-beta.5 is the APPROVED, documented
# production jar for this whole pipeline (Confluence DV space, page
# 4807622662, APPROVED 2026-09-14, DI-3865) — verified byte-for-byte against
# the official GitHub release asset. Using a different, unapproved jar
# version to build a resource this pipeline depends on would itself be a
# provenance problem, independent of anything else in this file.

# Decision 2 — Use -ref_genome_version 38, never 37. Two independent reasons:
#   (a) The approved jar's own documented behaviour (same Confluence page,
#       4807622662): "-ref_genome_version 38 is required -- 37 causes an
#       IndexOutOfBoundsException in this build." The historical build could
#       only use -ref_genome_version 37 because it used a *different*,
#       unapproved jar that happened to tolerate it.
#   (b) Even where 37 doesn't crash, it's a formatting trick, not a genome
#       assembly choice: NormalisationFileBuilder's FileWriter reconstructs
#       each expected chromosome key via
#       RefGenomeVersion.versionedChromosome() -- hmftools' own convention,
#       V37 -> no-chr, V38 -> chr-prefixed -- and looks that key up against
#       the BED-derived region map. Passing 38 against a no-chr BED silently
#       skips every chromosome (zero output rows, no error); passing 37
#       against a chr-prefixed BED does the same. The BED's own chr
#       convention must match whatever -ref_genome_version's own convention
#       produces.

# Decision 3 — chr-prefix the backbone BED ourselves (from the live no-chr
# source, backbone_padded_150bp_nochr.bed) rather than reuse a separate
# no-chr GC profile as the historical build did. Two reasons:
#   (a) GcProfileCache.findGcProfile() does a raw, unnormalised string-map
#       lookup keyed by whatever the GC profile file literally says -- it is
#       NOT chr-prefix agnostic (unlike the COBALT-ratio-matching path,
#       which normalises via HumanChromosome.fromString() and tolerates
#       either convention). Pairing a no-chr BED with the live chr-prefixed
#       GC_profile.1000bp.38.cnp throws a NullPointerException in this exact
#       method -- confirmed by direct reproduction with the real jar.
#   (b) The no-chr GC profile the historical build used
#       (GC_profile.1000bp.38.nochr.cnp) was independently confirmed
#       byte-identical to the live chr-prefixed copy except for the literal
#       "chr" prefix on every line (md5 match after re-stripping, zero diff
#       lines) -- so there was never a need for two separate copies. Simplest
#       correct fix: chr-prefix the BED, keep the one GC profile already in
#       resource_ids.env.

# Decision 4 — no Gender column in the sample manifest; gender for all 41
# samples is inferred by NormalisationFileBuilder itself from real AMBER BAF
# chrX heterozygosity data (com.hartwig.hmftools.common.amber.AmberGender),
# using whichever pseudo-autosomal-region (PAR) boundary set matches
# -ref_genome_version. The historical build (V37 PAR boundaries against
# genuinely GRCh38-aligned AMBER BAF data) was a real, if theoretical,
# coordinate mismatch -- confirmed via a full 41-sample controlled A/B rerun
# to have changed ZERO gender calls and ZERO output values for this cohort
# (holding the jar version constant). V38 is still the correct choice
# because it's genuinely consistent, not because the V37 mismatch caused
# visible harm here.

# Decision 5 — this cohort manifest (41 EF v1 samples) is frozen, versioned
# data, not re-derived. It matches the samples used for the currently-live
# production file, letting a rebuild be checked against that file directly.
# A different/larger training cohort is a deliberate, separate decision
# (new panel validation), not something this script infers on its own.


def sh(cmd, **kw):
    print(f"+ {' '.join(str(c) for c in cmd)}", file=sys.stderr)
    return subprocess.run(cmd, shell=False, check=True, **kw)


def dx_download(project_file_id, out_path):
    sh(["dx", "download", str(project_file_id), "-o", str(out_path), "-f"])


def chr_prefix_bed(src_path, dest_path):
    """Prepend 'chr' to every line's chromosome column. This BED has no
    header/comment lines, so a blind line-start substitution is safe and
    unambiguous -- confirmed against the real file (57,295 lines in, 57,295
    lines out, only column 1 changes)."""
    src_path, dest_path = Path(src_path), Path(dest_path)
    with src_path.open() as fin, dest_path.open("w") as fout:
        for line in fin:
            if line.strip():
                fout.write("chr" + line)
            else:
                fout.write(line)


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
    cobalt_jar, backbone_bed_chr, gc_profile, ref_genome_version,
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
            "-target_regions_bed", Path(backbone_bed_chr).name,
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
    with out_path.open() as f:
        lines = sum(1 for _ in f)
    actual_md5 = hashlib.md5(out_path.read_bytes()).hexdigest()
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

    if args.workdir:
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

    manifest = load_cohort_manifest(args.cohort_manifest)

    cobalt_jar = workdir / "cobalt.jar"
    backbone_bed_nochr = workdir / "backbone_padded_150bp_nochr.bed"
    backbone_bed_chr = workdir / "backbone_padded_150bp_chr.bed"
    gc_profile = workdir / "GC_profile.1000bp.38.cnp"
    cobalt_dir = workdir / "cobalt_bootstrap_all"
    amber_dir = workdir / "amber_flat"
    sample_ids_csv = workdir / "sample_ids.csv"
    out_path = (workdir / args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if not args.skip_download:
        dx_download(args.cobalt_jar_source, cobalt_jar)
        dx_download(args.backbone_bed_source, backbone_bed_nochr)
        dx_download(args.gc_profile_source, gc_profile)
        download_cohort_files(manifest, cobalt_dir, amber_dir)
    else:
        for p in (cobalt_jar, backbone_bed_nochr, gc_profile, cobalt_dir, amber_dir):
            if not p.exists():
                sys.exit(f"--skip-download given but {p} is missing")

    chr_prefix_bed(backbone_bed_nochr, backbone_bed_chr)
    build_sample_ids_csv(manifest, sample_ids_csv)

    run_normalisation_file_builder(
        cobalt_jar, backbone_bed_chr, gc_profile, args.ref_genome_version,
        cobalt_dir, amber_dir, sample_ids_csv, out_path, workdir,
    )

    verify_output(out_path, skip=args.skip_checksum_assert)
    print(f"Built {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
