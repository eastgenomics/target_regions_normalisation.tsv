# target_regions_normalisation.tsv

Reproducible build script for `target_regions_normalisation_v1.1.0.tsv` — the
COBALT panel-normalisation file used by `eggd_cgp-cobalt` to correct
per-window read-depth ratios for CGP backbone target-region biases.

**Status:** the v1.1.0 controlled document
v1.1.0 ([DV 4875714639](https://cuhbioinformatics.atlassian.net/wiki/spaces/DV/pages/4875714639))
is now the approved version for use in production, and can be found at the following location (001_Reference:/bed_files/b38/atlas_cnv/target_regions_normalisation_v1.1.0.tsv), with file id (file-JFv3QYj433GkYzp9zfGKkGpY). v1.0.0
([DV 4805623939](https://cuhbioinformatics.atlassian.net/wiki/spaces/DV/pages/4805623939))
has been removed from the 001_Reference folder.

## Versions

| Version | Training cohort | Data | md5 |
|---|---|---|---|
| 1.1.0 (this release) | 123 CNVkit static PoN members (FFPE tumours) | non-UMI production BAMs | `1c5c976586199830c509ac5dc19d8491` |
| 1.0.0 | 41 EF v1 samples (run 0734 Val1) | UMI-consensus BAMs | `ab8a12cf4fa129e7c227c6fdb3d0612f` |

Both are 74,502 data rows over the same windows. v1.1.0 is built from the data
production now uses (UMIs masked at demultiplexing, duplicates marked; DI-4114)
and from the same sample type (FFPE). Jira: DI-4246.

## What this builds

`target_regions_normalisation_v1.1.0.tsv` is built from a fixed 123-sample
training cohort, by running Hartwig's own
`com.hartwig.hmftools.cobalt.norm.NormalisationFileBuilder` (a class shipped
inside `cobalt.jar`) — this is Hartwig's documented ["Panel resources
training procedure"](https://github.com/hartwigmedical/hmftools/blob/master/pipeline/README_TARGETED.md),
not a bespoke process, executed against our own cohort's COBALT ratio and
AMBER BAF outputs.

This is a **one-off, panel-setup-time build**, not a routine or per-sample
pipeline step. It's not expected to ever need rebuilding with this exact
tool version — the next real rebuild will happen against a mature (non-beta)
COBALT release, which is a bigger, separate decision than anything this
script controls.

The script:

1. Downloads `cobalt.jar` (the approved production release, COBALT
   3.0-beta.5), the chr-prefixed backbone target-regions BED, and the GC
   profile from DNAnexus, plus the 123-sample training cohort's COBALT ratio
   and AMBER BAF files (file IDs frozen in `training_cohort_manifest.tsv`).
2. Builds a `SampleId`-only sample list — no Gender column, so gender is
   inferred by `NormalisationFileBuilder` itself from real AMBER BAF data.
3. Runs `NormalisationFileBuilder` and asserts the result against a pinned
   checksum.

The Confluence controlled documents record each version's provenance and
testing:
[target_regions_normalisation_v1.1.0.tsv](https://cuhbioinformatics.atlassian.net/wiki/spaces/DV/pages/4875714639)
and, for v1.0.0 (including the investigation behind Decisions 1-4),
[target_regions_normalisation.tsv](https://cuhbioinformatics.atlassian.net/wiki/spaces/DV/pages/4805623939).

## Decisions

Every manual, judgement-based resolution from the investigation that
produced this script, with the full rationale and evidence. The
`DECISIONS` block in `build_target_regions_normalisation.py` records the
same six items in short form, pointing back here.

### Decision 1 — use the approved COBALT jar

Use COBALT 3.0-beta.5 (`file-J893p9Q470j4zY3zzpVBjP11`), not the older jar
(`file-J88y0Jj4VP405z36Jq1zyQQj`, `hmf-common-cobalt-2.2`) that the
historical build actually used. 3.0-beta.5 is the APPROVED, documented
production jar for this whole pipeline (Confluence DV space, page
4807622662, APPROVED 2026-09-14, DI-3865) — verified byte-for-byte against
the official GitHub release asset. Using a different, unapproved jar
version to build a resource this pipeline depends on would itself be a
provenance problem, independent of anything else in this file.

### Decision 2 — always `-ref_genome_version 38`, never `37`

Two independent reasons:

- The approved jar's own documented behaviour (same Confluence page,
  4807622662): "`-ref_genome_version 38` is required -- `37` causes an
  `IndexOutOfBoundsException` in this build." The historical build could
  only use `-ref_genome_version 37` because it used a *different*,
  unapproved jar that happened to tolerate it.
- Even where `37` doesn't crash, it's a formatting trick, not a genome
  assembly choice: `NormalisationFileBuilder`'s `FileWriter` reconstructs
  each expected chromosome key via `RefGenomeVersion.versionedChromosome()`
  -- hmftools' own convention, V37 → no-chr, V38 → chr-prefixed -- and
  looks that key up against the BED-derived region map. Passing `38`
  against a no-chr BED silently skips every chromosome (zero output rows,
  no error); passing `37` against a chr-prefixed BED does the same. The
  BED's own chr convention must match whatever `-ref_genome_version`'s own
  convention produces.

### Decision 3 — use a permanently-stored chr-prefixed BED

`backbone_padded_150bp.bed` is padded directly from the original
chr-prefixed raw source (`targets_covered_CNV_WG_s50kb_TE-96770770_hg38_annot.bed`,
the same vendor/design-provided SNP-position list `backbone_padded_150bp_nochr.bed`
was itself padded from) via `bedtools slop -b 150`, so it's chr-prefixed
from the start — no in-script chr-prefixing step, and no need to reuse a
separate no-chr GC profile, as the historical build did. Two reasons this
matters:

- `GcProfileCache.findGcProfile()` does a raw, unnormalised string-map
  lookup keyed by whatever the GC profile file literally says -- it is
  NOT chr-prefix agnostic (unlike the COBALT-ratio-matching path, which
  normalises via `HumanChromosome.fromString()` and tolerates either
  convention). Pairing a no-chr BED with the live chr-prefixed
  `GC_profile.1000bp.38.cnp` throws a `NullPointerException` in this exact
  method -- confirmed by direct reproduction with the real jar.
- The no-chr GC profile the historical build used
  (`GC_profile.1000bp.38.nochr.cnp`) was independently confirmed
  byte-identical to the live chr-prefixed copy except for the literal
  "chr" prefix on every line (md5 match after re-stripping, zero diff
  lines) -- so there was never a need for two separate copies.

This chr-prefixed BED is also used directly by PURPLE's own
`target_regions_bed` input elsewhere in the pipeline (unrelated to this
script, but the same underlying resource). PURPLE's own BED parsing was
traced to `hmf-common`'s `BedLine.chromosome()`
(`HumanChromosome.fromString()`, chr-prefix agnostic) and confirmed via a
real A/B PURPLE run on a production sample that swapping between the
no-chr and chr-prefixed forms of this BED produces byte-identical output
across every downstream file (purity, CNV somatic/gene, target-region CN,
QC) — so this change didn't require re-wiring PURPLE.

### Decision 4 — no Gender column in the sample manifest

Gender for all 123 samples is inferred by `NormalisationFileBuilder` itself
from real AMBER BAF chrX heterozygosity data
(`com.hartwig.hmftools.common.amber.AmberGender`), using whichever
pseudo-autosomal-region (PAR) boundary set matches `-ref_genome_version`.
The historical build (V37 PAR boundaries against genuinely
GRCh38-aligned AMBER BAF data) was a real, if theoretical, coordinate
mismatch -- confirmed via a full 41-sample controlled A/B rerun (the v1.0.0 cohort) to have
changed ZERO gender calls and ZERO output values for that cohort (holding
the jar version constant). V38 is still the correct choice because it's
genuinely consistent, not because the V37 mismatch caused visible harm
here.

### Decision 5 — the training cohort is 123 of the 124 FFPE CNVkit PoN members

The cohort is the members of the CNVkit static PoN (DI-4093), under the
same rules: FFPE tumours only (production is FFPE), one sequencing per
specimen, median depth of at least 100x. It excludes the commercial
controls, the fresh-frozen WGS truth tumours, their dilutions and the blood
germline normals, so none of the PURPLE clinical-validation truth samples is
in the training set. 123 samples is six times the 20 that Hartwig's training
procedure recommends.

One PoN member is excluded: 26127S0031, probably 47,XXY. It carries a Y
chromosome (laboratory record male; chrY reads at male level) and two X
chromosomes (62 of 84 informative chrX sites heterozygous; female chrX
dosage). `NormalisationFileBuilder` models sex chromosomes for two sexes
only: it divides each sample's chrX ratio by 1 (female) or 0.5 (male) and
leaves female samples out of chrY. Neither model fits this sample's sex
chromosomes, so it is left out rather than entered under either sex. The
file built with it included differs from this one by a median 0.0014 log2 on
chrX (maximum 0.032; no row switches between masked and unmasked), not at all
on chrY, and by a median 0.0013 log2 on the autosomes, where 21 of 70,020
rows switch between masked and unmasked.

The AMBER BAF files are the fgbio-clipped AMBER outputs (`eggd_cgp-amber`
1.0.1 on `ClipBam --clip-overlapping-reads` BAMs), the production AMBER
input. They are used only for sex inference (Decision 4). The AMBER sex
call PURPLE reports for these samples is female for 90 and male for 33, in
agreement with the laboratory record wherever it gives a sex.

The manifest is frozen, versioned data, not re-derived. Any change to the
cohort is a new version of the file.

### Decision 6 — production COBALT ratio files as training input

The COBALT ratio files are the production outputs of `eggd_cgp-cobalt`
1.0.2 (COBALT 3.0-beta.5) on the unclipped non-UMI BAMs. That run uses the
normalisation file of the time, but `NormalisationFileBuilder`'s
calculation uses only two columns of each ratio file, `tumorReadDepth` and
`tumorGCContent` (`cobalt.norm.DataLoader`, `Normaliser`); `tumorGCRatio`,
the column the normalisation file changes, is read but not used. COBALT 3.0-beta.5 run without a
normalisation file on the BAMs of two cohort members (one female, one male)
gave a `tumorReadDepth` identical to the production file in all 3,088,257
windows, and a `tumorGCContent` identical in all 74,502 rows the builder
reads.

The ratio and BAF files were extracted from the production `cobalt.tar.gz`
and `amber.tar.gz` outputs, named by specimen, and stored in
`project-J88p7V0470jPZ2Vv1VqB3BFz:/target_regions_normalisation_training/`
(`cobalt/`, `amber/`).

See the Confluence controlled documents for
[v1.1.0](https://cuhbioinformatics.atlassian.net/wiki/spaces/DV/pages/4875714639)
and [v1.0.0](https://cuhbioinformatics.atlassian.net/wiki/spaces/DV/pages/4805623939)
for the fuller narrative behind each of these, including what was tried
and ruled out.

## Usage — recreating the file

```bash
python3 build_target_regions_normalisation.py
```

The output defaults to `target_regions_normalisation_v1.1.0.tsv` in the work
directory.

Requires: `dxpy` (the DNAnexus Python API, authenticated via the same
credentials as the `dx` CLI), `java` (a JRE — no Docker needed,
`NormalisationFileBuilder` is a plain JAR invocation), Python 3 stdlib
otherwise.

Run without `--workdir` for a fresh build — the script creates a private,
freshly-generated temp directory (`tempfile.mkdtemp()`) and prints its path
when it's created, so there's never a stale-artifact risk between runs. Pass an
explicit `--workdir` to reuse a location (required for `--skip-download`,
since there's nothing to reuse otherwise). Either way, this is the exact
same command whether you're recreating the file for the first time or
independently verifying reproducibility — there is no separate test mode.
A real, fresh run downloads the jar, BED, GC profile, and all 246
cohort-level files from scratch, runs the real `NormalisationFileBuilder`,
and the script asserts the result against the pinned checksum
(`EXPECTED_OUTPUT_MD5`/`EXPECTED_OUTPUT_LINES`), exiting non-zero on any
mismatch — so a clean exit *is* the reproducibility proof. v1.1.0 was
confirmed this way by two independent from-scratch runs (about 12 minutes
each: about 6 minutes downloading the 246 cohort files, about 6 minutes in
`NormalisationFileBuilder`).

## Tests

```bash
pip install pytest dxpy
pytest test_build_target_regions_normalisation.py -v
```

Unit tests cover all pure logic (manifest loading, sample-list
construction, checksum verification) with synthetic fixtures — no DNAnexus
or Java access needed. The actual `NormalisationFileBuilder` invocation is
exercised by the reproducibility runs instead, not mocked here.

## Files

| File | Purpose |
|---|---|
| `build_target_regions_normalisation.py` | Main build script |
| `training_cohort_manifest.tsv` | The frozen 123-sample training cohort — SampleId, RatioFileId, AmberBafFileId |
| `test_build_target_regions_normalisation.py` | Unit tests |
