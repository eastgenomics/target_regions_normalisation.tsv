# target_regions_normalisation.tsv

Reproducible build script for `target_regions_normalisation.tsv` — the COBALT
panel-normalisation file used by `eggd_cgp-cobalt` to correct per-window
read-depth ratios for CGP backbone target-region biases.

## What this builds

`target_regions_normalisation.tsv` is built once, from a fixed 41-sample
training cohort ("EF v1"), by running Hartwig's own
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
   3.0-beta.5), the backbone target-regions BED, and the GC profile from
   DNAnexus, plus the 41-sample training cohort's COBALT ratio and AMBER BAF
   files (file IDs frozen in `training_cohort_manifest.tsv`).
2. chr-prefixes the backbone BED (it's sourced as no-chr) — see *Why
   chr-prefix?* below.
3. Builds a `SampleId`-only sample list — no Gender column, so gender is
   inferred by `NormalisationFileBuilder` itself from real AMBER BAF data.
4. Runs `NormalisationFileBuilder` and asserts the result against a pinned
   checksum.

Full narrative (why each decision was made, what was tried, the historical
build's actual mechanism and its problems) is in the Confluence controlled
document:
[target_regions_normalisation.tsv](https://cuhbioinformatics.atlassian.net/wiki/spaces/DV/pages/4805623939).

## Decisions

Every manual, judgement-based resolution from the investigation that
produced this script, with the full rationale and evidence. The
`DECISIONS` block in `build_target_regions_normalisation.py` records the
same five items in short form, pointing back here.

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

### Decision 3 — chr-prefix the backbone BED ourselves

The backbone BED is chr-prefixed here (from the live no-chr source,
`backbone_padded_150bp_nochr.bed`) rather than reusing a separate no-chr GC
profile, as the historical build did. Two reasons:

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
  lines) -- so there was never a need for two separate copies. Simplest
  correct fix: chr-prefix the BED, keep the one GC profile already in
  `resource_ids.env`.

### Decision 4 — no Gender column in the sample manifest

Gender for all 41 samples is inferred by `NormalisationFileBuilder` itself
from real AMBER BAF chrX heterozygosity data
(`com.hartwig.hmftools.common.amber.AmberGender`), using whichever
pseudo-autosomal-region (PAR) boundary set matches `-ref_genome_version`.
The historical build (V37 PAR boundaries against genuinely
GRCh38-aligned AMBER BAF data) was a real, if theoretical, coordinate
mismatch -- confirmed via a full 41-sample controlled A/B rerun to have
changed ZERO gender calls and ZERO output values for this cohort (holding
the jar version constant). V38 is still the correct choice because it's
genuinely consistent, not because the V37 mismatch caused visible harm
here.

### Decision 5 — the 41-sample cohort manifest is frozen

This cohort manifest (41 EF v1 samples) is frozen, versioned data, not
re-derived. It matches the samples used for the currently-live production
file, letting a rebuild be checked against that file directly. A
different/larger training cohort is a deliberate, separate decision (new
panel validation), not something this script infers on its own.

See the [Confluence controlled
document](https://cuhbioinformatics.atlassian.net/wiki/spaces/DV/pages/4805623939)
for the fuller narrative behind each of these, including what was tried
and ruled out.

## Usage — recreating the file

```bash
python3 build_target_regions_normalisation.py --out target_regions_normalisation.tsv
```

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
A real, fresh run downloads the jar, BED, GC profile, and all 82
cohort-level files from scratch, runs the real `NormalisationFileBuilder`,
and the script asserts the result against the pinned checksum
(`EXPECTED_OUTPUT_MD5`/`EXPECTED_OUTPUT_LINES`), exiting non-zero on any
mismatch — so a clean exit *is* the reproducibility proof. Confirmed this
way across multiple independent from-scratch runs during development
(~4.5 minutes each, dominated by downloading 82 cohort files).

## Tests

```bash
pip install pytest
pytest test_build_target_regions_normalisation.py -v
```

Unit tests cover all pure logic (chr-prefixing, manifest loading, sample-list
construction, checksum verification) with synthetic fixtures — no DNAnexus
or Java access needed. The actual `NormalisationFileBuilder` invocation is
exercised by the reproducibility runs instead, not mocked here.

## Files

| File | Purpose |
|---|---|
| `build_target_regions_normalisation.py` | Main build script |
| `training_cohort_manifest.tsv` | The frozen 41-sample EF v1 training cohort — SampleId, RatioFileId, AmberBafFileId |
| `test_build_target_regions_normalisation.py` | Unit tests |
