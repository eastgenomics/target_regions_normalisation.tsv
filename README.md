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

## Why chr-prefix the BED, and why `-ref_genome_version 38`?

The backbone BED is sourced as no-chr
(`backbone_padded_150bp_nochr.bed`), but `NormalisationFileBuilder`'s
chromosome-matching for the GC profile is a raw, unnormalised string
lookup — not chr-prefix agnostic — so it must match the GC profile's own
convention (chr-prefixed) exactly, or it crashes. Separately, the file
writer reconstructs each chromosome key from `-ref_genome_version`'s own
naming convention (hmftools' standard: V37 → no-chr, V38 → chr-prefixed) and
silently produces zero output rows if that doesn't match the BED's
convention. `-ref_genome_version 38` is also the only value the approved
production jar (COBALT 3.0-beta.5) accepts at all — `37` throws an
`IndexOutOfBoundsException` in this build. See the `DECISIONS` block in
`build_target_regions_normalisation.py` for the full evidence trail,
including a full 41-sample controlled A/B rebuild confirming this produces
identical output to the historical file (once the historical build's own
different, unapproved jar version is held constant) — the only real
difference between the two approaches is which jar version is used, not
the chr-prefix/ref-genome-version fix itself.

## Usage — recreating the file

```bash
python3 build_target_regions_normalisation.py --out target_regions_normalisation.tsv
```

Requires: `dx` (authenticated DNAnexus CLI), `java` (a JRE — no Docker
needed, `NormalisationFileBuilder` is a plain JAR invocation), Python 3
stdlib only otherwise.

Run without `--workdir` for a fresh build — the script creates a private,
freshly-generated temp directory (`tempfile.mkdtemp()`) and prints its path
on completion, so there's never a stale-artifact risk between runs. Pass an
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
