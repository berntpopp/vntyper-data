# VNtyper reference data releases

Immutable metadata and GitHub Release assets for [VNtyper](https://github.com/hassansaei/VNtyper)
MUC1 reference bundles. **No reference data is committed here** — no chromosome FASTAs, no
BWA indexes, no generated motif files. Git holds the release specifications, the
non-derivable seeds, and the workflow that builds a release.

## Why this repository exists

VNtyper used to assemble its references at Docker image build time: six chromosome-1
genomes downloaded from UCSC, NCBI and Ensembl and BWA-indexed in-image. That was the
single most expensive step of the build, it depended on four third-party hosts staying up
and byte-stable, and nothing pinned the *content* of what was downloaded — two images built
a year apart could contain different sequence under the same name.

Reference sets are now built once, verified, and published here as versioned release
bundles. `vntyper install-references` and `Dockerfile.base` fetch them. See
[hassansaei/VNtyper#217](https://github.com/hassansaei/VNtyper/issues/217).

## Release lifecycle

VNtyper software is released from
[hassansaei/VNtyper](https://github.com/hassansaei/VNtyper). Reference data is released
here with tags of the form `refs-vN`.

Each committed JSON file under `releases/` pins:

- every upstream source URL, with its SHA-256 digest, and Ensembl pinned to an explicit
  release rather than `current`;
- the SHA-256 of every seed in `seeds/`;
- every derivation, as the command that produces it and the SHA-256 it must produce;
- the VNtyper commit the builder was taken from, and the BWA and samtools versions used.

The release workflow invokes a full-SHA-pinned reusable workflow from the software
repository. It creates a **draft** release only. Assets are verified before a maintainer
publishes the draft, and published releases are immutable — a builder change means a new
`refs-vN`, never a re-cut of an existing tag.

## Install and verify

VNtyper resolves this repository by default:

```bash
vntyper install-references --output-dir reference --references hg38
```

Each asset's SHA-256 is committed in VNtyper's own `install_references_config.json`, not
taken from the `SHA256SUMS` published beside the assets — a checksum file hosted next to
the files it describes cannot be its own root of trust. The published `SHA256SUMS` is for
humans reviewing a draft:

```bash
gh release download refs-v1 --repo berntpopp/vntyper-data
sha256sum --check --ignore-missing SHA256SUMS
```

To build from upstream instead of fetching a bundle — the path this repository's own build
workflow runs, so it stays exercised:

```bash
vntyper install-references --output-dir reference --from-source
```

## Publication policy

Every release must contain, per assembly, a bundle carrying the chromosome FASTA, its BWA
index and a `BUILD_INFO.json`; the common MUC1 bundle; a top-level `SHA256SUMS`;
`release-manifest.json`; and `verification-report.json`. Publish only after the full
artifact set and its remote checksums have been reviewed.

## Status

**`refs-v2`** is the current release, consumed by
[VNtyper v2.0.21](https://github.com/hassansaei/VNtyper/releases/tag/v2.0.21) and later.
It corrects the GRCh38 adVNTR model, which described 840 bp of a repeat array GRCh38
carries at 3,525 bp — adVNTR derives its read-fetch window from the model's own content,
so it saw 24% of the locus. Measured over 400 simulated samples, detection goes from
132/200 carriers to 176/200. **Requires adVNTR >= 2.0.4.**

Only the MUC1 common asset changed. The six per-assembly bundles are the `refs-v1`
artifacts byte-for-byte, republished under the new tag keeping their original file names
and digests, so their own `BUILD_INFO.json` and manifests stay truthful about when and
from what they were built. `verification-report.json` records that check.

The models are now derived rather than shipped as an opaque blob:
`seeds/derive_advntr_muc1_model.py` builds them from chr1 plus pinned array bounds, and
reproduces the previous hg19 database content byte-for-byte — 9 differing bytes, all
SQLite header, zero content bytes. See
[hassansaei/VNtyper#268](https://github.com/hassansaei/VNtyper/issues/268) and #1.

`refs-v2` was cut published rather than as a draft, skipping the review step described
above. That deviation and the verification available for post-hoc review are recorded in
the release notes. `refs-v1` remains published and untouched, so installations pinned to
it are unaffected.

`refs-v1` is specified in `releases/refs-v1.json` and built by
`.github/workflows/release-data.yml`, which dispatches VNtyper's reusable builder pinned at
`699b81e` — the merge commit of
[hassansaei/VNtyper#239](https://github.com/hassansaei/VNtyper/pull/239) and #240.

The release is cut as a **draft**. A maintainer reviews `verification-report.json` and checks
`SHA256SUMS` against the remote assets before publishing. Published releases are immutable:
a correction is a new tag, never a re-cut of an existing one.

**The trust anchor is VNtyper, not this repository.** Every source URL and digest in a spec
here must match the values committed in VNtyper's `install_references_config.json`, and the
build refuses to start if they disagree. To move an upstream, change VNtyper first in a
reviewed commit, then match the spec to it.
