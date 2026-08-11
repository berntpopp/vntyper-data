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

`refs-v1` is not yet published. It is blocked on the builder workflow landing in VNtyper —
see the design and plan in that repository under `docs/superpowers/`.
