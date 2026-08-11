# Attribution

This repository distributes derived chromosome-1 sequence and MUC1-scoped reference files
for [VNtyper](https://github.com/hassansaei/VNtyper). Each release specification under
`releases/` records the exact upstream source URL and SHA-256 digest used for its build, so
every published byte is traceable to the release it was derived from.

The MIT licence in this repository covers **this repository's own scripts and metadata**.
It does not extend to the upstream sequence data, which remains subject to the terms below.

## Upstream sources

- **UCSC Genome Browser** — hg19 and hg38 chromosome 1.
  <https://hgdownload.soe.ucsc.edu/> · [Conditions of use](https://genome.ucsc.edu/conditions.html)
  The human genome assembly sequence is freely available for any use.
- **NCBI RefSeq** — GRCh37.p13 and GRCh38.p14 chromosome 1.
  <https://www.ncbi.nlm.nih.gov/refseq/> · [NLM copyright and use](https://www.ncbi.nlm.nih.gov/home/about/policies/)
  NCBI places no restrictions on the use or distribution of its data.
- **Ensembl** — GRCh37 and GRCh38 chromosome 1, pinned to an explicit release, never `current`.
  <https://www.ensembl.org/> · Ensembl data is released under a
  [no-restriction policy](https://www.ensembl.org/info/about/legal/disclaimer.html).

## Derived files

- `muc1_region_hg19.fa` and `muc1_region_hg38.fa` are `samtools faidx` cuts of UCSC chr1
  at the coordinates recorded in each release spec. Both are reproducible byte-for-byte
  from the pinned source; the build asserts this against a committed digest rather than
  trusting it.
- BWA indexes are built with the BWA version recorded in `BUILD_INFO.json`.

## Seeds

`seeds/` holds the artefacts that have no upstream and no derivation script, and so must
be carried:

| file | provenance |
|---|---|
| `MUC1_motifs_Rev_com.fa` | hand-curated MUC1 VNTR motifs |
| `code-adVNTR_RUs.fa` | hand-curated adVNTR repeat units |
| `vntr_db_advntr.zip` | adVNTR databases with all non-MUC1 entries removed, and the MUC1 entry added to the hg38 database with its start position located by UCSC BLAT |
| `filter_config.json` | disallowed motif pairings for `generate_vntr_reference.py` |
| `generate_vntr_reference.py` | derives `All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa` from the two files above |

adVNTR itself: Bakhtiari M *et al.* Targeted genotyping of variable number tandem repeats
with adVNTR. *Genome Research*. 2018;28(11):1709–1719.
