# Standard MUC1 VNTR Length Estimation Model (GRCh38)

## Overview

This directory contains the standard 13-feature Bayesian Ridge model (`grch38-standard-length-model-v1.json`) used by VNtyper for mapping-based MUC1 VNTR repeat count estimation from short-read sequencing data.

- **Model File**: `grch38-standard-length-model-v1.json`
- **Model Version**: `standard13-bayesian-v1`
- **Assembly**: `GRCh38` (contigs `chr1`, `1`)
- **Target**: Complete diploid repeat count (`complete` count convention)
- **Model Digest (SHA-256)**: `0d910946ab14d38cf5b3bd33a0b16e1b985c8377a9d2b466d98693a324ca3b51`

---

## Invariable Repeat Units & Count Convention

The MUC1 coding VNTR region contains terminal invariable units flanking the central variable canonical repeats:
- **Pre-repeats (Units 1–5)**: 5 invariable units (300 bp, `chr1:155191940–155192239` on GRCh38).
- **After-repeats (Units 6–9)**: 4 invariable units (240 bp, `chr1:155188487–155188726` on GRCh38).
- **Total Invariable Units**: $5 + 4 = 9$ units per allele ($18$ diploid units).

Ground truth repeat counts obtained by PacBio long-read circular consensus sequencing (Prague laboratory, Martina Živná / Stanislav Kmoch) report:
$$\text{Length of allele} = \text{\# of pre-repeats} + \text{\# of repeats} + \text{\# of after-repeats}$$
Therefore, reported ground truth counts **already include** the 9 invariable units per allele. This model predicts the **complete** diploid repeat count ($Y_{\text{complete}} = \text{short} + \text{long}$).
*(To obtain canonical core repeats alone, subtract 18: $Y_{\text{canonical}} = Y_{\text{complete}} - 18$.)*

---

## Training Dataset & Specifications

- **Sample Count**: $N = 76$ independent clinical exome specimens with confirmed PacBio paired numeric truth.
  - German cohort: $N = 42$ samples.
  - French cohort: $N = 34$ samples.
- **Enrichment / Library Prep**: Twist Comprehensive Exome capture.
- **Sequencing Technology**: Illumina paired-end short-read sequencing (100–150 bp).
- **Ground Truth Technology**: PacBio Sequel single-molecule real-time (SMRT) circular consensus sequencing (CCS), with consensus reconstruction of both alleles.
- **Locus Definition**: GRCh38 `chr1:155188296–155192429` (4,133 bp spanning flanks, invariable units, and central core array).
- **Quality Control**: 100% of training samples pass frozen QC thresholds (coverage $> 90\%$, invariant mean depth $> 10\times$, supporting fragments $> 100$, non-zero reads).

---

## Statistical Methodology

The model extracts 13 closed, reproducible alignment and pileup features from the indexed locus:
1. `A`: Depth ratio of core variable repeats to terminal invariant units.
2. `F`: Depth ratio of the VNTR array to 190 bp flanks.
3. `core_zero_fraction`: Fraction of zero-depth positions in the core array.
4. `core_depth_cv`: Coefficient of variation of per-base depth across core repeats.
5. `core_bin_cv`: Coefficient of variation across 8 equal bins of the core array.
6. `core_half_log_ratio`: Log ratio of first-half to second-half core depth.
7. `invariant_end_log_ratio`: Log ratio of 5' to 3' invariant depth.
8. `flank_end_log_ratio`: Log ratio of left to right flank depth.
9. `log_invariant_depth`: $\ln(1 + \text{depth})$ of invariant terminal units.
10. `mapq_zero_fraction`: Fraction of locus reads with $\text{MAPQ} = 0$.
11. `mean_mapq`: Mean mapping quality across locus reads.
12. `soft_clipped_read_fraction`: Fraction of reads with soft-clipped CIGAR operations.
13. `query_sequence_gc_fraction`: GC base fraction of reads at the locus.

Parameters are fitted via Bayesian Ridge regression (`max_iter=1000`) with standardized features, then collapsed into exact raw-space coefficients and intercept.

---

## Validation & Benchmark Results

### 1. Leave-One-Out Cross-Validation ($N = 76$)

| Model | MAE (repeats) | Median AE | RMSE | Bias | $R^2$ | Within $\pm 20\%$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Training Mean Baseline** | 19.24 | 17.73 | 24.72 | 0.00 | -0.027 | 72.4% |
| **Linear Ratio $A$ Alone** | 18.17 | 15.20 | 22.06 | -0.09 | 0.182 | 73.7% |
| **Standard 13-Feature Bayesian Ridge** | **10.80** | **8.62** | **14.26** | **+0.07** | **0.659** | **88.2%** |

### 2. Cross-Cohort Transfer (German $\to$ French)
When trained strictly on the German cohort ($N=42$) and evaluated on the held-out French cohort ($N=34$):
- **MAE**: **10.12 repeats**
- **Median AE**: **7.51 repeats**
- **RMSE**: **13.49 repeats**
- **Bias**: **+2.77 repeats**
- **$R^2$**: **0.712**
- **Agreement within $\pm 20\%$**: **91.2%**
