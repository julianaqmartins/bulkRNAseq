# STAR + Salmon Outputs

This folder contains the outputs from `nf-core/rnaseq` when run with the `--aligner star_salmon` option. It includes transcript-level and gene-level quantifications.

## Contents

- `salmon.merged.gene_counts.tsv`: Raw estimated gene counts from Salmon (aggregated to gene level).
- `salmon.merged.gene_tpm.tsv`: TPM-normalized expression values (gene level).
- Individual sample folders: Contain sorted BAM files, quality metrics, and intermediate outputs.

## Note

This repository includes **truncated placeholder versions** of the count files (`.tsv`). The full files are too large to upload to GitHub but are generated automatically by the pipeline.

To regenerate the full files, rerun the `nf-core/rnaseq` pipeline with the `--aligner star_salmon` option.
