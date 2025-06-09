# DESeq2 Differential Expression Outputs

This folder contains results from DESeq2 differential expression analysis, run on gene-level counts from `nf-core/rnaseq` (STAR + Salmon).

## Contents

Each comparison (e.g. `H1_AAVS1 vs H1_WT`) produces:

- `*_ALL.csv`: All genes with DESeq2 statistics.
- `*_FILTERED.csv`: Significantly differentially expressed genes (FDR < 0.05, |log2FC| > 1).
- `*_volcano.pdf`: Annotated volcano plots for visualization.

### Example Files

- `H1_H1_AAVS1_vs_H1_WT_ALL.csv`  
- `H1_H1_AAVS1_vs_H1_WT_FILTERED.csv`  
- `H1_H1_AAVS1_vs_H1_WT_volcano.pdf`  
- ...and so on for each group-condition pair.

## Note

Only a few example/truncated files may be committed to GitHub for illustration. Full DESeq2 outputs are generated locally using:

```bash
Rscript deseq2/deseq2_gene_names.R \
  results/star_salmon/salmon.merged.gene_counts.tsv \
  deseq2/metadata_stem_1.csv \
  results/deseq2_outputs
