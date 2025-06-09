# Example Figures

This folder includes representative figures generated from the DESeq2 analysis:

## 🔬 Volcano Plot

**File:** `example_volcano_H1_AAVS1_vs_WT.pdf`

- Differential expression analysis between edited (H1_AAVS1 via Cas9-HDR) and wild-type (H1_WT) samples.
- Highlighted genes meet FDR < 0.05 and |log2FC| > 1 thresholds.
- Top 10 most significant genes are labeled.

## 📊 PCA Plot

**File:** `example_PCA_all_samples.pdf`

- Principal Component Analysis (PCA) of all RNA-seq samples.
- Samples grouped by condition.
- Labels correspond to experimental conditions.

---

These are for visualization only. Full analysis outputs are available in `results/deseq2_outputs/`.
