# Publication-quality DESeq2 analysis script using Salmon gene-level counts
# - Input: salmon.merged.gene_counts.tsv (with gene_id, gene_name, and counts)
# - Metadata: CSV with sample_id, condition, group_id
# - Output: Filtered + full results CSVs and volcano plots with labeled gene names

library(DESeq2)
library(tidyverse)
library(ggrepel)
library(EnhancedVolcano)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript deseq2_gene_names.R <gene_counts.tsv> <metadata.csv> <output_dir>")
}

counts_file <- args[1]
metadata_file <- args[2]
output_dir <- args[3]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load gene-level counts
counts_df <- read_tsv(counts_file)
counts_mat <- counts_df %>%
  column_to_rownames("gene_id")

# Preserve gene_name column separately
gene_names <- counts_mat$gene_name
counts_mat <- counts_mat %>% select(-gene_name)

# Convert column names to safe R format
colnames(counts_mat) <- make.names(colnames(counts_mat))

# Load metadata
metadata <- read_csv(metadata_file)
metadata$sample_id <- gsub("-", ".", metadata$sample_id)
metadata$sample_id <- ifelse(grepl("^[0-9]", metadata$sample_id), paste0("X", metadata$sample_id), metadata$sample_id)

# Ensure gene names are unique rownames
rownames(counts_mat) <- make.unique(gene_names)

# Loop through each group (e.g. H1, KOLF, WTC)
for (group in unique(metadata$group_id)) {
  cat("Processing stem cell line:", group, "\n")
  md_group <- metadata %>% filter(group_id == group)
  wt_condition <- unique(md_group$condition[str_detect(md_group$condition, "WT$")])
  
  if (length(wt_condition) != 1) {
    stop(paste("Expected exactly 1 WT condition, found:", length(wt_condition)))
  }

  samples_group <- md_group$sample_id
  counts_group <- counts_mat[, samples_group]
  
  for (cond in unique(md_group$condition)) {
    if (cond == wt_condition) next
    cat("  Comparing:", cond, "vs", wt_condition, "\n")

    samples_to_use <- md_group %>% filter(condition %in% c(cond, wt_condition))
    selected_counts <- counts_group[, samples_to_use$sample_id]

    coldata <- DataFrame(condition = factor(samples_to_use$condition, levels = c(wt_condition, cond)))
    rownames(coldata) <- samples_to_use$sample_id

    dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(selected_counts)), colData = coldata, design = ~condition)
    dds <- DESeq(dds)

    res <- results(dds, contrast = c("condition", cond, wt_condition))
    res_df <- as.data.frame(res) %>% 
      rownames_to_column("gene") %>% 
      arrange(padj)

    # Save results
    prefix <- paste(group, cond, "vs", wt_condition, sep = "_")
    write_csv(res_df, file.path(output_dir, paste0(prefix, "_ALL.csv")))

    res_filtered <- res_df %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
    write_csv(res_filtered, file.path(output_dir, paste0(prefix, "_FILTERED.csv")))

    # Volcano plot
    top_genes <- res_df %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(10) %>% pull(gene)
    pdf(file.path(output_dir, paste0(prefix, "_volcano.pdf")), width = 7, height = 6)
    EnhancedVolcano(res_df,
                    lab = res_df$gene,
                    selectLab = top_genes,
                    x = "log2FoldChange",
                    y = "padj",
                    pCutoff = 0.05,
                    FCcutoff = 1,
                    title = paste(group, ":", cond, "vs", wt_condition),
                    pointSize = 2.0,
                    labSize = 3.0,
                    drawConnectors = TRUE,
                    widthConnectors = 0.5)
    dev.off()
  }
}
