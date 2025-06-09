## 📂 results/

The `results/` directory contains the outputs of the entire bulk RNA-seq pipeline.

- `fastqc/`: Raw FASTQ quality reports
- `multiqc/`: Aggregated quality reports
- `star_salmon/`: STAR alignment and Salmon quantifications (gene-level and transcript-level)
- `deseq2_outputs/`: DESeq2 differential expression results (with plots)

All gene counts and differential expression analyses are reproducible with the included config and R scripts.
