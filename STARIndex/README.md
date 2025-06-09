# STARIndex

This folder is expected to contain the pre-built STAR genome index based on GENCODE v46 (GRCh38).
⚠️ Due to the large size of STAR index files, the contents are **not included** in this repository.

## 🔧 How to Generate the STAR Index

Use the STAR aligner to generate the index manually before running the pipeline.

### Requirements

- STAR 2.7.10b (quay.io/biocontainers/star:2.7.10b--h9ee0642_0)
- GENCODE v46 FASTA and GTF files

### Command

```bash
STAR --runMode genomeGenerate \ #Run STAR in genome indexing mode.
  --genomeDir genome/gencode_v46/STARIndex \ #Output directory for STAR index.
  --genomeFastaFiles genome/gencode_v46/GRCh38.primary_assembly.genome.fa \ #Your GRCh38 reference genome in FASTA format.
  --sjdbGTFfile genome/gencode_v46/gencode.v46.annotation.gtf \ #Your GTF annotation file.
  --sjdbOverhang 149 \ #sjdbOverhang should be set to (read length - 1). For 150 bp reads, use 149.
  --runThreadN 16 #Adjust --runThreadN based on your available CPU cores.
```
