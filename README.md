# Bulk RNA-seq Pipeline (nf-core/rnaseq)

This repository documents a working configuration of the `nf-core/rnaseq` pipeline (v3.12.0) using Docker containers and custom genome references.

## ✅ Requirements
- Docker
- Nextflow >= 21.04.0
- Google Cloud SDK (`gsutil`)
- ~500 GB free space on main disk or a mounted external volume

## ✅ Genome
Genome: GRCh38 (GENCODE v46)

Folder: `~/bulk_rnaseq_clean_run/genome/gencode_v46/`

## ✅ Pipeline Command

```bash
nextflow run nf-core/rnaseq -r 3.12.0 -profile docker \
  --input input/stem_sample_sheet.csv \
  --outdir results \
  --fasta ~/bulk_rnaseq_clean_run/genome/gencode_v46/GRCh38.primary_assembly.genome.fa \
  --gtf ~/bulk_rnaseq_clean_run/genome/gencode_v46/gencode.v46.annotation.gtf \
  --star_index ~/bulk_rnaseq_clean_run/genome/gencode_v46/STARIndex \
  -c config/docker_user.config \
  -with-report -with-timeline -with-trace# bulkRNAseq
