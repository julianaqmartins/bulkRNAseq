# MultiQC Report

This folder contains the aggregated quality control summary from [MultiQC](https://multiqc.info/), combining results from FastQC, STAR, Salmon, and other tools.

## Output

- `multiqc_report.html`: Interactive summary of the pipeline run
- `multiqc_data/`: Raw metrics (JSON, TSV) for each tool
- `multiqc_config.yaml`: MultiQC configuration

---

⚠️ These files are typically too large or browser-sensitive for GitHub, so only examples or summaries may be included.

To open the report:

```bash
open results/multiqc/multiqc_report.html
