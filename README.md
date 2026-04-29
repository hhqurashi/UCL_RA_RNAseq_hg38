# UCL RA RNA-seq hg38 Analysis

This repository contains R scripts used for downstream bulk RNA-seq analysis using the human `hg38` reference genome.

The workflow covers differential expression analysis, gene set enrichment analysis, and heatmap generation from processed RNA-seq count/expression outputs.

## Repository contents

```text
create_heatmap_updated_270126.R
run_DESeq2.R
run_gsea_fgsea.R
```

## Workflow summary

The scripts are used to:

1. Load processed RNA-seq count/expression data
2. Run differential expression analysis with DESeq2
3. Generate ranked gene lists from differential expression results
4. Run gene set enrichment analysis using fgsea
5. Create heatmaps for selected genes or significant results
6. Export plots and tables for downstream interpretation and reporting

## Script summary

### `run_DESeq2.R`

Runs differential expression analysis using DESeq2.

Expected outputs include:

- differential expression result tables
- normalised count tables
- significant gene lists
- diagnostic plots such as PCA/MA-style outputs

### `run_gsea_fgsea.R`

Runs gene set enrichment analysis using `fgsea`.

Expected outputs include:

- ranked gene lists
- enriched pathway tables
- enrichment plots or summary outputs

### `create_heatmap_updated_270126.R`

Generates heatmaps from processed RNA-seq expression or differential expression outputs.

Expected outputs include:

- heatmap figures
- filtered expression matrices
- visual summaries of selected genes/pathways

## Input

Expected inputs include:

- gene-level count matrices
- sample metadata
- DESeq2 result tables
- gene set files
- processed expression matrices

Raw sequencing data and large upstream pipeline outputs are not tracked in this repository.

## Output

Outputs are written to the configured results directories and include:

- CSV/TSV result tables
- normalised expression matrices
- enrichment results
- heatmap figures
- downstream plots for reporting

Large results files should remain outside GitHub unless they are small, non-sensitive, and useful as examples.
