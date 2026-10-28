# mirna-tf-mrna-cancer

**Transcription factors mediate microRNA-gene relationships in cancer**

MicroRNAs (miRNAs) regulate gene expression by binding to messenger RNAs (mRNAs), and the usual way to find miRNA targets is to look for miRNA–mRNA expression correlations. But miRNAs also target transcription factors (TFs), so many of those correlations may be *indirect*: the miRNA regulates a TF, and the TF regulates the gene.

This project tests that idea using matched mRNA and miRNA sequencing data from The Cancer Genome Atlas (TCGA). For each miRNA it:

1. finds the genes whose expression correlates significantly with the miRNA (Spearman, BH-adjusted *p* < 0.05);
2. finds TFs the miRNA is experimentally validated to target ([multiMiR](https://doi.org/10.1093/nar/gku631), [Human TFs](http://humantfs.ccbr.utoronto.ca));
3. finds each TF's target genes from cancer cell-line ChIP-seq ([ChIP-Atlas](https://chip-atlas.org));
4. removes the miRNA's own direct targets and uses over-representation analysis ([clusterProfiler](https://doi.org/10.1016/j.xinn.2021.100141)) to test whether the remaining correlated genes are enriched for that TF's targets.

Significant miRNA–TF pairs point to candidate miRNA → TF → mRNA regulatory relationships.

## History

This started as my Master's thesis at Monash University (submitted in June 2024) with [Dr. Miles C. Andrews](https://www.monash.edu/medicine/translational/cancer-medicine/research/andrews-group) as my supervisor and [Dr. Peinan Zhao](https://scholar.google.com/citations?user=3t8PeMwAAAAJ&hl=en) as my co-supervisor. As a part of the thesis, we analysed melanoma only (**TCGA-SKCM**). We then extended the analysis to breast cancer (**TCGA-BRCA**) to test whether these indirect, TF-mediated miRNA–mRNA relationships generalise beyond melanoma.

We developed an approach for inferring indirect, TF-mediated relationships between miRNAs and genes, and identified candidate miRNA–TF–gene triplets that warrant experimental validation. 


## Running the analysis

Open `R/main.R`, set the options at the top and source it:

```r
CANCER     <- "BRCA"          # "BRCA", "SKCM", or c("BRCA", "SKCM")
REGENERATE <- character(0)    # read saved outputs; "all" or e.g. c("ora") to recompute
```

Study-specific settings are in `R/config_brca.R` and `R/config_skcm.R`; shared settings are in `R/config_common.R`. Saved intermediate results are read from `data_files/<study>/` unless a step is listed in `REGENERATE`.


**Data:** the raw TCGA files (`*.RDS`), filtered expression matrices, correlation tables and ChIP-Atlas TF target files are large and not tracked in git. They can be regenerated with `REGENERATE <- c("tcga_data", "chip_atlas", "all")` (this downloads several GB).


## Poster

I presented this work as a poster at the [ABACBS 2025 conference](https://static1.squarespace.com/static/5423875be4b03f0c482a58c4/t/692282310b37fa4406e88c30/1763869233702/ABACBS+FINAL+program.pdf) in Adelaide. The full-resolution poster is available as a [PDF](poster/ABACBS_Poster.pdf).

![Poster presented at ABACBS 2025, Adelaide](poster/abacbs_poster.JPG)
