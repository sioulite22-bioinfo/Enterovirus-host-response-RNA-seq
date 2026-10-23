# Host blood transcriptional response to entero/rhinovirus infection: a re-analysis with cell-composition checks

A reproducible re-analysis of public whole-blood RNA-seq data (GEO: [GSE157240](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157240)) asking whether entero/rhinovirus infection shows evidence of **host metabolic reprogramming**, and how much of any apparent metabolic signal survives controls for batch and blood cell composition.

> **Status:** exploratory analysis of public data by a single author. This is a re-analysis, not original research, and it makes no novel or validated biological claims. All findings below are hypotheses limited by the design described under *Limitations*.

## Summary of findings

1. **A strong interferon response** (Hallmark interferon alpha and gamma, NES about 3.1 to 3.2 versus healthy controls) in entero/rhinovirus patients. It persists after adjusting for a blood cell-composition score. It is much **weaker than in influenza**, in both the all-batch analysis and within a single batch.
2. **No robust change in core energy metabolism.** Glycolysis, fatty acid metabolism and oxidative phosphorylation were not consistently changed.
3. **Several apparent metabolic signals did not hold up on inspection:**
   - *Heme metabolism* was strongly enriched, but the leading-edge genes include erythroid-lineage genes (for example GYPB, ANK1, CA1, ERMAP), which suggests immature red-cell transcripts, not immune-cell metabolism. This interpretation was not tested further.
   - *"Selenoamino acid metabolism"* (Reactome) was strongly down, but **96% of its leading-edge genes were RPL/RPS ribosomal protein genes**, so it reports a ribosomal/translation signal, not amino-acid metabolism.
   - *Phospholipid and sphingolipid metabolism* were enriched before, and **not after**, adjustment for cell composition.
   - A respiratory electron transport/OXPHOS difference between entero/rhinovirus and influenza in the all-batch analysis **did not replicate** within a single batch.
4. **Cell composition dominates the signal.** Neutrophil/myeloid markers (CD177, ELANE, S100A8, MPO) were up and lymphocyte markers (CD3E, CD3D, MS4A1, CD79A) were down in infected patients.

## Data

| Item | Detail |
|---|---|
| Source | GEO series GSE157240 (Tsalik et al., *Front Immunol* 2021;12:741837, PMID 34777354) |
| Material | Peripheral blood RNA-seq, Illumina HiSeq 4000 |
| Counts used | NCBI-generated raw counts, GRCh38.p13 (`GSE157240_raw_counts_GRCh38.p13_NCBI.tsv.gz`) |
| Annotation | `Human.GRCh38.p13.annot.tsv.gz` (Annotation Release 109.20190905) |
| Metadata | Series matrix (`GSE157240_series_matrix.txt.gz`) |

The GEO page lists 187 samples; the NCBI count matrix contains **142**, and only these were analysed.

| Group | n |
|---|---|
| Healthy controls | 13 |
| Entero/Rhinovirus | 20 |
| Influenza | 51 |
| Dengue | 14 |
| Metapneumovirus | 13 |
| DNA viruses (adenovirus, CMV, EBV, HSV) | 15 |
| Other respiratory RNA viruses (parainfluenza, RSV) | 16 |

Raw data are **not** included in this repository. Download the three files from the GEO page above and place them in `data/`.

## Methods

1. **Batch inspection.** Sample titles carry a prefix (DU09, DU10, DU14, DU16) that appears to mark cohorts or batches. PCA of variance-stabilised counts (DESeq2 `vst`) was used to look for batch structure. All 13 healthy controls are in DU09, so entero/rhinovirus versus healthy was restricted to **DU09 (6 versus 13 samples)**.
2. **Differential expression.** DESeq2 with `~ sex + group`. Genes were pre-filtered (at least 10 counts in at least 5 samples). Entero/rhinovirus versus influenza was run on all batches with `~ batch + sex + group`, and repeated within DU09 only.
3. **Gene-set enrichment.** `fgsea` on a ranked list of Wald statistics (protein-coding genes, one value per gene symbol), using MSigDB Hallmark sets and Reactome sets whose names match metabolism-related terms (131 sets).
4. **Cell-composition check.** Marker genes for neutrophil/myeloid, lymphocyte and erythroid cells were inspected. A composition score (mean z-scored neutrophil/myeloid markers minus mean z-scored lymphocyte markers) was added as a covariate and the enrichment re-run. This is a crude adjustment, not a proper deconvolution.
5. **Leading-edge inspection** of pathways of interest, to check what the enrichment is actually driven by.

## Key results

| Pathway | EV/RV vs Healthy (DU09) | Flu vs Healthy (DU09) | EV/RV vs Flu (all batches) | EV/RV vs Flu (DU09 only) |
|---|---|---|---|---|
| Interferon alpha response | NES +3.2 \* | NES +3.5 \* | NES -3.6 \* | NES -3.7 \* |
| Interferon gamma response | NES +3.1 \* | NES +3.4 \* | NES -3.6 \* | NES -3.5 \* |
| Heme metabolism | NES +2.7 \* | NES +2.3 \* | NES +1.4 \* | NES +2.1 \* |
| Glycolysis | +1.3 (ns) | +1.3 (ns) | +1.1 (ns) | +1.4 (p = 0.049) |
| Fatty acid metabolism | +1.3 (ns) | +1.5 \* | -1.1 (ns) | +1.2 (ns) |
| Oxidative phosphorylation | +1.0 (ns) | +1.2 (ns) | -1.6 \* | +0.8 (ns) |
| Phospholipid metabolism (unadjusted) | +1.6 \* | +1.1 (ns) | +1.7 \* | +1.8 \* |

\* adjusted p < 0.05. Full table: `results/summary_key_pathways.csv`.

After composition adjustment (DU09), phospholipid and sphingolipid metabolism were no longer enriched (adjusted p between 0.60 and 0.94 in both comparisons), while interferon gamma response remained strongly enriched (adjusted p below 1e-36).

Figures (in `results/figures/`):
- `volcano_EVRV_vs_Healthy.png`: differential expression, entero/rhinovirus versus healthy (DU09).
- `fgsea_selected_hallmark.png`: selected Hallmark pathways (unadjusted).
- `summary_pathways_heatmap.png`: key pathways across the four comparisons. **Unadjusted for composition.**

## Limitations

- **Whole blood is a mixture of cell types.** Infection changes cell proportions, so many pathway changes may reflect composition. The marker-based adjustment is crude and, with these small groups, is strongly correlated with infection status, so it may also remove real signal.
- **Small and unbalanced groups.** The entero/rhinovirus versus healthy comparison has 6 versus 13 samples, with only one male among the 6 (sex is in the model).
- **"Entero/Rhinovirus" is a single label** and cannot be resolved into enterovirus species or genotypes.
- **Batch is partly confounded with group.** Healthy controls come only from DU09, and the entero/rhinovirus samples are concentrated in DU14 and DU09. The prefix-to-batch interpretation is inferred from sample titles and is not documented in the GEO metadata. Site or country information was not available to check further.
- **RNA levels are not metabolite levels.** No metabolomic or flux data were analysed, and no genome-scale metabolic modelling was performed. Pathway enrichment is not a flux prediction.
- **Multiple testing and exploratory scope.** Many comparisons were run, and results should be treated as hypothesis-generating.
- 45 of the 187 listed samples are absent from the NCBI count matrix and were not analysed.

## Reproducing the analysis

Requirements: R (4.5.2 was used) with the packages `tidyverse`, `GEOquery`, `DESeq2`, `fgsea`, `msigdbr`.

```r
install.packages(c("BiocManager", "tidyverse"))
BiocManager::install(c("GEOquery", "DESeq2", "fgsea", "msigdbr"))
```

1. Download the three GEO files listed under *Data* into `data/`.
2. Set `base_dir` at the top of `scripts/01_analysis.R` if your project folder is not auto-detected.
3. Run `scripts/01_analysis.R`. Outputs are written to `results/` and `results/figures/`.

## Repository layout

```
.
├── README.md
├── scripts/
│   └── 01_analysis.R
├── data/               # not tracked; download from GEO
└── results/
    ├── figures/
    └── *.csv

## Results and Visualizations

### Differential expression

![Volcano plot](results/figures/volcano_EVRV_vs_Healthy.png)

### Hallmark pathway enrichment

![Hallmark enrichment](results/figures/fgsea_selected_hallmark.png)

### Pathway comparison

![Pathway enrichment heatmap](results/figures/summary_pathways_heatmap.png)

```

## Reference

Tsalik EL, Fiorino C, Aqeel A, Liu Y, et al. The Host Response to Viral Infections Reveals Common and Virus-Specific Signatures in the Peripheral Blood. *Front Immunol* 2021;12:741837. PMID: 34777354.

## Author

Mohammed Barkathullah S, [github.com/sioulite22-bioinfo](https://github.com/sioulite22-bioinfo)

## Context

This project was done as independent practice in host-response transcriptomics and as an exercise in checking whether pathway-level metabolic signals are robust to batch and cell-composition effects. Extending it to cell-type-resolved data (for example single-cell or organoid infection data) and to genome-scale metabolic modelling is the natural next step.
