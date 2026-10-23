# 01_analysis.R
# Whole-blood RNA-seq (GSE157240): entero/rhinovirus vs healthy and vs influenza,
# with batch and cell-composition checks.
# Run top to bottom. Expects the three GEO files in <base_dir>/data/.

suppressPackageStartupMessages({
  library(tidyverse)
  library(GEOquery)
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
})

# ---------------------------------------------------------------- 0. paths ----
base_dir <- if (dir.exists("D:/ev_project")) "D:/ev_project" else "/mnt/d/ev_project"
data_dir <- file.path(base_dir, "data")
res_dir  <- file.path(base_dir, "results")
fig_dir  <- file.path(res_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------ 1. load data ----
counts <- read.delim(file.path(data_dir, "GSE157240_raw_counts_GRCh38.p13_NCBI.tsv.gz"),
                     row.names = 1, check.names = FALSE)
annot  <- read.delim(file.path(data_dir, "Human.GRCh38.p13.annot.tsv.gz"))
gse    <- getGEO(filename = file.path(data_dir, "GSE157240_series_matrix.txt.gz"))
meta   <- pData(gse)

# keep only samples that have counts (142 of the 187 listed on GEO)
meta142       <- meta[colnames(counts), ]
meta142$group <- meta142$`status:ch1`
meta142$sex   <- factor(meta142$`Sex:ch1`)
meta142$batch <- substr(meta142$title, 1, 4)   # DU09 / DU10 / DU14 / DU16 (inferred from titles)

lookup <- c(
  "healthy_ctrl" = "Healthy",
  "Entero/Rhinovirus infected patient" = "EV_RV",
  "Influenza virus infected patient" = "Flu",
  "Dengue virus infected patient" = "Dengue",
  "Metapneumovirus infected patient" = "hMPV",
  "DNA virus infected patient (Adenovirus, Cytomegalovirus, Ebstein-Barr virus, Herpes Simplex virus)" = "DNA_virus",
  "other respiratory RNA virus infected patient (Parainfluenza virus and Respiratory Syncytial virus)" = "Other_RNA"
)
meta142$grp <- unname(lookup[meta142$group])

print(table(meta142$grp))
print(table(meta142$batch, meta142$grp))

# ---------------------------------------------------------- 2. helpers -------
run_fgsea <- function(pc_df, sets) {
  rk <- pc_df %>% group_by(Symbol) %>%
    slice_max(abs(stat), n = 1, with_ties = FALSE) %>% ungroup()
  v <- setNames(rk$stat, rk$Symbol)
  set.seed(1)
  out <- fgsea(pathways = sets, stats = v, minSize = 15, maxSize = 500)
  out[order(out$padj), ]
}

# DESeq2 on a subset of samples; returns protein-coding results with symbols
do_de <- function(meta_sub, lev, design_formula) {
  meta_sub <- droplevels(meta_sub)
  meta_sub$grp <- factor(meta_sub$grp, levels = lev)
  d <- DESeqDataSetFromMatrix(as.matrix(counts[, rownames(meta_sub)]),
                              meta_sub, design = design_formula)
  d <- d[rowSums(counts(d) >= 10) >= 5, ]
  d <- DESeq(d)
  r <- results(d, contrast = c("grp", lev[2], lev[1]))
  as.data.frame(r) %>%
    rownames_to_column("GeneID") %>%
    mutate(Symbol   = annot$Symbol[match(GeneID, annot$GeneID)],
           GeneType = annot$GeneType[match(GeneID, annot$GeneID)]) %>%
    filter(GeneType == "protein-coding", !is.na(padj), !is.na(Symbol), Symbol != "")
}

# ----------------------------------------------------------- 3. gene sets ----
h <- tryCatch(msigdbr(species = "Homo sapiens", collection = "H"),
              error = function(e) msigdbr(species = "Homo sapiens", category = "H"))
h_list <- split(h$gene_symbol, h$gs_name)

rx <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"),
  error = function(e) msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
)
rx_list  <- split(rx$gene_symbol, rx$gs_name)
metab_rx <- rx_list[grepl("METABOLISM|GLYCOLYSIS|CITRIC|RESPIRATORY_ELECTRON|BETA_OXIDATION|FATTY_ACID|LIPID|CHOLESTEROL|GLUTAM|NUCLEOTIDE|PENTOSE",
                          names(rx_list))]
cat("Reactome metabolism-related sets used:", length(metab_rx), "\n")
sets_all <- c(h_list, metab_rx)

# --------------------------------------------- 4. batch / group PCA (all) ----
dds_all <- DESeqDataSetFromMatrix(as.matrix(counts), meta142, design = ~ 1)
dds_all <- dds_all[rowSums(counts(dds_all) >= 10) >= 10, ]
vsd_all <- vst(dds_all, blind = TRUE)

p_pca_batch <- plotPCA(vsd_all, intgroup = "batch") + ggtitle("PCA by batch")
p_pca_grp   <- plotPCA(vsd_all, intgroup = "grp")   + ggtitle("PCA by infection group")
ggsave(file.path(fig_dir, "pca_batch.png"), p_pca_batch, width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "pca_group.png"), p_pca_grp,   width = 6, height = 5, dpi = 300)

# ------------------------- 5. composition markers and composition score ------
neut    <- c("CD177", "S100A8", "S100A9", "MPO", "ELANE")   # neutrophil / myeloid
lymph   <- c("CD3E", "CD3D", "MS4A1", "CD79A")              # T and B cells
markers <- c(neut, lymph, "NKG7", "HBB", "HBA1", "ALAS2")   # + NK and erythroid

expr_all <- assay(vsd_all)
sym_all  <- annot$Symbol[match(rownames(expr_all), annot$GeneID)]
z_all    <- t(scale(t(expr_all)))
meta142$comp <- colMeans(z_all[sym_all %in% neut,  , drop = FALSE]) -
                colMeans(z_all[sym_all %in% lymph, , drop = FALSE])

# ------------------------------------- 6. EV_RV vs Healthy (DU09 only) -------
du09 <- meta142$batch == "DU09"
m1 <- meta142[du09 & meta142$grp %in% c("Healthy", "EV_RV"), ]
print(table(droplevels(m1)$grp, droplevels(m1)$sex))

pc_ev_h <- do_de(m1, c("Healthy", "EV_RV"), ~ sex + grp)
cat("EV_RV vs Healthy: protein-coding genes with padj<0.05 & |log2FC|>1:",
    sum(pc_ev_h$padj < 0.05 & abs(pc_ev_h$log2FoldChange) > 1), "\n")

print(pc_ev_h %>% filter(Symbol %in% markers) %>%
        dplyr::select(Symbol, log2FoldChange, padj) %>% arrange(desc(log2FoldChange)))

# volcano
pc_ev_h$sig <- ifelse(pc_ev_h$padj < 0.05 & pc_ev_h$log2FoldChange >  1, "Up",
               ifelse(pc_ev_h$padj < 0.05 & pc_ev_h$log2FoldChange < -1, "Down", "NS"))
p_vol <- ggplot(pc_ev_h, aes(log2FoldChange, -log10(padj), colour = sig)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_colour_manual(values = c(Down = "steelblue", NS = "grey70", Up = "firebrick")) +
  labs(title = "Entero/Rhinovirus vs Healthy (DU09)",
       x = "log2 fold change", y = "-log10 adjusted p") +
  theme_bw()
ggsave(file.path(fig_dir, "volcano_EVRV_vs_Healthy.png"), p_vol, width = 7, height = 5, dpi = 300)

# enrichment (unadjusted)
f1 <- run_fgsea(pc_ev_h, sets_all)

# composition-adjusted version
pc_ev_h_adj <- do_de(meta142[du09 & meta142$grp %in% c("Healthy", "EV_RV"), ],
                     c("Healthy", "EV_RV"), ~ sex + comp + grp)
fa2 <- run_fgsea(pc_ev_h_adj, sets_all)

# what drives heme metabolism, and the "amino acid metabolism" hit?
le_heme <- f1$leadingEdge[[which(f1$pathway == "HALLMARK_HEME_METABOLISM")]]
print(head(pc_ev_h %>% filter(Symbol %in% le_heme) %>% arrange(padj) %>%
             dplyr::select(Symbol, log2FoldChange, padj), 25))

le_se <- f1$leadingEdge[[which(f1$pathway == "REACTOME_SELENOAMINO_ACID_METABOLISM")]]
cat("Fraction of selenoamino-acid leading edge that are RPL/RPS genes:",
    mean(grepl("^RP[LS]", le_se)), "\n")

# ----------------------------- 7. Flu vs Healthy (DU09) --------------------
pc_flu_h <- do_de(meta142[du09 & meta142$grp %in% c("Healthy", "Flu"), ],
                  c("Healthy", "Flu"), ~ sex + grp)
f3 <- run_fgsea(pc_flu_h, sets_all)

# --------------------- 8. EV_RV vs Flu: all batches, and DU09 only ---------
m2 <- droplevels(meta142[meta142$grp %in% c("EV_RV", "Flu"), ])
m2$batch <- factor(m2$batch)
print(table(m2$grp, m2$batch))
pc_ev_flu_all <- do_de(m2, c("Flu", "EV_RV"), ~ batch + sex + grp)
cat("EV_RV vs Flu (all batches): genes with padj<0.05 & |log2FC|>1:",
    sum(pc_ev_flu_all$padj < 0.05 & abs(pc_ev_flu_all$log2FoldChange) > 1), "\n")
f2 <- run_fgsea(pc_ev_flu_all, sets_all)

pc_ev_flu9 <- do_de(meta142[du09 & meta142$grp %in% c("Flu", "EV_RV"), ],
                    c("Flu", "EV_RV"), ~ sex + grp)
f4 <- run_fgsea(pc_ev_flu9, sets_all)

pc_ev_flu9_adj <- do_de(meta142[du09 & meta142$grp %in% c("Flu", "EV_RV"), ],
                        c("Flu", "EV_RV"), ~ sex + comp + grp)
fa1 <- run_fgsea(pc_ev_flu9_adj, sets_all)

# ----------------------------------- 9. summary table and heatmap ----------
key <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
         "HALLMARK_GLYCOLYSIS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
         "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_MTORC1_SIGNALING",
         "HALLMARK_HEME_METABOLISM", "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
         "REACTOME_PHOSPHOLIPID_METABOLISM", "REACTOME_SPHINGOLIPID_METABOLISM")

pick <- function(f, nm) {
  d <- as.data.frame(f)[, c("pathway", "NES", "padj")]
  d <- d[d$pathway %in% key, ]
  d$NES <- round(d$NES, 2); d$padj <- signif(d$padj, 2)
  names(d)[2:3] <- paste0(c("NES_", "padj_"), nm)
  d
}
tab <- Reduce(function(a, b) merge(a, b, by = "pathway", all = TRUE),
              list(pick(f1, "EV_vs_H"), pick(f3, "Flu_vs_H"),
                   pick(f2, "EV_vs_Flu_all"), pick(f4, "EV_vs_Flu_DU09")))
print(tab, row.names = FALSE)
write.csv(tab, file.path(res_dir, "summary_key_pathways.csv"), row.names = FALSE)

# composition-adjusted results for the same key pathways
adj_tab <- merge(pick(fa2, "EV_vs_H_adj"), pick(fa1, "EV_vs_Flu_DU09_adj"),
                 by = "pathway", all = TRUE)
print(adj_tab, row.names = FALSE)
write.csv(adj_tab, file.path(res_dir, "summary_key_pathways_composition_adjusted.csv"),
          row.names = FALSE)

cmp_list <- list(EV_vs_Healthy = f1, Flu_vs_Healthy = f3,
                 EV_vs_Flu_all = f2, EV_vs_Flu_DU09 = f4)
long <- bind_rows(lapply(names(cmp_list), function(n) {
  d <- as.data.frame(cmp_list[[n]])[, c("pathway", "NES", "padj")]
  d$comparison <- n; d
})) %>%
  filter(pathway %in% key) %>%
  mutate(label = ifelse(padj < 0.05, "*", ""),
         pathway = gsub("HALLMARK_|REACTOME_", "", pathway))

p_sum <- ggplot(long, aes(comparison, pathway, fill = NES)) +
  geom_tile(colour = "white") + geom_text(aes(label = label), size = 6) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = NULL, y = NULL,
       title = "Key pathways across comparisons (* padj < 0.05; not composition-adjusted)")
ggsave(file.path(fig_dir, "summary_pathways_heatmap.png"), p_sum, width = 8, height = 5, dpi = 300)

# ------------------------------------------------------------- 10. save ------
write.csv(pc_ev_h,       file.path(res_dir, "DE_EVRV_vs_Healthy_DU09.csv"),   row.names = FALSE)
write.csv(pc_ev_flu_all, file.path(res_dir, "DE_EVRV_vs_Flu_allbatches.csv"), row.names = FALSE)
write.csv(pc_ev_flu9,    file.path(res_dir, "DE_EVRV_vs_Flu_DU09.csv"),       row.names = FALSE)

sink(file.path(res_dir, "sessionInfo.txt")); print(sessionInfo()); sink()
cat("Done. Results in", res_dir, "\n")
