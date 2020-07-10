#DE of genes in brain organoids upon infection with WT or L454W measles infection
#naliebe@uw.edu

setwd("~/Desktop/measles/")

library(tidyverse)
library(DESeq2)
library(pheatmap)

metadata <- read_csv("~/Desktop/measles/metadata.csv")

AA <- read_tsv("~/Desktop/measles/AA_counts.txt")
FB <- read_tsv("~/Desktop/measles/FB_counts.txt")
MA <- read_tsv("~/Desktop/measles/MA_counts.txt")

counts <- inner_join(AA, FB, by = "gene")
counts <- inner_join(counts, MA, by = "gene")

counts <- counts[, c("gene", metadata$Sample)]

counts <- column_to_rownames(counts, "gene")
metadata <- column_to_rownames(metadata, "Sample")

#normalization and results
dds <-
  DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ Run + Organoid + Infection
  ) #last is for DEG, earlier are as counfounders (conf1 + conf2 + forDEG)
dds$Infection <-
  factor(dds$Infection, levels = c("none", "WT", "L454W")) #first position is reference; if not explicit, will be alphabetical
keep <-
  rowSums(counts(dds)) >= 15  #pre-filter to remove any genes without avg 1 count per sample
dds <- dds[keep, ]
dds <- DESeq(dds, parallel = FALSE, quiet = FALSE)
norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
write.csv(norm_counts, "~/Desktop/measles/measles_norm_counts.csv")

#contrasts
res_WT_vs_none <-
  as.data.frame(results(dds, contrast = c("Infection", "WT", "none")))
write.csv(results(dds, contrast = c("Infection", "WT", "none")),
          "~/Desktop/measles/results_WTvsNone.csv")
res_WT_vs_none_Sig <-
  res_WT_vs_none[which(res_WT_vs_none$padj < 0.1),]
res_WT_vs_none_Sig <- as.data.frame(res_WT_vs_none_Sig)
res_WT_vs_none_Sig <- rownames_to_column(res_WT_vs_none_Sig)
res_WT_vs_none_Sig <- dplyr::arrange(res_WT_vs_none_Sig, padj)
write.csv(res_WT_vs_none_Sig, "~/Desktop/measles/sig_WTvsNone.csv")

res_L454W_vs_none <-
  as.data.frame(results(dds, contrast = c("Infection", "L454W", "none")))
write.csv(results(dds, contrast = c("Infection", "L454W", "none")),
          "~/Desktop/measles/results_L454WvsNone.csv")
res_L454W_vs_none_Sig <-
  res_L454W_vs_none[which(res_L454W_vs_none$padj < 0.1),]
res_L454W_vs_none_Sig <- as.data.frame(res_L454W_vs_none_Sig)
res_L454W_vs_none_Sig <- rownames_to_column(res_L454W_vs_none_Sig)
res_L454W_vs_none_Sig <- dplyr::arrange(res_L454W_vs_none_Sig, padj)
write.csv(res_L454W_vs_none_Sig,
          "~/Desktop/measles/sig_L454WvsNone.csv")

res_L454W_vs_WT <-
  as.data.frame(results(dds, contrast = c("Infection", "L454W", "WT")))
write.csv(results(dds, contrast = c("Infection", "L454W", "WT")),
          "~/Desktop/measles/results_L454WvsWT.csv")
res_L454W_vs_WT_Sig <-
  res_L454W_vs_WT[which(res_L454W_vs_WT$padj < 0.1),]
res_L454W_vs_WT_Sig <- as.data.frame(res_L454W_vs_WT_Sig)
res_L454W_vs_WT_Sig <- rownames_to_column(res_L454W_vs_WT_Sig)
res_L454W_vs_WT_Sig <- dplyr::arrange(res_L454W_vs_WT_Sig, padj)
write.csv(res_L454W_vs_WT_Sig, "~/Desktop/measles/sig_L454WvsWT.csv")

#heatmap of top 50 most significant genes in L454W vs uninfected
df <- metadata[, c("Organoid", "Infection")]
var <- res_L454W_vs_none_Sig$rowname[1:50]
de <- norm_counts[var, ] + 1 #avoids divide by 0 issue
fc <- de / rowMeans(de)
l2fc <- log2(fc)

png(
  filename = "~/Desktop/measles/measles_top50DE_heatmap.png",
  width = 6,
  height = 8,
  units = "in",
  res = 300
)
print(pheatmap(l2fc, annotation_col = df, cluster_cols = TRUE)) #creates dendrogram for samples and genes. if cluster_cols=FALSE it does not cluster
dev.off()
