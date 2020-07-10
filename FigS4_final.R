#Heatmap of samples with different mutations
#naliebe@uw.edu

library(tidyverse)
library(DESeq2)
library(pheatmap)

metadata <- read_csv("~/Desktop/S4_metadata.csv")

FB <- read_tsv("~/Desktop/measles/FB_counts.txt")
MA <- read_tsv("~/Desktop/measles/MA_counts.txt")

counts <- inner_join(FB, MA, by = "gene")

counts <- counts[, c("gene", metadata$Sample)]

counts <- column_to_rownames(counts, "gene")
metadata <- column_to_rownames(metadata, "Sample")

dds <-
  DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ Organoid + Infection
  ) #includes organoid of sample origin as confounder
keep <-
  rowSums(counts(dds)) >= 10  #pre-filter to remove any genes without avg 1 count per sample
dds <- dds[keep, ]
dds <-
  estimateSizeFactors(dds) #since not doing DE here, just need normalized counts
norm_counts <- as.data.frame(counts(dds, normalized = TRUE))

#read in sig genes from F454WvsNone to get list of top 50 from the larger experiment in main text, then will subset current norm counts for heatmap
L454W <- read_csv("~/Desktop/measles/sig_L454WvsNone.csv")
genes <- L454W$rowname[1:50]

#Heatmap
df <- metadata[, c("Organoid", "Infection")]
de <- norm_counts[genes, ] + 1 #avoids divide by 0 issue
fc <- de / rowMeans(de)
l2fc <- log2(fc)

png(
  filename = "~/Desktop/S4_measles_mutants_top50DE_heatmap.png",
  width = 6,
  height = 8,
  units = "in",
  res = 300
)
print(pheatmap(l2fc, annotation_col = df, cluster_cols = TRUE)) #creates dendrogram for samples and genes. if cluster_cols=FALSE it does not cluster
dev.off()
