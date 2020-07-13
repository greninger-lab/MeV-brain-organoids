#Heatmap of samples with different mutations
#naliebe@uw.edu

library(tidyverse)
library(DESeq2)
library(pheatmap)

setwd("/Users/vikas/Downloads/Fig_S4/with_new_mutant")

metadata <-
  read_csv("/Users/vikas/Downloads/Fig_S4/with_new_mutant/S4_metadata.csv")

FB <- read_tsv("/Users/vikas/Downloads/Fig_S4/FB_counts.txt")
MA <- read_tsv("/Users/vikas/Downloads/Fig_S4/MA_counts.txt")

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
  rowSums(counts(dds)) >= 12  #pre-filter to remove any genes without avg 1 count per sample
dds <- dds[keep, ]
dds <-
  estimateSizeFactors(dds) #since not doing DE here, just need normalized counts
norm_counts <- as.data.frame(counts(dds, normalized = TRUE))

#read in sig genes from F454WvsNone to get list of top 50 from the larger experiment in main text, then will subset current norm counts for heatmap
L454W <-
  read_csv("/Users/vikas/Downloads/Fig_S4/sig_L454WvsNone.csv")
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



metadata <- rownames_to_column(metadata)

metadata$log_rpm <- log10(metadata$RPM + 1)
new_order <-
  c(
    "FB011",
    "FB009",
    "MA007",
    "MA009",
    "FB001",
    "FB003",
    "MA008",
    "MA010",
    "MA006",
    "MA005",
    "FB012",
    "FB004"
  )


reordered_metadata <- metadata %>%
  dplyr::slice(match(new_order, rowname))
reordered_metadata$order <- c(1:length(reordered_metadata$rowname))
reordered_metadata$rowname <-
  factor(reordered_metadata$rowname, levels = reordered_metadata$rowname)

barplot = ggplot(reordered_metadata,
                 aes (
                   x = rowname,
                   y = log_rpm,
                   order = as.factor(order)
                 )) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("log10(RPM + 1)") +
  geom_text(aes(label = round(RPM, digits = 1)), nudge_y = .1) +
  scale_y_continuous(expand = c(0, 0))
barplot

ggsave(plot = barplot,
       'S4_barplot.pdf',
       height = 2,
       width = 4)
