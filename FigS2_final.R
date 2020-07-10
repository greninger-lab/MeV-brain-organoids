#go from transcript output from Kallisto to DE and heatmap
#naliebe@uw.edu

library(tidyverse)
library(biomaRt)
library(EnsDb.Mmusculus.v79)
library(tximport)
library(ensembldb)
library(forestmangr)
library(DESeq2)
library(tidyverse)
library(pheatmap)

#Part One: first collapse ensembl transcripts to gene, then convert ensembl gene ID to MGI gene name

setwd("~/Desktop/")

f_path <- "~/Desktop/CM_transcript_tsvs"
filenames <-
  list.files(f_path, pattern = '*_abundance.tsv', full.names = TRUE) #defines filenames as the same as from tsv file

# Annotates by transcript
txdb <- EnsDb.Mmusculus.v79 #makes connection to Ensembl database
tx2gene <- transcripts(txdb, return.type = "DataFrame")
tx2gene <-
  as.data.frame(tx2gene[, c("tx_id", "gene_id")]) #these two lines make a table of transcript ids and associated gene IDs
txi <-
  tximport(
    as.character(filenames),
    type = "kallisto",
    tx2gene = tx2gene,
    ignoreTxVersion = TRUE
  )
df <-
  txi$counts #creates matrix of counts for all samples; associates them with just the counts from each tsv
filenames_short <- list.files(f_path, pattern = "*_abundance.tsv")
filenames_short <- substr(filenames_short, 1, 5)
colnames(df) <-
  filenames_short #adds column names in same order as files
df <- as.data.frame(df)

mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("mmusculus_gene_ensembl", mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filter = "ensembl_gene_id",
  values = rownames(df),
  uniqueRows = TRUE
)

df$mgi_symbol <- rownames(df)
merged <-
  merge(df, annotLookup, by.x = "mgi_symbol", by.y = "ensembl_gene_id")
merged$mgi_symbol <- NULL

# Merges counts by gene symbol
counts <-
  aggregate(merged[, c(1:ncol(merged) - 1)],
            by = list(Category = merged$mgi_symbol),
            FUN = sum)
colnames(counts)[1] <- 'gene'

#tidying
counts <-
  round_df(counts, digits = 0, rf = "round") #rounds counts data off
counts <-
  counts[-1, ] #for some reason there was no gene name associated with the first row.......

# Writes count table
write.csv(
  counts,
  "~/Desktop/mouse_measles_raw_counts.csv",
  quote = FALSE,
  row.names = FALSE
)

##################
#Part Two: Differential Expression Analysis
metadata <- read_csv("~/Desktop/mouse_metadata.csv")
metadata <- column_to_rownames(metadata, "Sample")

#counts from annotation, above:
counts <-
  rownames_to_column(counts) #tidying, R was throwing error with next line
counts <- column_to_rownames(counts, "gene")

#if just reading in raw counts from csv:
#counts <- read_csv("~/Desktop/mouse_measles_raw_counts.csv")
#counts <- column_to_rownames(counts, "gene")

#normalization and results
dds <-
  DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = ~ Site + Infection
  ) #includes site expt was performed as confounder
dds$Infection <-
  factor(dds$Infection, levels = c("none", "WT", "L454W")) #first position is reference; if not explicit, will be alphabetical
keep <-
  rowSums(counts(dds)) >= 10  #pre-filter to remove any genes without avg 1 count per sample
dds <- dds[keep, ]
dds <- DESeq(dds, parallel = FALSE, quiet = FALSE)
norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
write.csv(norm_counts, "~/Desktop/mouse_measles_norm_counts.csv")

#contrasts

res_WT_vs_none <-
  as.data.frame(results(dds, contrast = c("Infection", "WT", "none")))
write.csv(results(dds, contrast = c("Infection", "WT", "none")),
          "~/Desktop/mouse_results_WTvsNone.csv")
res_WT_vs_none_Sig <-
  res_WT_vs_none[which(res_WT_vs_none$padj < 0.1),]
res_WT_vs_none_Sig <- as.data.frame(res_WT_vs_none_Sig)
res_WT_vs_none_Sig <- rownames_to_column(res_WT_vs_none_Sig)
res_WT_vs_none_Sig <- dplyr::arrange(res_WT_vs_none_Sig, padj)
write.csv(res_WT_vs_none_Sig, "~/Desktop/mouse_sig_WTvsNone.csv")

res_L454W_vs_none <-
  as.data.frame(results(dds, contrast = c("Infection", "L454W", "none")))
write.csv(results(dds, contrast = c("Infection", "L454W", "none")),
          "~/Desktop/mouse_results_L454WvsNone.csv")
res_L454W_vs_none_Sig <-
  res_L454W_vs_none[which(res_L454W_vs_none$padj < 0.1),]
res_L454W_vs_none_Sig <- as.data.frame(res_L454W_vs_none_Sig)
res_L454W_vs_none_Sig <- rownames_to_column(res_L454W_vs_none_Sig)
res_L454W_vs_none_Sig <- dplyr::arrange(res_L454W_vs_none_Sig, padj)
write.csv(res_L454W_vs_none_Sig,
          "~/Desktop/mouse_sig_L454WvsNone.csv")

res_L454W_vs_WT <-
  as.data.frame(results(dds, contrast = c("Infection", "L454W", "WT")))
write.csv(results(dds, contrast = c("Infection", "L454W", "WT")),
          "~/Desktop/mouse_results_L454WvsWT.csv")
res_L454W_vs_WT_Sig <-
  res_L454W_vs_WT[which(res_L454W_vs_WT$padj < 0.1),]
res_L454W_vs_WT_Sig <- as.data.frame(res_L454W_vs_WT_Sig)
res_L454W_vs_WT_Sig <- rownames_to_column(res_L454W_vs_WT_Sig)
res_L454W_vs_WT_Sig <- dplyr::arrange(res_L454W_vs_WT_Sig, padj)
write.csv(res_L454W_vs_WT_Sig, "~/Desktop/mouse_sig_L454WvsWT.csv")


#heatmap of top 50 most significant genes in L454W vs uninfected
df2 <- metadata[, c("Infection", "Site")]
var <- res_L454W_vs_none_Sig$rowname[1:50]
de <- norm_counts[var, ] + 1 #avoids divide by 0 issue
fc <- de / rowMeans(de)
l2fc <- log2(fc)

#for heatmap without "Site" - run these lines:
df2 <- rownames_to_column(df2)
df2 <- df2[, -3]
df2 <- column_to_rownames(df2)

png(
  filename = "~/Desktop/mouse_top50DE_heatmap.png",
  width = 6,
  height = 8,
  units = "in",
  res = 300
)
print(pheatmap(l2fc, annotation_col = df2, cluster_cols = TRUE)) #creates dendrogram for samples and genes. if cluster_cols=FALSE it does not cluster
dev.off()
