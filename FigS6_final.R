#Statistical enrichment test for L454W vs uninfected in KEGG Pathways
#naliebe@uw.edu

library(tidyverse)
library(ReactomePA)
library(org.Hs.eg.db)

setwd("~/Desktop/")

DE <- read_csv("~/Desktop/measles/sig_L454WvsNone.csv")
DE$abs <- abs(DE$log2FoldChange)
DE <-
  dplyr::filter(DE, abs > 1 &
                  padj < 0.01) #selects only genes with greater than 2-fold change and adjusted P of 0.01

hs <- org.Hs.eg.db
entrez <- AnnotationDbi::select(
  hs,
  keys = DE$rowname,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "SYMBOL"
)
entrez <-
  distinct(entrez, SYMBOL, .keep_all = TRUE) #removes duplicate entries if there was multiple mapping
DE$entrez <- entrez$ENTREZID

genes <- DE$entrez
kegg <- enrichPathway(gene = genes,
                      pvalueCutoff = 0.05,
                      readable = T)

png(
  filename = "~/Desktop/4S_kegg_top20.png",
  width = 12,
  height = 8,
  units = "in",
  res = 300
)
dotplot(kegg, showCategory = 20)
dev.off()
