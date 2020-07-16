## This code takes the BrainSpan gene FPKM dataset and our attempts to map the measles RNA-seq to GRCh37/hg19 and create an FPKM dataset
## and then generates a correlation matrix for a subset of genes
## It is modeled after Figure 1D in PMID 28009303
## Alex Greninger Vikas Peddu 10/20/19

## Front Matter
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(ggplot2)
library(reshape2)
library(dplyr)
library(zFPKM)
library(DESeq2)
library(pheatmap)
library(tidyverse)

##Preprocessing

setwd("/Users/gerbix/Downloads/brainspan")

## Read in the measles and BrainSpan gene FPKM datasets for GenCode v10 as well as the row and column metadata for BrainSpan
## These files are available https://www.brainspan.org/api/v2/well_known_file_download/267666525


### VIKAS: all_RPKMs.csv comes is made from the featurecounts output (grch37) of the samples
### This is created from a folder of featurecounts output txt files using featurecounts_to_rpkm.r

measles <-
  read.table(
    "/Users/gerbix/Downloads/brainspan/output/all_RPKMs.csv",
    sep = ',',
    header = TRUE,
    stringsAsFactors = FALSE
  )
measles$X <- NULL
measles[grepl("MF", colnames(measles))] <- NULL

#removing extra samples
measles <-
  measles[-which(names(measles) %in% c('AA001','AA002','AA004','AA003', "AA009", "AA010", "AA011", "AA015"))]


measles_header <- colnames(measles)[2:ncol(measles)]
#measles_header <- c("MA001","MA002","MA003","MA004","MA005","MA006","MA007","MA008","MA009","MA010")
expression_matrix <-
  read.table(
    "expression_matrix.csv",
    sep = ',',
    row.names = NULL,
    stringsAsFactors = FALSE
  )
column_metadata <-
  read.table(
    "columns_metadata.csv",
    sep = ',',
    header = TRUE,
    stringsAsFactors = FALSE
  )

## Create a variable our_column within the column_metadata that links the column names in expression matrix to the column metadata
column_metadata$our_column <- colnames(expression_matrix)[2:525]
row_metadata <-
  read.table(
    "rows_metadata.csv",
    sep = ',',
    header = TRUE,
    stringsAsFactors = FALSE
  )

## Add the row metadata as front matter.  We will get rid of it soon though.
expression_column <- cbind(row_metadata, expression_matrix)

## Remove row numbers and fix the gene symbol column name to allow matching based on it
#expression_column$V1 <- NULL
names(measles)[names(measles) == 'Gene_symbol'] <- 'ensembl_gene_id'


measles$ensembl_gene_id <-
  sapply(strsplit(as.character(measles$ensembl_gene_id), "\\."), `[`, 1)

## We will do an inner join of the BrainSpan and measles FPKM data to only include common genes
inner_data <-
  inner_join(expression_column, measles, by = c('ensembl_gene_id'))


## create a new dataframe just_inner_data that does not have the gene data in the front for correlation matrix
just_inner_data <- inner_data[, c(6:544)]

##Filter reads based on sum total of fpkm counts across all 534 experiments

fpkm_sum_filter_threshold <- 1000
expressed <-
  just_inner_data[rowSums(just_inner_data[, -1]) > fpkm_sum_filter_threshold, ]

## Transform inner_data by log2 after adding a pseudocount to every value to prevent -Inf
inner_data_log2 <- log2(expressed + 1)

#Check to see how the filtering affects the count
count(expressed$V2 == 0)
count(expressed$V2 > 0)


##check a histogram of how the values look
ggplot(inner_data_log2, aes(x = inner_data_log2$MA001)) + geom_histogram(binwidth =
                                                                           .01)


## Make correlation matrices based on the gene FPKM data
cormat_inner_log2 <- round(cor(inner_data_log2), 3)


## Remove the measles FPKM columns
cormat_inner_log2 <-
  cormat_inner_log2[, -c(526:ncol(cormat_inner_log2))]

#Make a dataframe version of the correlation matrix for ease of sorting
cormat_inner_log2_df <- as.data.frame(cormat_inner_log2)

#Melt the correlation matrix for tidy ggplot data
melted_cormat_inner_log2 <- melt(cormat_inner_log2)

#Display
ggplot(data = melted_cormat_inner_log2[melted_cormat_inner_log2$Var1 %in% measles_header, ], aes(x =
                                                                                                   Var1, y = Var2, fill = value)) +
  geom_tile()

pheatmap(
  cormat_inner_log2[c(526:539), ],
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  width = 8,
  height = 4,
  fontsize = 4
)

#write.table(cormat_inner_log2, "correlation_table.csv", sep = ',')

cormat_inner_log2_df_MA001 <-
  cormat_inner_log2_df[, order(cormat_inner_log2_df[530, ], decreasing = T)]


top100_list <- colnames(cormat_inner_log2_df_MA001[1:100])

subset_top100 <-
  column_metadata[match(top100_list, column_metadata$our_column), ][c(1:100), ]

labels_top100 <-
  paste(
    as.character(subset_top100$structure_acronym),
    as.character(subset_top100$gender),
    as.character(subset_top100$age)
  )


### VIKAS: This is the heatmap figure
library(tidyverse)

for_heatmap <- cormat_inner_log2_df_MA001[c(526:539), c(1:100)]
name_order <-
  c(
    "MA001",
    "MA002",
    "FB018",
    #"AA004",
    "MA007",
    "FB009",
    #"AA001",
    "MA005",
    #"AA002",
    "MA009",
    "MA003",
    "MA004",
    "FB017",
    "MA008",
    "FB001",
    "MA006",
    "MA010"
  )
# A = AA, B = MA, C = FB
differentiation<- data.frame(
  Differentiation = c(
    "B",
    "B",
    "C",
    #"A",
    "B",
    "C",
    #"A",
    "B",
    #"A",
    "B",
    "B",
    "B",
    "C",
    "B",
    "C",
    "B",
    "B"
    )
  )
rownames(differentiation) = name_order



temp <- tibble::rownames_to_column(for_heatmap)

reordered_for_heatmap <- temp %>%
  dplyr::slice(match(name_order, rowname))
rownames(reordered_for_heatmap) <- reordered_for_heatmap$rowname
reordered_for_heatmap$rowname <- NULL


# png(
#   filename = "brainspan_labeled_heatmap_final.png",
#   width = 8,
#   height = 6,
#   units = "in",
#   res = 300
# )
pheatmap(
  reordered_for_heatmap,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  labels_col = labels_top100,
  #annotation_row = differentiation,
  labels_row = c(
    "FA10 uninfected (1)",
    "FA10 uninfected (2)",
    "FA10 uninfected (3)",
    "FA10 + MeV wt (1)",
    "FA10 + MeV wt (2)",
    "FA10 + MeV L454W",
    "FA10 + MeV L454W/E455G",
    "FA11 uninfected (1)",
    "FA11 uninfected (2)",
    "FA11 uninfected (3)",
    "FA11 + MeV wt (1)",
    "FA11 + MeV wt (2)",
    "FA11 + MeV L454W",
    "FA11 + MeV L454W/E455G"
  ),
  width = 11,
  filename = "brainspan_labeled_heatmap_final.png",
  height = 8,
  fontsize = 7,
  fontsize_col = 5
)

dev.off()
#save_pheatmap_pdf(brainspan_labeled_heatmap, "brainspan_labeled_heatmap.pdf")
