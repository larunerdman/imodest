
#####
#   Generate 2D matrix, each value is a lncRNA isoform expression for an individual
#       row x col = individual x lncRNA isoform

library(tidyr)    # separate
library(optparse) # make_option
library(dplyr, warn.conflicts = F)    # select
library(stringr)
library(readr)    # read_tsv
library(tibble)
library(assertthat)




####
writeLines("10 - Input specification")
################################################################

option_list = list(
    make_option(c("--CANCER_TYPE"), type = "character", default = "KIRC",
                help = "Cancer set abbreviation name", metavar = "character"),
    make_option(c("--RNA_EXPRESSION_MATRIX"), type = "character",
                default = "/Users/MJ/Desktop/KIRC_RNASeq_gene.rds",
                help = "RNA expression data by isoform", metavar = "character"),
    make_option(c("--LNCRNA_ID"), type = "character",
                default = "/Users/MJ/Desktop/lncRNA_preprocess/10_lncRNA_ids.rds",
                help = "RNA expression data by isoform", metavar = "character"),
    make_option(c("--LNCRNA_EXPRESSION_MATRIX"), type = "character",
                default = "/Users/MJ/Desktop/lncRNA_preprocess/res.rds",
                help = "RNA expression data by isoform", metavar = "character")
);

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)


####
writeLines("20 - Output specification")
################################################################
OUTPUT <- str_c("/Volumes/USB128/lncRNA_preprocess_output/", opt$CANCER_TYPE, "_lncRNA_expression.rds")


####
writeLines("40 - Collect lncRNA isoform expression")
################################################################
RNA_expression <- read_rds(opt$RNA_EXPRESSION_MATRIX)
LncRNA_ids <- sort(read_rds(opt$LNCRNA_ID))

# Select availible lncRNA isoform in both target data and expression data
LncRNA_Expression <- RNA_expression[rownames(RNA_expression) %in% LncRNA_ids,]

assert_that(!is.unsorted(rownames(RNA_expression)))
assert_that(!is.unsorted(colnames(RNA_expression)))
assert_that(!is.unsorted(rownames(LncRNA_Expression)))
assert_that(!is.unsorted(colnames(LncRNA_Expression)))


#### 5) Write to local
################################################################
writeLines("50 - output lncRNA expression data to local")
write_rds(t(LncRNA_Expression), opt$LNCRNA_EXPRESSION_MATRIX, compress = "gz")
