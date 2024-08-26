
########
# Genereate lncRNA expression, extend from 10_lncRNA_preprocess


library(methods)
library(readr)
library(optparse) # make_option
library(dplyr)
library(stringr)
library(assertthat)
####
writeLines("10 - Input specification")
################################################################

option_list = list(
    make_option(c("--CANCER_TYPE"), type = "character", default = "COAD",
                help = "Cancer set abbreviation name", metavar = "character"),
    make_option(c("--RNA_EXPRESSION_MATRIX"), type = "character",
                default = "/Volumes/USB128/RNA_preprocess_output/COAD_RNASeq_gene.rds",
                help = "RNA expression data by isoform", metavar = "character"),
    make_option(c("--LNCRNA_IDS"), type = "character",
                default = "/Users/MJ/Desktop/lncRNA_preprocess/10_lncRNA_ids.rds",
                help = "RNA expression data by isoform", metavar = "character"),
    make_option(c("--LNCRNA_EXPRESSION_MATRIX"), type = "character",
                default = "/Users/MJ/Desktop/res.rds",
                help = "RNA expression data by isoform", metavar = "character")
);

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

Gene_Expression <- read_rds(opt$RNA_EXPRESSION_MATRIX)


LncRNA_Expression_2d <- Gene_Expression[rownames(Gene_Expression) %in% read_rds(opt$LNCRNA_IDS),] %>% t()


assert_that(!any(duplicated(colnames(LncRNA_Expression_2d))))
assert_that(!any(is.na(colnames(LncRNA_Expression_2d))))
assert_that(!any(is.na(rownames(LncRNA_Expression_2d))))
assert_that(!any(is.na(LncRNA_Expression_2d)))

LncRNA_Expression_2d <- LncRNA_Expression_2d[order(rownames(LncRNA_Expression_2d)),]
LncRNA_Expression_2d <- LncRNA_Expression_2d[,order(colnames(LncRNA_Expression_2d))]

write_rds(LncRNA_Expression_2d, opt$LNCRNA_EXPRESSION_MATRIX, compress = "gz")
