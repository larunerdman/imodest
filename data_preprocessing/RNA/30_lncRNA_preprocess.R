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
    make_option(c("--CANCER_TYPE"), type = "character", default = "KIRC",
                help = "Cancer set abbreviation name", metavar = "character"),
    make_option(c("--RNA_EXPRESSION_MATRIX"), type = "character",
                default = "/Volumes/USB128/RNA_preprocess_output/KIRC_RNASeq_gene.rds",
                help = "RNA expression data by isoform", metavar = "character"),
    make_option(c("--LNCRNA_EXPRESSION_PROFILE"), type = "character",
                default = "/Users/MJ/Desktop/lncRNA_preprocess/TCGA_lncRNA_expression/TCGA-KIRC-rnaexpr/TCGA-KIRC-rnaexpr.tsv",
                help = "RNA expression data by isoform", metavar = "character"),
    make_option(c("--LNCRNA_EXPRESSION_MATRIX"), type = "character",
                default = "/Users/MJ/Desktop/res.rds",
                help = "RNA expression data by isoform", metavar = "character")
);

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

Gene_Expression <- read_rds(opt$RNA_EXPRESSION_MATRIX)


LncRNA_Expression <- read_tsv(opt$LNCRNA_EXPRESSION_PROFILE, col_types = cols(
    .default = col_double(),
    Gene_ID = col_character()
))
LncRNA_Expression_2d <- LncRNA_Expression %>% select(-Gene_ID) %>% as.matrix.data.frame() %>% t()
colnames(LncRNA_Expression_2d) <- LncRNA_Expression$Gene_ID %>% str_extract("^ENSGR?\\d+")
rownames(LncRNA_Expression_2d) <- rownames(LncRNA_Expression_2d) %>% str_replace_all("-", ".")

LncRNA_Expression_2d <- LncRNA_Expression_2d[order(rownames(LncRNA_Expression_2d)),]



Patient_Sample_Mapping <- tibble(lncRNA_sample_id = rownames(LncRNA_Expression_2d)) %>%
    mutate(patient_id = str_extract(lncRNA_sample_id, "TCGA........$"), tumor = str_detect(lncRNA_sample_id, ".Tumor.")) %>%
    inner_join(tibble(sample_id = colnames(Gene_Expression)) %>%
                   mutate(patient_id = str_extract(sample_id, "^TCGA........"), tumor = str_detect(sample_id, "0..$")),
               by = c("patient_id", "tumor")) %>%
    arrange(lncRNA_sample_id)

assert_that(!any(duplicated(Patient_Sample_Mapping$sample_id)))

# find the index of x for which y contian, repeation is allowed
whichRepeat <- function(x, y) {
    result <- c()

    i <- integer(0)
    for (i in seq_along(x)) {
        occurence <- length(which(y == x[i]))
        if (occurence > 0) {
            result <- c(result, rep(i, times = occurence))
        }
    }

    result <- sort(result)

    return(result)
}

LncRNA_Expression_2d <- LncRNA_Expression_2d[whichRepeat(rownames(LncRNA_Expression_2d), Patient_Sample_Mapping$lncRNA_sample_id),]

rownames(LncRNA_Expression_2d) <- Patient_Sample_Mapping$sample_id

assert_that(!any(duplicated(colnames(LncRNA_Expression_2d))))
assert_that(!any(is.na(colnames(LncRNA_Expression_2d))))
assert_that(!any(is.na(rownames(LncRNA_Expression_2d))))
assert_that(!any(is.na(LncRNA_Expression_2d)))

LncRNA_Expression_2d <- LncRNA_Expression_2d[order(rownames(LncRNA_Expression_2d)),]
LncRNA_Expression_2d <- LncRNA_Expression_2d[,order(colnames(LncRNA_Expression_2d))]

write_rds(LncRNA_Expression_2d, opt$LNCRNA_EXPRESSION_MATRIX, compress = "gz")
