

########
# Collect lncRNA target matrix for genes/exons in each chromosome
#
library(methods)
library(readr)
library(dplyr)
library(optparse)
library(stringr)
library(Matrix)
library(assertthat)
library(compiler)

option_list = list(
    make_option(c("--INPUT_DIR"), type = "character",
                default = "/Users/MJ/Desktop/lncRNA/",
                help = "", metavar = "character"),
    make_option(c("--EXON_SPECIFIC"), type = "integer",
                default = 1,
                help = "", metavar = "character"),
    make_option(c("--OUTPUT_DIR"), type = "character",
                default = "/Users/MJ/Desktop/lncRNA_preprocess/",
                help = "", metavar = "character")
)


opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

full_join_matrix <- function(m1, m2) {

    m3 <- sparseMatrix(i = c(summary(m1)$i, summary(m2)$i + nrow(m1)),
                       j = c(summary(m1)$j, summary(m2)$j + ncol(m1)),
                       x = c(summary(m1)$x, summary(m2)$x),
                       dims = dim(m1) + dim(m2),
                       dimnames = list(c(rownames(m1), rownames(m2)),
                                       c(colnames(m1), colnames(m2))))
    assert_that(!any(duplicated(rownames(m3))))
    assert_that(!any(duplicated(colnames(m3))))
    return(m3)
}
full_join_matrix <- cmpfun(full_join_matrix)

if (opt$EXON_SPECIFIC == 0) {
    lncRNA_target_gene_files <- list.files(opt$INPUT_DIR, full.names = T) %>% str_subset("gene")

    for (i in seq_along(lncRNA_target_gene_files)) {
        print(i)
        if (i == 1) {
            lncRNA_Target_Gene <- read_rds(lncRNA_target_gene_files[i])
        } else {
            lncRNA_Target_Gene <- full_join_matrix(lncRNA_Target_Gene, read_rds(lncRNA_target_gene_files[i]))
        }
    }

    lncRNA_Target_Gene <- lncRNA_Target_Gene[order(rownames(lncRNA_Target_Gene)),]
    lncRNA_Target_Gene <- lncRNA_Target_Gene[,order(colnames(lncRNA_Target_Gene))]

    assert_that(!is.unsorted(rownames(lncRNA_Target_Gene)))
    assert_that(!is.unsorted(colnames(lncRNA_Target_Gene)))

    write_rds(lncRNA_Target_Gene, str_c(opt$OUTPUT_DIR, "50_lncRNA_target_score_gene.rds"), compress = "gz")

} else {
    lncRNA_target_exon_files <- list.files(opt$INPUT_DIR, full.names = T) %>% str_subset("exon")

    for (i in seq_along(lncRNA_target_exon_files)) {
        print(i)
        if (i == 1) {
            lncRNA_Target_Exon <- read_rds(lncRNA_target_exon_files[i])
        } else {
            lncRNA_Target_Exon <- full_join_matrix(lncRNA_Target_Exon, read_rds(lncRNA_target_exon_files[i]))
        }
    }

    lncRNA_Target_Exon <- lncRNA_Target_Exon[order(rownames(lncRNA_Target_Exon)),]
    lncRNA_Target_Exon <- lncRNA_Target_Exon[,order(colnames(lncRNA_Target_Exon))]

    assert_that(!is.unsorted(rownames(lncRNA_Target_Exon)))
    assert_that(!is.unsorted(colnames(lncRNA_Target_Exon)))

    write_rds(lncRNA_Target_Exon, str_c(opt$OUTPUT_DIR, "50_lncRNA_target_score_exon.rds"), compress = "gz")

}


