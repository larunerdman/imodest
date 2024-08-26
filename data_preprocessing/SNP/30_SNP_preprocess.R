

########
# Collect SNP target matrix for genes/exons from all chromosomes
#

library(methods)
library(readr)
library(dplyr)
library(optparse)
library(stringr)
library(Matrix)
library(assertthat)


option_list = list(
    make_option(c("--INPUT_DIR"), type = "character",
                default = "/Users/MJ/Desktop/SNP_preprocess/temp/",
                help = "", metavar = "character"),
    make_option(c("--EXON_SPECIFIC"), type = "integer",
                default = 0,
                help = "", metavar = "character"),
    make_option(c("--SNP_TARGET_MATRIX"), type = "character",
                default = "/Users/MJ/Desktop/SNP_preprocess/50_SNP_target_matrix_gene.rds",
                help = "", metavar = "character")
)


print(opt <- parse_args(OptionParser(option_list = option_list)))

if (opt$EXON_SPECIFIC == 0) {
    SNP_target_gene_files <- list.files(opt$INPUT_DIR, full.names = T) %>% str_subset("target_matrix_gene")

    assert_that(length(SNP_target_gene_files) == 22)
    for (i in seq_along(SNP_target_gene_files)) {
        print(SNP_target_gene_files[i])
        SNP_target_gene_chr_i <- read_rds(SNP_target_gene_files[i])
        if (i == 1) {
            SNP_target_gene <- SNP_target_gene_chr_i
            num_targets <- nrow(SNP_target_gene_chr_i)
            num_probes <- ncol(SNP_target_gene_chr_i)
        } else {
            assert_that(num_targets == nrow(SNP_target_gene_chr_i))
            assert_that(num_probes == ncol(SNP_target_gene_chr_i))
            assert_that(all(rownames(SNP_target_gene) == rownames(SNP_target_gene_chr_i)))
            assert_that(all(colnames(SNP_target_gene) == colnames(SNP_target_gene_chr_i)))
            SNP_target_gene <- SNP_target_gene | SNP_target_gene_chr_i
        }
    }


    SNP_target_gene <- SNP_target_gene[order(rownames(SNP_target_gene)),]
    SNP_target_gene <- SNP_target_gene[,order(colnames(SNP_target_gene))]

    assert_that(!is.unsorted(rownames(SNP_target_gene)))
    assert_that(!is.unsorted(colnames(SNP_target_gene)))

    write_rds(SNP_target_gene, opt$SNP_TARGET_MATRIX, compress = "gz")

} else {
    SNP_target_exon_files <- list.files(opt$INPUT_DIR, full.names = T) %>% str_subset("target_matrix_exon")
    assert_that(length(SNP_target_gene_files) == 22)
    for (i in seq_along(SNP_target_exon_files)) {
        print(SNP_target_exon_files[i])
        SNP_target_exon_chr_i <- read_rds(SNP_target_exon_files[i])
        if (i == 1) {
            SNP_target_exon <- SNP_target_exon_chr_i
            num_targets <- nrow(SNP_target_exon_chr_i)
            num_probes <- ncol(SNP_target_exon_chr_i)
        } else {
            assert_that(num_targets == nrow(SNP_target_exon_chr_i))
            assert_that(num_probes == ncol(SNP_target_exon_chr_i))
            assert_that(all(rownames(SNP_target_exon) == rownames(SNP_target_exon_chr_i)))
            assert_that(all(colnames(SNP_target_exon) == colnames(SNP_target_exon_chr_i)))
            SNP_target_exon <- SNP_target_exon | SNP_target_exon_chr_i
        }
    }


    SNP_target_exon <- SNP_target_exon[order(rownames(SNP_target_exon)),]
    SNP_target_exon <- SNP_target_exon[,order(colnames(SNP_target_exon))]

    assert_that(!is.unsorted(rownames(SNP_target_exon)))
    assert_that(!is.unsorted(colnames(SNP_target_exon)))

    write_rds(SNP_target_exon, opt$SNP_TARGET_MATRIX, compress = "gz")
}


