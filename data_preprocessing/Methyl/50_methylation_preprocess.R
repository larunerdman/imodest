

########
# Collect methylation target matrix for genes/exons in one chromosome
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
                default = "/Users/MJ/Desktop/methylation_preprocess/temp/",
                help = "", metavar = "character"),
    make_option(c("--EXON_SPECIFIC"), type = "integer",
                default = 0,
                help = "", metavar = "character"),
    make_option(c("--METHYLATION_TARGET_MATRIX"), type = "character",
                default = "/Users/MJ/Desktop/methylation_preprocess/50_methylation_target_matrix_gene.rds",
                help = "", metavar = "character")
)


print(opt <- parse_args(OptionParser(option_list = option_list)))

if (opt$EXON_SPECIFIC == 0) {
    methylation_target_gene_files <- list.files(opt$INPUT_DIR, full.names = T) %>% str_subset("target_matrix_gene")

    for (i in seq_along(methylation_target_gene_files)) {
        print(methylation_target_gene_files[i])
        methylation_target_gene_chr_i <- read_rds(methylation_target_gene_files[i])
        if (i == 1) {
            methylation_target_gene <- methylation_target_gene_chr_i
            num_targets <- nrow(methylation_target_gene_chr_i)
            num_probes <- ncol(methylation_target_gene_chr_i)
        } else {
            assert_that(num_targets == nrow(methylation_target_gene_chr_i))
            assert_that(num_probes == ncol(methylation_target_gene_chr_i))
            assert_that(all(rownames(methylation_target_gene) == rownames(methylation_target_gene_chr_i)))
            assert_that(all(colnames(methylation_target_gene) == colnames(methylation_target_gene_chr_i)))
            methylation_target_gene <- methylation_target_gene | methylation_target_gene_chr_i
        }
    }


    methylation_target_gene <- methylation_target_gene[order(rownames(methylation_target_gene)),]
    methylation_target_gene <- methylation_target_gene[,order(colnames(methylation_target_gene))]

    assert_that(!is.unsorted(rownames(methylation_target_gene)))
    assert_that(!is.unsorted(colnames(methylation_target_gene)))

    write_rds(methylation_target_gene, opt$METHYLATION_TARGET_MATRIX, compress = "gz")

} else {
    methylation_target_exon_files <- list.files(opt$INPUT_DIR, full.names = T) %>% str_subset("target_matrix_exon")

    for (i in seq_along(methylation_target_exon_files)) {
        print(methylation_target_exon_files[i])
        methylation_target_exon_chr_i <- read_rds(methylation_target_exon_files[i])
        if (i == 1) {
            methylation_target_exon <- methylation_target_exon_chr_i
            num_targets <- nrow(methylation_target_exon_chr_i)
            num_probes <- ncol(methylation_target_exon_chr_i)
        } else {
            assert_that(num_targets == nrow(methylation_target_exon_chr_i))
            assert_that(num_probes == ncol(methylation_target_exon_chr_i))
            assert_that(all(rownames(methylation_target_exon) == rownames(methylation_target_exon_chr_i)))
            assert_that(all(colnames(methylation_target_exon) == colnames(methylation_target_exon_chr_i)))
            methylation_target_exon <- methylation_target_exon | methylation_target_exon_chr_i
        }
    }


    methylation_target_exon <- methylation_target_exon[order(rownames(methylation_target_exon)),]
    methylation_target_exon <- methylation_target_exon[,order(colnames(methylation_target_exon))]

    assert_that(!is.unsorted(rownames(methylation_target_exon)))
    assert_that(!is.unsorted(colnames(methylation_target_exon)))

    write_rds(methylation_target_exon, opt$METHYLATION_TARGET_MATRIX, compress = "gz")
}


