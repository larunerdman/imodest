
# Collect corrected methylation matrix for normal sample and tumor sample
## prerequisite: 35
library(methods)
library(compiler)
library(optparse)

library(readr)
library(stringr)

library(assertthat)

library(dplyr)

library(Matrix)

library(iterators)
library(foreach)
library(parallel)
library(doParallel)

library(peer)
library(ggplot2)

logging <- function(..., sep = " ", level = 1) {
    writeLines(str_c(str_c(str_c(rep("-", level), collapse = ""), ">", sep = ""),..., sep = sep))

    if (level <= 3) {
        writeLines(str_c(str_c(rep(" ", level+2), collapse = ""),
                         "[", Sys.time(), "]"))
    }


    if (level <= 2) {
        for (i in seq_len(5)) {
            gc()
        }
    }

}
logging <- cmpfun(logging)

# Read in options from command line
option_list = list(
    make_option(c("--CANCER"), type = "character",
                default = "HNSC",
                help="full path to the methylation data", metavar = "character"),
    make_option(c("--METHYLATION_MATRIX_CORRECTED_1"), type = "character",
                default = "/Volumes/USB128/MoR/methylation_preprocess_output/HNSC_methylation_matrix_corrected.rds",
                help="full path to the tumor clinic file", metavar = "character"),
    make_option(c("--METHYLATION_MATRIX_CORRECTED_2"), type = "character",
                default = "/Volumes/USB128/MoR/methylation_preprocess_output/HNSC_methylation_matrix_corrected.rds",
                help="full path to the tumor clinic file", metavar = "character"),
    make_option(c("--METHYLATION_MATRIX_CORRECTED"), type = "character",
                default = "/Volumes/USB128/MoR/methylation_preprocess_output/HNSC_methylation_matrix_corrected.rds",
                help="full path to the tumor clinic file", metavar = "character")
)

print(opt <- parse_args(OptionParser(option_list=option_list)))
methylation_matrix_normal <- NULL
methylation_matrix_tumor <- NULL


if (file.exists(opt$METHYLATION_MATRIX_CORRECTED_2)) {
    logging("Reading normal samples")
    methylation_matrix_normal <- read_rds(opt$METHYLATION_MATRIX_CORRECTED_2)
}



if (file.exists(opt$METHYLATION_MATRIX_CORRECTED_1)) {
    logging("Reading tumor samples")
    methylation_matrix_tumor <- read_rds(opt$METHYLATION_MATRIX_CORRECTED_1)
}



assert_that(!is.null(methylation_matrix_tumor))

if (!is.null(methylation_matrix_normal)) {
    logging("Normal sampel exists")
    common_ids <- intersect(colnames(methylation_matrix_tumor),
                            colnames(methylation_matrix_normal))
    methylation_matrix_corrected <- rbind(methylation_matrix_tumor[,colnames(methylation_matrix_tumor) %in% common_ids],
                                          methylation_matrix_normal[,colnames(methylation_matrix_normal) %in% common_ids])
    assert_that(length(common_ids) != 0)
    methylation_matrix_corrected <- methylation_matrix_corrected[order(rownames(methylation_matrix_corrected)),]
} else {
    logging("Normal sampel not exist")
    methylation_matrix_corrected <- methylation_matrix_tumor
}

logging("Write to local")
write_rds(methylation_matrix_corrected,
          opt$METHYLATION_MATRIX_CORRECTED,
          compress = "gz")

