
# Correct cnv data

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

# Read in options from command line
option_list = list(
    make_option(c("--CANCER"), type = "character",
                default = "LIHC",
                help="full path to the cnv data", metavar = "character"),
    make_option(c("--CNV_MATRIX_RAW"), type = "character",
                default = "/Volumes/USB128/MoR/cnv_preprocess_output/LIHC_cnv_matrix_gene.rds",
                help="full path to the cnv data", metavar = "character"),
    make_option(c("--CLINICAL_DATA"), type = "character",
                default = "/Volumes/USB128/MoR/clinical_preprocess_output/LIHC_clinical_data.rds",
                help="full path to the clinical data", metavar = "character"),

    make_option(c("--CNV_MATRIX_CORRECTED"), type = "character",
                default = "/Volumes/USB128/MoR/cnv_preprocess_output/LIHC_cnv_matrix_gene_corrected.rds",
                help="full path to the tumor clinic file", metavar = "character"),

    make_option(c("--HELPER"), type = "character",
                default = "/Users/MJ/Desktop/helper/correct_data.R",
                help="full path to the tumor clinic file", metavar = "character"),

    make_option(c("--TEMP_DIR"), type = "character",
                default = "/Users/MJ/Desktop/temp/plot/",
                help = "full path to the tumor clinic file", metavar = "character"),

    make_option(c("--NUM_CORES"), type = "integer",
                default = 2,
                help="", metavar = "character")
)

print(opt <- parse_args(OptionParser(option_list=option_list)))

registerDoParallel(opt$NUM_CORES)
source(opt$HELPER)

cnv_matrix_raw <- read_rds(opt$CNV_MATRIX_RAW)

assert_that(!any(is.na(cnv_matrix_raw)))
clinical_data <- read_rds(opt$CLINICAL_DATA)

logging("Processing normal samples")

cnv_matrix_normal <- correct_data(data_matrix = cnv_matrix_raw,
                                  clinical_data = clinical_data,
                                  data_set = 2,
                                  result_analysis = T,
                                  temp_file_prefix = str_c(opt$CANCER, "_cnv_normal_"),
                                  temp_dir = opt$TEMP_DIR,
                                  do_peer = F,
                                  do_combat = F,
                                  do_lm_clinical = T)


logging("Processing tumor samples")
cnv_matrix_tumor <- correct_data(data_matrix = cnv_matrix_raw,
                                 clinical_data = clinical_data,
                                 data_set = 1,
                                 result_analysis = T,
                                 temp_file_prefix = str_c(opt$CANCER, "_cnv_tumor_"),
                                 temp_dir = opt$TEMP_DIR,
                                 do_peer = F,
                                 do_combat = F,
                                 do_lm_clinical = T)


assert_that(!is.null(cnv_matrix_tumor))

if (!is.null(cnv_matrix_normal)) {
    common_ids <- intersect(colnames(cnv_matrix_tumor),
                            colnames(cnv_matrix_normal))
    cnv_matrix_corrected <- rbind(cnv_matrix_tumor[,colnames(cnv_matrix_tumor) %in% common_ids],
                                  cnv_matrix_normal[,colnames(cnv_matrix_normal) %in% common_ids])
    assert_that(length(common_ids) != 0)
    cnv_matrix_corrected <- cnv_matrix_corrected[order(rownames(cnv_matrix_corrected)),]
} else {
    cnv_matrix_corrected <- cnv_matrix_tumor
}

write_rds(cnv_matrix_corrected,
          opt$CNV_MATRIX_CORRECTED,
          compress = "gz")


