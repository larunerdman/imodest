
# Correct methylation raw matrix
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
                default = "LUAD",
                help="full path to the methylation data", metavar = "character"),
    make_option(c("--METHYLATION_MATRIX_RAW"), type = "character",
                default = "/Volumes/USB128/MoR/methylation_preprocess_output/LUAD_methylation_matrix.rds",
                help="full path to the methylation data", metavar = "character"),
    make_option(c("--CLINICAL_DATA"), type = "character",
                default = "/Volumes/USB128/MoR/clinical_preprocess_output/LUAD_clinical_data.rds",
                help="full path to the clinical data", metavar = "character"),

    make_option(c("--METHYLATION_MATRIX_CORRECTED"), type = "character",
                default = "/Volumes/USB128/MoR/methylation_preprocess_output/LUAD_methylation_matrix_corrected.rds",
                help="full path to the tumor clinic file", metavar = "character"),

    make_option(c("--TEMP_DIR"), type = "character",
                default = "/Users/MJ/Desktop/temp/plot/",
                help = "full path to the tumor clinic file", metavar = "character"),

    make_option(c("--DATA_SET"), type = "integer",
                default = 2,
                help = "", metavar = "character"),

    make_option(c("--PEER_HIDDENS_SAMPLES_RATIO"), type = "double",
                default = 0,
                help = "", metavar = "character"),
    make_option(c("--HELPER"), type = "character",
                default = "/Users/MJ/Desktop/helper/correct_data.R",
                help="full path to the tumor clinic file", metavar = "character"),
    make_option(c("--NUM_CORES"), type = "integer",
                default = 1,
                help="", metavar = "character")
)

print(opt <- parse_args(OptionParser(option_list=option_list)))

registerDoParallel(opt$NUM_CORES)
source(opt$HELPER)

methylation_matrix_raw <- read_rds(opt$METHYLATION_MATRIX_RAW)

clinical_data <- read_rds(opt$CLINICAL_DATA)

if (opt$DATA_SET == 2) {
    logging("Processing normal samples")
    methylation_matrix_normal <- correct_data(data_matrix = methylation_matrix_raw,
                                              clinical_data = clinical_data,
                                              data_set = opt$DATA_SET,
                                              result_analysis = F,
                                              do_logit = T,
                                              temp_file_prefix = str_c(opt$CANCER, "_methylation_normal_"),
                                              temp_dir = opt$TEMP_DIR,
                                              peer_hidden_samples_ratio = opt$PEER_HIDDENS_SAMPLES_RATIO)
    if (!is.null(methylation_matrix_normal)) {
        write_rds(methylation_matrix_normal,
                  opt$METHYLATION_MATRIX_CORRECTED,
                  compress = "gz")
    }

} else {
    logging("Processing tumor samples")
    methylation_matrix_tumor <- correct_data(data_matrix = methylation_matrix_raw,
                                             clinical_data = clinical_data,
                                             data_set = opt$DATA_SET,
                                             result_analysis = F,
                                             do_logit = T,
                                             temp_file_prefix = str_c(opt$CANCER, "_methylation_tumor_"),
                                             temp_dir = opt$TEMP_DIR,
                                             peer_hidden_samples_ratio = opt$PEER_HIDDENS_SAMPLES_RATIO)
    if (!is.null(methylation_matrix_tumor)) {
        write_rds(methylation_matrix_tumor,
                  opt$METHYLATION_MATRIX_CORRECTED,
                  compress = "gz")
    }

}



