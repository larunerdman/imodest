
########
#   Generate miRNA isoform expression matrix

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

writeLines("10 - Input Specification")
###############################################################################

option_list <- list(
    make_option(c("--CANCER"), type = "character",
                default = "CESC",
                help = "", metavar = "character"),

    make_option(c("--MIRNA_EXPRESSION_DATA"), type = "character",
                default = "/Volumes/USB128/MoR/miRNA_preprocess_input/CESC.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_isoform_expression__data.data.txt",
                help = "", metavar = "character"),

    make_option(c("--CLINICAL_DATA"), type = "character",
                default = "/Volumes/USB128/MoR/clinical_preprocess_output/CESC_clinical_data.rds",
                help = "full path to the clinical data", metavar = "character"),

    make_option(c("--TEMP_DIR"), type = "character",
                default = "/Users/MJ/Desktop/temp/plot/",
                help = "full path to the tumor clinic file", metavar = "character"),
    make_option(c("--MIRNA_EXPRESSION_MATRIX_CORRECTED"), type = "character",
                default = "/Volumes/USB128/MoR/miRNA_preprocess_output/CESC_miRNA_matrix_corrected.rds",
                help = "full path to the tumor clinic file", metavar = "character"),

    make_option(c("--PEER_HIDDENS_SAMPLES_RATIO"), type = "double",
                default = 0.25,
                help = "", metavar = "character"),

    make_option(c("--HELPER"), type = "character",
                default = "/Users/MJ/Desktop/helper/correct_data.R",
                help = "full path to the tumor clinic file", metavar = "character"),

    make_option(c("--NUM_CORES"), type = "integer",
                default = 1,
                help = "", metavar = "character")
)

print(opt <- parse_args(OptionParser(option_list=option_list)))

registerDoParallel(opt$NUM_CORES)
source(opt$HELPER)

####
writeLines("30 - preprocess miRNAseq data")
################################################################################
# Column name setting
miRNA_isoform_expression <- read_tsv(opt$MIRNA_EXPRESSION_DATA,
                                     col_types = cols(
                                         SampleId = col_character(),
                                         miRNA_ID = col_character(),
                                         isoform_coords = col_character(),
                                         read_count = col_integer(),
                                         reads_per_million_miRNA_mapped = col_double(),
                                         `cross-mapped` = col_character(),
                                         miRNA_region = col_character()
                                     )) %>%
    tidyr::separate(miRNA_region, c("miRNA_region", "miRNA_id"), sep = ",") %>%
    select(sample_id = 1, miRNA_id = 8, RPM = 5) %>%
    filter(!is.na(miRNA_id)) %>%
    mutate(sample_id = str_extract(sample_id, "TCGA-..-....-...") %>% str_replace_all("-", ".")) %>%
    group_by(sample_id, miRNA_id) %>%
    summarise(sum_RPM = sum(RPM)) %>% ungroup() %>%
    tidyr::spread(miRNA_id, sum_RPM, fill = 0, convert = T)

# Turn expression data into matrix
miRNA_expression_matrix_raw <- select(miRNA_isoform_expression, -1) %>% as.matrix.data.frame()
rownames(miRNA_expression_matrix_raw) <- miRNA_isoform_expression$sample_id

miRNA_expression_matrix_raw <- miRNA_expression_matrix_raw[order(rownames(miRNA_expression_matrix_raw)), ]
miRNA_expression_matrix_raw <- miRNA_expression_matrix_raw[,order(colnames(miRNA_expression_matrix_raw))]


clinical_data <- read_rds(opt$CLINICAL_DATA)

logging("Processing normal samples")
miRNA_expression_matrix_normal <- correct_data(data_matrix = miRNA_expression_matrix_raw,
                                               clinical_data = clinical_data,
                                               data_set = 2,
                                               result_analysis = T,
                                               temp_file_prefix = str_c(opt$CANCER, "_miRNA_expression_normal_"),
                                               temp_dir = opt$TEMP_DIR,
                                               peer_hidden_samples_ratio = opt$PEER_HIDDENS_SAMPLES_RATIO)

logging("Processing tumor samples")
miRNA_expression_matrix_tumor <- correct_data(data_matrix = miRNA_expression_matrix_raw,
                                              clinical_data = clinical_data,
                                              data_set = 1,
                                              result_analysis = T,
                                              temp_file_prefix = str_c(opt$CANCER, "_miRNA_expression_tumor_"),
                                              temp_dir = opt$TEMP_DIR,
                                              peer_hidden_samples_ratio = opt$PEER_HIDDENS_SAMPLES_RATIO)


assert_that(!is.null(miRNA_expression_matrix_tumor))

if (!is.null(miRNA_expression_matrix_normal)) {
    common_ids <- intersect(colnames(miRNA_expression_matrix_tumor),
                            colnames(miRNA_expression_matrix_normal))
    miRNA_expression_matrix_corrected <- rbind(miRNA_expression_matrix_tumor[,colnames(miRNA_expression_matrix_tumor) %in% common_ids],
                                               miRNA_expression_matrix_normal[,colnames(miRNA_expression_matrix_normal) %in% common_ids])
    assert_that(length(common_ids) != 0)
    miRNA_expression_matrix_corrected <- miRNA_expression_matrix_corrected[order(rownames(miRNA_expression_matrix_corrected)),]
} else {
    miRNA_expression_matrix_corrected <- miRNA_expression_matrix_tumor
}

write_rds(miRNA_expression_matrix_corrected,
          opt$MIRNA_EXPRESSION_MATRIX_CORRECTED,
          compress = "gz")

