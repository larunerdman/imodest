

########
# Preprocess metyhalation data from firebrowse and generaete 2d methylation matrix
#   output format: row x col = sample_id x methylation probes


library(readr)
library(optparse)
library(stringr)
library(dplyr)

option_list <- list(
    make_option(c("--CANCER_TYPE"),
                type = "character",
                default = "OV",
                help = "Cancer set abbreviation name", metavar = "character"),
    make_option(c("--METHYLATION_INPUT_FOLDER"),
                type = "character",
                default = "/Users/MJ/Desktop/",
                help = "CNV", metavar = "character"),
    make_option(c("--METHYLATION_OUTPUT_FOLDER"),
                type = "character",
                default = "/Volumes/USB128/methylation_preprocess_output/",
                help = "", metavar = "character")
);

opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

meyhylation_raw_file <- str_c(opt$METHYLATION_INPUT_FOLDER, opt$CANCER_TYPE,".methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt")
Sample_ids <- read_lines(meyhylation_raw_file, n_max = 1) %>%
    str_split(pattern = "\\t", simplify = T) %>%
    str_replace_all("-", "\\.") %>%
    str_extract("^TCGA............") %>%
    str_subset("^TCGA.+") %>%
    unique()


####
#f <- function(x, pos) x %>% select(c(methylation_probe_id = 1, seq(from = 2, length(Sample_ids) * 4 + 1, 4)))
methylation_data <- read_tsv(meyhylation_raw_file, skip = 2, col_names = F) %>% select(c(methylation_probe_id = 1, seq(from = 2, length(Sample_ids) * 4 + 1, 4)))

methylation_data_2d <- methylation_data %>% select(-methylation_probe_id) %>% as.matrix.data.frame()
rownames(methylation_data_2d) <- methylation_data$methylation_probe_id
colnames(methylation_data_2d) <- Sample_ids
methylation_data_2d <- methylation_data_2d[complete.cases(methylation_data_2d),]

methylation_data_2d <- methylation_data_2d[order(rownames(methylation_data_2d)),]
methylation_data_2d <- methylation_data_2d[,order(colnames(methylation_data_2d))]

write_rds(t(methylation_data_2d), str_c(opt$METHYLATION_OUTPUT_FOLDER, opt$CANCER_TYPE, "_methylation_matrix.rds"), compress = "gz")
