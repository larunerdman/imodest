
####
#   Generate matrix of dosage file
#   Output file format:
#       Size: num_sample x num_SNP

library(methods)
library(optparse) # make_option
library(stringr)  # str_replace, str_c
library(assertthat)
library(readr)
library(dplyr)
library(tidyr)
library(Matrix)

option_list = list(
    make_option(c("--BIM_FILE"), type = "character", default = "/Users/MJ/Desktop/temp/KIRC-cancer-postqc-postqc2-postqc-somup-chr22.bim",
                help = "", metavar = "character"),
    make_option(c("--DOSAGE_FILE"), type = "character", default = "/Users/MJ/Desktop/temp/KIRC-cancer-postqc-postqc2-postqc-somup-chr22.dosage",
                help = "", metavar = "character"),
    make_option(c("--SAMPLE_FILE"), type = "character", default = "/Users/MJ/Desktop/temp/KIRC.sample",
                help = "", metavar = "character"),
    make_option(c("--SNP_MATRIX"), type = "character", default = "/Users/MJ/Desktop/temp/res.rds",
                help = "", metavar = "character")
)


print(opt <- parse_args(OptionParser(option_list = option_list)))

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
writeLines("####\n 10 Preprocess SNP sample id")
sample_ids <- read_tsv(opt$SAMPLE_FILE, col_names = F, col_types = cols(X1 = col_character()))$X1

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
writeLines("####\n 20 Preprocess SNP summary")

SNP_data <- read_delim(opt$BIM_FILE, col_names = F, delim = " ",
                       col_types = cols(
                           X1 = col_integer(),
                           X2 = col_character(),
                           X3 = col_integer(),
                           X4 = col_integer(),
                           X5 = col_character(),
                           X6 = col_character()
                       )) %>%
    select(chr = 1, coord = 3, major = 5, minor = 6) %>%
    mutate(snp_id = str_c("chr", chr, "_", coord, "_", major, "_", minor))

total_size <- nrow(SNP_data)
max_size <- 8000
skip <- 0
num_rounds <- total_size %/% max_size + 1

for (i in seq_len(num_rounds)){
    writeLines(str_c("skip: ", skip))
    SNP_data <- read_delim(opt$BIM_FILE, col_names = F, delim = " ",
                           col_types = cols(
                               X1 = col_integer(),
                               X2 = col_character(),
                               X3 = col_integer(),
                               X4 = col_integer(),
                               X5 = col_character(),
                               X6 = col_character()
                           ),
                           skip = skip,
                           n_max = max_size) %>%
        select(chr = 1, coord = 3, major = 5, minor = 6) %>%
        mutate(snp_id = str_c("chr", chr, "_", coord, "_", major, "_", minor))

    SNP_dosage <- read_delim(opt$DOSAGE_FILE, delim = " ", col_names = F,
                             col_types = cols(
                                 .default = col_double()
                             ),
                             skip = skip,
                             n_max = max_size)

    if (any(duplicated(sample_ids))) {
        SNP_dosage <- SNP_dosage %>% select(-which(duplicated(sample_ids)))
        colnames(SNP_dosage) <- sample_ids[-which(duplicated(sample_ids))]
    } else {
        colnames(SNP_dosage) <- sample_ids
    }


    SNP_dosage_2d <- SNP_dosage %>% as.matrix() %>% t()
    colnames(SNP_dosage_2d) <- SNP_data$snp_id

    SNP_dosage_2d <- SNP_dosage_2d[order(rownames(SNP_dosage_2d)),]
    SNP_dosage_2d <- SNP_dosage_2d[,order(colnames(SNP_dosage_2d))]

    assert_that(!is.unsorted(rownames(SNP_dosage_2d)))
    assert_that(!any(duplicated(rownames(SNP_dosage_2d))))
    assert_that(!is.unsorted(colnames(SNP_dosage_2d)))
    assert_that(!any(duplicated(colnames(SNP_dosage_2d))))

    output_name <- str_replace(opt$SNP_MATRIX, ".rds$", str_c("_", i,".rds"))
    write_rds(SNP_dosage_2d, output_name, compress = "gz")

    skip <- skip + max_size
}

for (i in seq_len(num_rounds)) {
    output_name <- str_replace(opt$SNP_MATRIX, ".rds$", str_c("_", i,".rds"))
    if (i == 1) {
        SNP_dosage_2d <- read_rds(output_name)
    } else {
        SNP_dosage_2d <- cbind(SNP_dosage_2d, read_rds(output_name))
    }
}


SNP_dosage_2d <- SNP_dosage_2d[,order(colnames(SNP_dosage_2d))]
assert_that(!is.unsorted(rownames(SNP_dosage_2d)))
assert_that(!any(duplicated(rownames(SNP_dosage_2d))))
assert_that(!is.unsorted(colnames(SNP_dosage_2d)))
assert_that(!any(duplicated(colnames(SNP_dosage_2d))))

write_rds(SNP_dosage_2d, opt$SNP_MATRIX, compress = "gz")
