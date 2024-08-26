library(readr)
library(dplyr, warn.conflicts = F)
library(stringr)

setwd("/Volumes/USB128/SNP_preprocess_output")
# BRCA-pruned-snps.txt* KIRC-pruned-snps.txt* LGG-pruned-snps.txt*  LUAD-pruned-snps.txt* THCA-pruned-snps.txt*

# a) Cancer type (ex. LGG, BRCA, KIRC, LUAD, THCA)
CANCER_TYPE <- "LGG"

# b) SNP pruned data
SNP_DATA <- str_c(CANCER_TYPE, "-pruned-snps.txt")

# c) sample id from isoform expression data
SAMPLE_ID <- colnames(read_rds(str_c("../RNA_preprocess_output/", CANCER_TYPE,"_RNASeq_isoform.rds")))

###############################################################################
# 2) Output Specification
###############################################################################
writeLines("Step 2 - Output Specification")
# a) Genome-wide model:
#       2D matrix -- row: SNP probe, col: Sample ID
OUTPUT <- str_c(CANCER_TYPE, "_SNP_expression.rds")

###############################################################################
# 3) Prepross SNP data
###############################################################################
writeLines("Step 3  - Preprocess SNP data")
SNP_data <- read_lines(SNP_DATA) %>% str_trim() %>% str_split("\t")

# Collect sample Id information
num_sample <- length(SNP_data) - 1
sample_ID <- rep(as.character(NA), num_sample)
for (i in (seq_len(num_sample))) {
    sample_ID[i] <- SNP_data[[i+1]][1]
}

# Simplify sample ID formatting
if (CANCER_TYPE == "LUAD") {
    sample_ID <- sample_ID %>%
    str_extract("TCGA-..-....-...") %>%
    str_replace_all("-", ".")
}


# Filter missing sample ID
SNP_data <- SNP_data[c(TRUE,!is.na(sample_ID))]
sample_ID <- sample_ID[!is.na(sample_ID)]
num_sample <- length(sample_ID)

# Collect SNP ID
num_SNP <- length(SNP_data[[1]]) - 6
Old_num <- num_SNP
SNP_ID <- SNP_data[[1]][7:(num_SNP + 6)]

# Transform raw data into 2D matrix
SNP_data_2D <- matrix(as.integer(NA), num_SNP, num_sample, dimnames = list(SNP_ID, sample_ID))

# Collect SNP value, starting from second row, and starting from the 4th column
for (i in seq_len(num_sample)) {
    SNP_value <- SNP_data[[i+1]]
    SNP_data_2D[,i] <- type.convert(SNP_value[7:(num_SNP + 6)], na.strings = c(""))
}

SNP_data_2D <- SNP_data_2D[,colnames(SNP_data_2D) %in% SAMPLE_ID]
SNP_data_2D <- SNP_data_2D[, order(colnames(SNP_data_2D))]
SNP_data_2D <- SNP_data_2D[order(rownames(SNP_data_2D)),]
SNP_data_2D_final <- SNP_data_2D[complete.cases(SNP_data_2D),]

###############################################################################
# 4) Output SNP data
###############################################################################
print("Step 4 - Output SNP data")
write_rds(t(SNP_data_2D_final), OUTPUT)
