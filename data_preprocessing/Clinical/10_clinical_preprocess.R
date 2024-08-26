
#### clinical data preprocessing

# Utill
library(methods)
library(lazyeval)
library(gtools)
library(optparse)
library(compiler)
library(stringr)
library(assertthat, warn.conflicts = F)

# tibble minipulation
library(readr)
library(tibble,warn.conflicts = F)
library(dplyr, warn.conflicts = F)


# Parallel
library(iterators)
library(parallel)
library(foreach)
library(doParallel)


####################################################
writeLines("####\n 10 - Input spcfication")
####################################################

option_list <- list(
    make_option(c("--CANCER_TYPE"), type = "character", default="LIHC",
                help="cancer type", metavar = "character"),

    make_option(c("--PICKED_CLINICAL_DATA"), type = "character",
                default="/Volumes/USB128/MoR/clinical_data_input/LIHC.clin.picked.txt",
                help = "clinical data", metavar = "character"),
    make_option(c("--SNP_PCA"), type = "character",
                default="/Volumes/USB128/MoR/clinical_data_input/all-TCGA-postqc-ur-ols-rmvd-pruned-PCS.txt",
                help = "SNP pc", metavar = "character"),

    make_option(c("--CLINICAL_DATA"), type = "character",
                default = "/Volumes/USB128/MoR/clinical_preprocess_output/LIHC_clinical_data.rds",
                help = "", metavar = "character"),
    make_option(c("--GENERAL_CLINICAL_MATRIX"), type = "character",
                default = "/Volumes/USB128/MoR/clinical_preprocess_output/LIHC_general_clinical_matrix.rds",
                help = "", metavar = "character"),
    make_option(c("--TUMOR_CLINICAL_MATRIX"), type = "character",
                default = "/Volumes/USB128/MoR/clinical_preprocess_output/LIHC_tumor_clinical_matrix.rds",
                help = "", metavar = "character"),
    make_option(c("--NORMAL_CLINICAL_MATRIX"), type = "character",
                default = "/Volumes/USB128/MoR/clinical_preprocess_output/LIHC_tumor_clinical_matrix.rds",
                help = "", metavar = "character")
)

print(opt <- parse_args(OptionParser(option_list=option_list)))

registerDoParallel(opt$NUM_CORES)

####################################################
writeLines("####\n 10 - Read in clincial data to construct covaraite matrix")
####################################################
clinical <- as_tibble(t(read_tsv(opt$PICKED_CLINICAL_DATA, col_names = F,
                                 col_types = cols(
                                     .default = col_character()
                                 ))))

colnames(clinical) <- clinical[1,]
clinical <- clinical[-1,]

clinical_other <- clinical[,!str_detect(colnames(clinical),"sample.")]
clinical_sample_1 <- bind_cols(clinical[,str_detect(colnames(clinical),"sample\\.")], clinical_other)
clinical_sample_2 <- bind_cols(clinical[,str_detect(colnames(clinical),"sample-2")], clinical_other)
clinical_sample_3 <- bind_cols(clinical[,str_detect(colnames(clinical),"sample-3")], clinical_other)

colnames(clinical_sample_2) <- str_replace(colnames(clinical_sample_2), "sample-2", "sample")
clinical_full <- bind_rows(clinical_sample_1, clinical_sample_2)

colnames(clinical_sample_3) <- str_replace(colnames(clinical_sample_3), "sample-3", "sample")
clinical_full <- bind_rows(clinical_full, clinical_sample_3)



####################################################
writeLines("####\n 20 - Preprocess bcr_sample_barcode")
####################################################

column_bcr_sample_barcode <- str_which(colnames(clinical_full), "sample.bcr_sample_barcode")
assert_that(length(column_bcr_sample_barcode) == 1)

colnames(clinical_full)[column_bcr_sample_barcode] <- "sample_id"

clinical_full <- clinical_full %>% mutate(sample_id = str_to_upper(sample_id) %>%
                                              str_replace_all("-", "\\.") %>%
                                              str_extract("^TCGA.{12}"))
####################################################
writeLines("####\n 30 - Preprocess gender")
####################################################
column_gender <- str_which(colnames(clinical_full), "gender")
assert_that(length(column_gender) == 1)
colnames(clinical_full)[column_gender] <- "gender"



####################################################
writeLines("####\n 40 - Preprocess age")
####################################################
column_age <- str_which(colnames(clinical_full), "age_at_initial_pathologic_diagnosis")
assert_that(length(column_age) == 1)
colnames(clinical_full)[column_age] <- "age_at_initial_pathologic_diagnosis"
clinical_full <- clinical_full %>%
    mutate(age_at_initial_pathologic_diagnosis = age_at_initial_pathologic_diagnosis %>% as.integer())


####################################################
writeLines("####\n 50 - Preprocess tumor_nuclei_percentage")
####################################################
column_tumor_nuclei_percent <- str_which(colnames(clinical_full), "tumor_nuclei_percent")

if (length(column_tumor_nuclei_percent) == 1) {
    colnames(clinical_full)[column_tumor_nuclei_percent] <- "tumor_nuclei_percent"
    clinical_full <- clinical_full %>% mutate(tumor_nuclei_percent = tumor_nuclei_percent %>% as.integer())
    clinical_full$tumor_nuclei_percent[str_which(clinical_full$sample_id, "1..$")] <- 0
}

####################################################
writeLines("####\n 60 - Preprocess histological_type")
####################################################
column_histological_type <- str_which(colnames(clinical_full), "histological_type$")

if (length(column_histological_type) == 1) {
    colnames(clinical_full)[column_histological_type] <- "histological_type"
    clinical_full <- clinical_full %>% mutate(histological_type = histological_type %>%
                                                  str_replace_na("") %>%
                                                  str_replace("mixed histology \\(please specify\\)", "") %>%
                                                  str_replace("other, specify", ""))
}

column_histological_type_other <- str_which(colnames(clinical_full), "histological_type_other$")

if (length(column_histological_type_other) == 1) {
    colnames(clinical_full)[column_histological_type_other] <- "histological_type_other"
    clinical_full <- clinical_full %>% mutate(histological_type = str_c(histological_type, str_replace_na(histological_type_other, ""))) %>%
        select(-histological_type_other)
}


####################################################
writeLines("####\n 70 - Preprocess site of disease")
####################################################
column_site_of_disease <- str_which(colnames(clinical_full), "site_of_disease")

if (length(column_site_of_disease) == 1) {
    colnames(clinical_full)[column_site_of_disease] <- "site_of_disease"
    clinical_full <- clinical_full %>% mutate(site_of_disease = site_of_disease %>% str_replace_na(""))
}


column_site_of_disease_description <- str_which(colnames(clinical_full), "site_of_disease_description")
if (length(column_site_of_disease_description) == 1) {
    colnames(clinical_full)[column_site_of_disease_description] <- "site_of_disease_description"
    clinical_full <- clinical_full %>%
        mutate(site_of_disease = str_c(site_of_disease, str_replace_na(site_of_disease_description, ""))) %>%
        select(-site_of_disease_description)
}


####################################################
writeLines("####\n 80 - Preprocess admin.batch_number")
####################################################
column_batch_number <- str_which(colnames(clinical_full), "batch_number")
if (length(column_batch_number) == 1) {
    colnames(clinical_full)[column_batch_number] <- "batch_number"
}

####################################################
writeLines("####\n 90 - Preprocess patient.icd_10")
####################################################
column_icd_10 <- str_which(colnames(clinical_full), "icd_10")
if (length(column_icd_10) == 1) {
    colnames(clinical_full)[column_icd_10] <- "icd_10"
}

####################################################
writeLines("####\n 100 - Preprocess patient.icd_o_3")
####################################################
column_icd_o_3_site <- str_which(colnames(clinical_full), "icd_o_3_site")
if (length(column_icd_o_3_site) == 1) {
    colnames(clinical_full)[column_icd_o_3_site] <- "icd_o_3_site"
}

column_icd_o_3_histology <- str_which(colnames(clinical_full), "icd_o_3_histology")
if (length(column_icd_o_3_histology) == 1) {
    colnames(clinical_full)[column_icd_o_3_histology] <- "icd_o_3_histology"
}

# #### 90 - Annotate sample types
# clinical_full <- clinical_full %>%
#     mutate(normal_sample = str_detect(sample_id, "1..$"),
#            tumor_sample = str_detect(sample_id, "0..$"))

####################################################
writeLines("####\n 110 - Preprocess SNP pc data")
####################################################
snp_pc <- read_delim(opt$SNP_PCA, delim = " ",
                     col_types = cols(
                         id = col_character(),
                         PC1 = col_double(),
                         PC2 = col_double()
                     )) %>%
    mutate(patient_id = id %>% str_sub(1, 12)) %>%
    select(-id)


####################################################
writeLines("####\n 120 - Annotate clinical data with SNP pc values")
####################################################
clinical_final <- clinical_full[complete.cases(clinical_full),] %>%
    arrange(sample_id) %>%
    mutate(patient_id = sample_id %>% str_extract("^TCGA.{8}")) %>%
    inner_join(snp_pc, by = c("patient_id")) %>%
    select(-patient_id) %>%
    arrange(sample_id)

write_rds(clinical_final, opt$CLINICAL_DATA, compress = "gz")


####################################################
writeLines("####\n 130 - Convert each categorical varible to indicator varible")
####################################################


covariate_summary <- clinical_final[,"sample_id"]
covariates_names <- c("sample_id")
COVARIATES <- c("sample_id")
for (i in seq_len(ncol(clinical_final))[-column_bcr_sample_barcode]) {
    varibles <- clinical_final[[i]]
    unique_varibles <- unique(varibles)
    unique_varibles <- unique_varibles[unique_varibles != "" | is.na(unique_varibles)]
    num_unique_varibles <- length(unique_varibles)

    if (!all(is.numeric(varibles)) && !all(is.logical(varibles))) {
        for (j in seq_along(unique_varibles)) {
            print(unique_varibles[j])
            print(sum(varibles == unique_varibles[j]))
            covariate_summary <- covariate_summary %>%
                bind_cols(tibble(varibles == unique_varibles[j]))

            covariates_names <- c(covariates_names,
                                  str_c(colnames(clinical_final)[i], j, sep = "_")) #unique_varibles[j], sep = "_"))
            COVARIATES <- c(COVARIATES, unique_varibles[j])
        }

    } else {
        covariate_summary <- covariate_summary %>% bind_cols(clinical_final[,i])
        covariates_names <- c(covariates_names, colnames(clinical_final)[i])
        COVARIATES <- c(COVARIATES, colnames(clinical_final)[i])
    }

}

# covariates_names <- covariates_names %>%
#     str_replace_all(" ", ".") %>%
#     str_replace_all("\\/", ".") %>%
#     str_replace_all("\\(", "") %>%
#     str_replace_all("\\)", "")

colnames(covariate_summary) <- covariates_names

covariate_summary[covariate_summary == T] <- 1

####################################################
writeLines("####\n 140 - Convert clinical data into 2d matrix")
####################################################
covariate_matrix <- covariate_summary %>% select(-sample_id) %>% as.matrix.data.frame()
rownames(covariate_matrix) <- covariate_summary$sample_id

covariate_matrix <- covariate_matrix[order(rownames(covariate_matrix)),]
covariate_matrix <- covariate_matrix[,order(colnames(covariate_matrix))]


####################################################
writeLines("####\n 150 Split feature into tumor nad normal sample")
####################################################
general_clinical_features <- c("gender",
                               "age_at_initial_pathologic_diagnosis",
                               "batch_number",
                               "tumor_nuclei_percent",
                               "PC1",
                               "PC2")
tumor_covariate_matrix <- covariate_matrix[str_detect(rownames(covariate_matrix), "0..$"),
                                           !str_detect(colnames(covariate_matrix),
                                                       str_c(general_clinical_features, collapse = "|"))]

normal_covariate_matrix <- covariate_matrix[str_detect(rownames(covariate_matrix), "1..$"),
                                           str_detect(colnames(covariate_matrix),
                                                       str_c(general_clinical_features, collapse = "|"))]

general_covariate_matrix <- covariate_matrix[,str_detect(colnames(covariate_matrix),
                                                         str_c(general_clinical_features, collapse = "|"))]

general_covariate_matrix <- cbind(general_covariate_matrix,
                                  tumor_sample = str_detect(rownames(covariate_matrix), "0..$"),
                                  normal_sample = str_detect(rownames(covariate_matrix), "1..$"))


####################################################
writeLines("####\n 160 Remove constant covariate column")
####################################################
tumor_covariate_matrix <- tumor_covariate_matrix[,!near(apply(tumor_covariate_matrix, 2, var),0)]
general_covariate_matrix <- general_covariate_matrix[,!near(apply(general_covariate_matrix, 2, var),0)]
normal_covariate_matrix <- normal_covariate_matrix[,!near(apply(normal_covariate_matrix, 2, var),0)]

####################################################
writeLines("####\n 170 Write to local")
####################################################
write_rds(general_covariate_matrix, opt$GENERAL_CLINICAL_MATRIX, compress = "gz")
write_rds(tumor_covariate_matrix, opt$TUMOR_CLINICAL_MATRIX, compress = "gz")
write_rds(normal_covariate_matrix, opt$NORMAL_CLINICAL_MATRIX, compress = "gz")
