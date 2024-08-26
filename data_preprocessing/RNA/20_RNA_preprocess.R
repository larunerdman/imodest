
#### Collec gene expression data and clinical data to be processed by data correction

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



####################################################
writeLines("####\n 10 - Input spcfication")
####################################################

option_list <- list(
    make_option(c("--CANCER"), type = "character", default="LIHC",
                help = "cancer type", metavar = "character"),

    # Input
    make_option(c("--CLINICAL_DATA_RAW"), type = "character",
                default="/Volumes/USB128/MoR/clinical_data_input/LIHC.clin.merged.txt",
                help = "", metavar = "character"),
    make_option(c("--RNA_INPUT_FOLDER"), type = "character",
                default = "/Volumes/USB128/MoR/RNA_preprocess_input/LIHC/",
                help = "", metavar = "character"),

    make_option(c("--CLINICAL_DATA"), type = "character",
                default = "/Volumes/USB128/MoR/clinical_preprocess_output/LIHC_clinical_data.rds",
                help = "full path to the clinical data", metavar = "character"),
    make_option(c("--TF_TARGET_MATRIX"), type = "character",
                default = "/Volumes/USB128/ref/target_matrix/10_TF_target_matrix_gene.rds",
                help = "", metavar = "character"),

    make_option(c("--PC_GENE_EXPRESSION_CORRECTED"), type = "character",
                default = "/Volumes/USB128/MoR/gene_expression_preprocess_output/LIHC_gene_expression_matrix_corrected.rds",
                help = "full path to the tumor clinic file", metavar = "character"),
    make_option(c("--LNCRNA_GENE_EXPRESSION_CORRECTED"), type = "character",
                default = "/Volumes/USB128/MoR/gene_expression_preprocess_output/LIHC_gene_expression_matrix_corrected.rds",
                help = "full path to the tumor clinic file", metavar = "character"),
    make_option(c("--TF_GENE_EXPRESSION_CORRECTED"), type = "character",
                default = "/Volumes/USB128/MoR/gene_expression_preprocess_output/LIHC_gene_expression_matrix_corrected.rds",
                help = "full path to the tumor clinic file", metavar = "character"),

    make_option(c("--HELPER"), type = "character",
                default = "/Users/MJ/Desktop/helper/correct_data.R",
                help = "full path to the tumor clinic file", metavar = "character"),

    make_option(c("--TEMP_DIR"), type = "character",
                default = "/Users/MJ/Desktop/temp/plot/",
                help = "full path to the tumor clinic file", metavar = "character"),

    make_option(c("--PEER_HIDDENS_SAMPLES_RATIO"), type = "double",
                default = 0.25,
                help = "", metavar = "character"),
    make_option(c("--NUM_CORES"), type = "integer",
                default = 2,
                help = "", metavar = "character"),
    make_option(c("--COME_FROM_TAR"), type = "integer",
                default = 0,
                help = "", metavar = "character")
)

print(opt <- parse_args(OptionParser(option_list=option_list)))

registerDoParallel(opt$NUM_CORES)
source(opt$HELPER)

####################################################
writeLines("####\n 20 - Preprocess clincal datta to obtain TCGA sample id to uuid mapping")
####################################################

if (!opt$COME_FROM_TAR) {
    clinical_data <- read_lines(opt$CLINICAL_DATA_RAW) %>%
        str_subset("(\\.bcr_aliquot_barcode)|(\\.bcr_aliquot_uuid)") %>%
        str_subset("^patient.samples.sample") %>%
        str_subset("portions.portion") %>%
        str_subset("analytes.analyte") %>%
        str_subset("aliquots.aliquot") %>%
        str_split("\\t")

    for (i in seq_len(length(clinical_data) / 2)) {
        assert_that(clinical_data[[2*i - 1]][1] %>%
                        str_replace("bcr_aliquot_barcode", "bcr_aliquot_uuid") == clinical_data[[2 * i]][1])
        aliquot_i <- tibble(bcr_aliquot_barcode = clinical_data[[2 * i - 1]][-1],
                            bcr_aliquot_uuid = clinical_data[[2 * i]][-1])

        if (i == 1) {
            aliquot_info <- aliquot_i
        } else {
            aliquot_info <- bind_rows(aliquot_info, aliquot_i)
        }
    }

    aliquot_info <- aliquot_info %>%
        filter(bcr_aliquot_barcode != "NA" & bcr_aliquot_uuid != "NA") %>%
        mutate(bcr_aliquot_barcode = bcr_aliquot_barcode %>% str_to_upper() %>% str_replace_all("-", ".")) %>%
        mutate(sample_id = bcr_aliquot_barcode %>% str_extract("TCGA.{12}")) %>%
        arrange(bcr_aliquot_uuid)

    assert_that(!any(duplicated(aliquot_info$bcr_aliquot_barcode)))
    assert_that(!is.unsorted(aliquot_info$bcr_aliquot_uuid))
}

####################################################
writeLines("####\n 30 - Read in protein coding gene expression")
####################################################

SIZE_FILTER <- str_detect((list.files(opt$RNA_INPUT_FOLDER, full.names = T) %>% str_subset("pc_gene"))[1],"_\\d+_pc_gene_expression.txt")

if (SIZE_FILTER) {
    pc_gene_expression_files <- tibble(pc_gene_file = list.files(opt$RNA_INPUT_FOLDER, full.names = T) %>% str_subset("pc_gene")) %>%
        mutate(aliquot_uuid = pc_gene_file %>% basename() %>% str_replace("_\\d+_pc_gene_expression.txt", ""),
               size = pc_gene_file %>% basename() %>% str_extract("_\\d+_pc_gene_expression.txt") %>% str_extract("\\d+") %>% as.integer()) %>%
        group_by(aliquot_uuid) %>% mutate(max_size = max(size)) %>% ungroup() %>%
        filter(size == max_size) %>%
        select(-size, -max_size) %>%
        arrange(aliquot_uuid)

    lncRNA_gene_expression_files <- tibble(lncRNA_gene_file = list.files(opt$RNA_INPUT_FOLDER, full.names = T) %>% str_subset("lncRNA_gene")) %>%
        mutate(aliquot_uuid = lncRNA_gene_file %>% basename() %>% str_replace("_\\d+_lncRNA_gene_expression.txt", ""),
               size = lncRNA_gene_file %>% basename() %>% str_extract("_\\d+_lncRNA_gene_expression.txt") %>% str_extract("\\d+") %>% as.integer()) %>%
        group_by(aliquot_uuid) %>% mutate(max_size = max(size)) %>% ungroup() %>%
        filter(size == max_size) %>%
        select(-size, -max_size) %>%
        arrange(aliquot_uuid)

} else {
    pc_gene_expression_files <- tibble(pc_gene_file = list.files(opt$RNA_INPUT_FOLDER, full.names = T) %>% str_subset("pc_gene")) %>%
        mutate(aliquot_uuid = pc_gene_file %>% basename() %>% str_replace("_pc_gene_expression.txt", "")) %>%
        arrange(aliquot_uuid)

    lncRNA_gene_expression_files <- tibble(lncRNA_gene_file = list.files(opt$RNA_INPUT_FOLDER, full.names = T) %>% str_subset("lncRNA_gene")) %>%
        mutate(aliquot_uuid = lncRNA_gene_file %>% basename() %>% str_replace("_lncRNA_gene_expression.txt", "")) %>%
        arrange(aliquot_uuid)
}

assert_that(nrow(pc_gene_expression_files) == nrow(lncRNA_gene_expression_files))
assert_that(all(pc_gene_expression_files$aliquot_uuid == lncRNA_gene_expression_files$aliquot_uuid))
if (!opt$COME_FROM_TAR) {
    gene_expression_files <- pc_gene_expression_files %>%
        inner_join(lncRNA_gene_expression_files, by = c("aliquot_uuid")) %>%
        inner_join(aliquot_info, by = c("aliquot_uuid" = "bcr_aliquot_uuid"))

    if (nrow(gene_expression_files) == 0) {
        pc_gene_expression_files <- pc_gene_expression_files %>% mutate(sample_id = aliquot_uuid %>% str_to_upper() %>% str_replace_all("-", ".") %>% str_extract("TCGA.{12}"))
        lncRNA_gene_expression_files <- lncRNA_gene_expression_files %>% mutate(sample_id = aliquot_uuid %>% str_to_upper() %>% str_replace_all("-", ".") %>% str_extract("TCGA.{12}"))
        gene_expression_files <- pc_gene_expression_files %>%
            inner_join(lncRNA_gene_expression_files, by = c("sample_id"))
        gene_expression_files <- gene_expression_files[complete.cases(gene_expression_files),]
    } else {
        assert_that(!any(duplicated(gene_expression_files$aliquot_uuid)))
        assert_that(!any(duplicated(gene_expression_files$aliquot_uuid)))
    }
} else {
    gene_expression_files <- pc_gene_expression_files %>%
        inner_join(lncRNA_gene_expression_files, by = c("aliquot_uuid")) %>%
        select(1,sample_id=2,3) %>%
        mutate(sample_id = sample_id %>% str_replace_all("-", "."))
}


assert_that(nrow(gene_expression_files) != 0)


# Enforce unique aliquot for each sample
gene_expression_files <- gene_expression_files %>% filter(!duplicated(sample_id))
assert_that(!any(duplicated(gene_expression_files$sample_id)))

# Collect gene expression for each aliquot
for (i in seq_len(nrow(gene_expression_files))) {
    pc_gene_expression_i <- read_tsv(gene_expression_files$pc_gene_file[i],
                                     col_types = cols(
                                         Name = col_character(),
                                         Length = col_double(),
                                         EffectiveLength = col_double(),
                                         TPM = col_double(),
                                         NumReads = col_double()
                                     ), progress = F) %>%
        select(1, 4) %>%
        setNames(c("gene_id", gene_expression_files$sample_id[i]))

    lncRNA_gene_expression_i <- read_tsv(gene_expression_files$lncRNA_gene_file[i],
                                         col_types = cols(
                                             Name = col_character(),
                                             Length = col_double(),
                                             EffectiveLength = col_double(),
                                             TPM = col_double(),
                                             NumReads = col_double()
                                         ), progress = F) %>%
        select(1, 4) %>%
        setNames(c("gene_id", gene_expression_files$sample_id[i]))

    gene_expression_i <- bind_rows(pc_gene_expression_i, lncRNA_gene_expression_i)

    num_genes <- nrow(gene_expression_i)
    pc_gene_ids <- pc_gene_expression_i$gene_id
    lncRNA_gene_ids <- lncRNA_gene_expression_i$gene_id

    if (i == 1) {
        gene_expression <- gene_expression_i
        NUM_GENES <- num_genes
        PC_GENE_IDS <- pc_gene_ids
        LNCRNA_GENE_IDS <- lncRNA_gene_ids

    } else {
        assert_that(NUM_GENES == num_genes)
        assert_that(all(PC_GENE_IDS == pc_gene_ids))
        assert_that(all(LNCRNA_GENE_IDS == lncRNA_gene_ids))
        assert_that(all(gene_expression_i$gene_id == gene_expression$gene_id))
        gene_expression <- gene_expression %>% bind_cols(gene_expression_i %>% select(2))
    }
}


####################################################
writeLines("####\n 40 - Convert gene expression to 2d matrix")
####################################################
gene_expression_matrix_raw <- gene_expression %>% select(-gene_id) %>% as.matrix() %>% t()
assert_that(all(colnames(gene_expression_matrix_raw) == gene_expression_files$sample_id))
colnames(gene_expression_matrix_raw) <- gene_expression$gene_id

gene_expression_matrix_raw <- gene_expression_matrix_raw[order(rownames(gene_expression_matrix_raw)),]
gene_expression_matrix_raw <- gene_expression_matrix_raw[,order(colnames(gene_expression_matrix_raw))]


####################################################
writeLines("####\n 50 - Correct data")
####################################################
clinical_data <- read_rds(opt$CLINICAL_DATA)

logging("Processing normal samples")
gene_expression_matrix_normal <- correct_data(data_matrix = gene_expression_matrix_raw,
                                              clinical_data = clinical_data,
                                              data_set = 2,
                                              result_analysis = T,
                                              temp_file_prefix = str_c(opt$CANCER, "_gene_expression_normal_"),
                                              temp_dir = opt$TEMP_DIR,
                                              peer_hidden_samples_ratio = opt$PEER_HIDDENS_SAMPLES_RATIO)

logging("Processing tumor samples")
gene_expression_matrix_tumor <- correct_data(data_matrix = gene_expression_matrix_raw,
                                             clinical_data = clinical_data,
                                             data_set = 1,
                                             result_analysis = T,
                                             temp_file_prefix = str_c(opt$CANCER, "_gene_expression_tumor_"),
                                             temp_dir = opt$TEMP_DIR,
                                             peer_hidden_samples_ratio = opt$PEER_HIDDENS_SAMPLES_RATIO)


assert_that(!is.null(gene_expression_matrix_tumor))

if (!is.null(gene_expression_matrix_normal)) {
    common_ids <- intersect(colnames(gene_expression_matrix_tumor),
                            colnames(gene_expression_matrix_normal))
    gene_expression_matrix_corrected <- rbind(gene_expression_matrix_tumor[,colnames(gene_expression_matrix_tumor) %in% common_ids],
                                              gene_expression_matrix_normal[,colnames(gene_expression_matrix_normal) %in% common_ids])
    assert_that(length(common_ids) != 0)
    gene_expression_matrix_corrected <- gene_expression_matrix_corrected[order(rownames(gene_expression_matrix_corrected)),]
} else {
    gene_expression_matrix_corrected <- gene_expression_matrix_tumor
}



lncRNA_gene_expression_matrix <- gene_expression_matrix_corrected[,colnames(gene_expression_matrix_corrected) %in% LNCRNA_GENE_IDS]
pc_gene_expression_matrix <- gene_expression_matrix_corrected[,colnames(gene_expression_matrix_corrected) %in% PC_GENE_IDS]
TF_target <- read_rds(opt$TF_TARGET_MATRIX)
TF_gene_expression <- gene_expression_matrix_corrected[,colnames(gene_expression_matrix_corrected) %in% rownames(TF_target)]


write_rds(pc_gene_expression_matrix, opt$PC_GENE_EXPRESSION_CORRECTED, compress = "gz")
write_rds(TF_gene_expression, opt$TF_GENE_EXPRESSION_CORRECTED, compress = "gz")
write_rds(lncRNA_gene_expression_matrix, opt$LNCRNA_GENE_EXPRESSION_CORRECTED, compress = "gz")
