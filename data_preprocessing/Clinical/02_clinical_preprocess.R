##### Normal sample collection

library(optparse)
library(readr)
library(stringr)
library(dplyr)

option_list <- list(
    # Model specification
    make_option(c("--CANCER"), type = "character", default = "LIHC",
                help = "Cancer set abbreviation name", metavar = "character"),

    # Output
    make_option(c("--OUTPUT"), type = "character",
                default = "~/Desktop/temp/LIHC_sample_ids.rds",
                help = "Output folder for the reference data", metavar = "character"),

    make_option(c("--RNA_EXPRESSION_DATA"), type = "character",
                default = "/Volumes/USB128/MoR/RNA_preprocess_output/LIHC_pc_gene_expression.rds",
                help = "RNA expression data in 2D matrix: target id x patient", metavar = "character"),

    # Predictor Matrix
    make_option(c("--CLINICAL_DATA_PICKED"), type = "character",
                default = "/Volumes/USB128/MoR/clinical_preprocess_output/LIHC_clinical_matrix.rds",
                help = "clinical features in 2D matrix: patient x clinical features", metavar = "character"),
    make_option(c("--MIRNA_EXPRESSION_DATA"), type = "character",
                default = "/Volumes/USB128/MoR/miRNA_preprocess_output/LIHC_miRNA_isoform_expression.rds",
                help = "miRNA expression in 2D matrix: patient x miRNA isoform", metavar = "character"),
    make_option(c("--LNCRNA_EXPRESSION_DATA"), type = "character",
                default = "/Volumes/USB128/MoR/lncRNA_preprocess_output/LIHC_lncRNA_gene_expression.rds",
                help = "lncRNA related gene expression in 2D matrix: patient x lncRNA gene", metavar = "character"),
    make_option(c("--TF_EXPRESSION_DATA"), type = "character",
                default = "/Volumes/USB128/MoR/TF_preprocess_output/LIHC_TF_gene_expression.rds",
                help = "TF related gene expression in 2D matrix: patient x TF gene", metavar = "character"),
    make_option(c("--METHYLATION_MATRIX"), type = "character",
                default = "/Volumes/USB128/MoR/methylation_preprocess_output/LIHC_methylation_matrix.rds",
                help = "Methylation 450K probe set in 2d matrix, patient x methylation probe", metavar = "character"),
    make_option(c("--CNV_MATRIX"), type = "character",
                default = "/Volumes/USB128/MoR/cnv_preprocess_output/LIHC_cnv_matrix_gene.rds",
                help = "cnv data in 2D matrix: patient x gene/exon cnv", metavar = "character"),
    make_option(c("--SNP_MATRIX_FOLDER"), type = "character",
                default = "/Volumes/USB128/MoR/SNP_preprocess_output/LIHC/",
                help = "Folder containing SNP data of all chromosomes in 2D matrix: patient x SNP", metavar = "character")

)

print(opt <- parse_args(OptionParser(option_list = option_list)))

opt <- parse_args(OptionParser(option_list=option_list))

# Read in data
# clinical_matrix <- read_rds(opt$CLINICAL_MATRIX)
sample_ids <- read_lines(opt$CLINICAL_DATA_PICKED) %>% str_split("\\t") %>% unlist() %>% str_subset("^tcga") %>% str_to_upper() %>% str_replace_all("-", ".")

miRNA_expression <- read_rds(opt$MIRNA_EXPRESSION_DATA)
methylation_matrix <- read_rds(opt$METHYLATION_MATRIX)
cnv_matrix <- read_rds(opt$CNV_MATRIX)

sample_ids <- sort(Reduce(
    intersect, list(
        sample_ids,
        rownames(cnv_matrix),
        rownames(miRNA_expression),
        rownames(methylation_matrix)
    )))

# write to local
print(length(sample_ids))
write_rds(sample_ids, opt$OUTPUT, compress = "gz")
