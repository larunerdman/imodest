

#####
# Generate lncRNA target matrix
#
# Note:
#   1) cis only: one the same strand of the chromsome
#   2) Use gencode.v19.long_noncoding_RNAs.gtf for lncRNA coordinate
#   3) Use Salmon quantified pc gene expression and lncRNA gene expression


library(methods)
library(optparse)
library(Matrix)
library(dplyr)
library(tidyr)
library(stringr)
library(assertthat)
library(readr)

option_list = list(
    make_option(c("--NUM_NUCLEOTIDES"),
                type = "integer",
                default = 1000000,
                help = "Exon speicifc or gene specific", metavar = "character"),
    make_option(c("--GENE_GTF"),
                type = "character",
                default = "/Volumes/USB128/ref/gencode.v19/gencode.v19.annotation.gtf",
                help = "",
                metavar = "character"),
    make_option(c("--LNCRNA_GTF"),
                type = "character",
                default = "/Volumes/USB128/ref/gencode.v19/gencode.v19.long_noncoding_RNAs.gtf",
                help = "",
                metavar = "character"),
    make_option(c("--LNCRNA_TARGET_MATRIX_GENE"),
                type = "character",
                default = "/Users/MJ/Desktop/lncRNA_preprocess/60_lncRNA_target_matrix_gene.rds",
                help = "", metavar = "character"),
    make_option(c("--LNCRNA_TARGET_MATRIX_EXON"),
                type = "character",
                default = "/Users/MJ/Desktop/lncRNA_preprocess/60_lncRNA_target_matrix_exon.rds",
                help = "", metavar = "character")
)

print(opt <- parse_args(OptionParser(option_list = option_list)))


####
writeLines("####\n 10 - preprocess gene gtf, generate gene annotation")

gene_data <- read_tsv(opt$GENE_GTF, skip = 5, col_names = F,
                      col_types = cols(
                          X1 = col_character(),
                          X2 = col_character(),
                          X3 = col_character(),
                          X4 = col_integer(),
                          X5 = col_integer(),
                          X6 = col_character(),
                          X7 = col_character(),
                          X8 = col_character(),
                          X9 = col_character()
                      )) %>%
    select(chr = 1, type = 3, coord_start = 4, coord_end = 5, strand = 7, description = 9) %>%
    filter(type == "gene") %>%
    mutate(gene_id = description %>% str_extract("gene_id \"[\\w|\\.]+") %>% str_replace("gene_id \"", ""),
           gene_type = description %>% str_extract("gene_type \"[\\w|\\.\\-]+") %>% str_replace("gene_type \"", "")) %>%
    filter(gene_type %in% c("protein_coding", "polymorphic_pseudogene", "IG_V_gene",
                            "TR_V_gene", "IG_C_gene", "IG_J_gene", "TR_J_gene", "TR_C_gene",
                            "IG_D_gene", "TR_D_gene", "lincRNA", "miRNA")) %>%
    arrange(gene_id)

assert_that(!any(is.na(gene_data$gene_id)))
assert_that(!any(duplicated(gene_data$gene_id)))

####
writeLines("####\n 20 - preprocess lncRNA data gtf, genereate lncRNA annotation")

lncRNA_data <- read_tsv(opt$LNCRNA_GTF, skip = 5, col_names = F) %>%
    select(chr = 1, type = 3, coord_start = 4, coord_end = 5, strand = 7, description = 9) %>%
    filter(type == "gene") %>%
    mutate(gene_id = description %>% str_extract("gene_id \"[\\w|\\.]+") %>% str_replace("gene_id \"", "")) %>%
    arrange(gene_id)
assert_that(!any(is.na(lncRNA_data$gene_id)))
assert_that(!any(duplicated(lncRNA_data$gene_id)))

####
writeLines("####\n 30 - Generate lncRNA target")

lncRNA_target_gene <- Matrix(F, nrow = length(lncRNA_data$gene_id) ,
                             ncol = length(gene_data$gene_id),
                             dimnames = list(lncRNA_data$gene_id, gene_data$gene_id), sparse = T)

for (i in seq_len(nrow(lncRNA_target_gene))) {
    target_ids <- which(gene_data$chr == lncRNA_data$chr[i] &
                            gene_data$strand == lncRNA_data$strand[i] &
                            gene_data$coord_start - lncRNA_data$coord_end[i] <= opt$NUM_NUCLEOTIDES &
                            lncRNA_data$coord_start[i] - gene_data$coord_end <= opt$NUM_NUCLEOTIDES &
                            pmin(lncRNA_data$coord_end[i], gene_data$coord_end) - pmax(lncRNA_data$coord_start[i], gene_data$coord_start) <= 0)
    if (length(target_ids) > 0) {
        lncRNA_target_gene[i, target_ids] <- T
    }

    if (i %% 500 == 499) {print(i)}
}

writeLines("lncRNA target gene, number lncRNA target a gene")
table(colSums(lncRNA_target_gene))

####
writeLines("####\n 40 - Write lncRNA target gene to local")
write_rds(lncRNA_target_gene, opt$LNCRNA_TARGET_MATRIX_GENE, compress = "gz")

# 57820 genes
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
# 1798 3708 5224 6351 6412 5868 5211 4782 3828 3174 2519 1961 1467 1203 1030  789
# 16   17   18   19   20   21   22   23   24   25   26   27   28   29   30   31
# 643  479  370  220  231  158   94   67   49   40   29   25   19   12   10   10
# 32   33   34   35
# 12   16    8    3
####
writeLines("####\n 50 - preprocess in gene gtf, generate exon annotation")

exon_data <- read_tsv(opt$GENE_GTF, skip=5, col_names = F,
                      col_types = cols(
                          X1 = col_character(),
                          X2 = col_character(),
                          X3 = col_character(),
                          X4 = col_integer(),
                          X5 = col_integer(),
                          X6 = col_character(),
                          X7 = col_character(),
                          X8 = col_character(),
                          X9 = col_character()
                      )) %>%
    select(chr=1, type=3, coord_start = 4, coord_end = 5, strand = 7, description = 9) %>%
    filter(type == "exon") %>%
    mutate(exon_id = str_c(chr, strand, coord_start, coord_end, sep = "_"),
           gene_id = description %>% str_extract("gene_id \"[\\w|\\.]+") %>% str_replace("gene_id \"", ""),
           gene_type = description %>% str_extract("gene_type \"[\\w|\\.\\-]+") %>% str_replace("gene_type \"", "")) %>%
    filter(gene_type %in% c("protein_coding", "polymorphic_pseudogene", "IG_V_gene",
                            "TR_V_gene", "IG_C_gene", "IG_J_gene", "TR_J_gene", "TR_C_gene",
                            "IG_D_gene", "TR_D_gene", "lincRNA", "miRNA")) %>%
    group_by(gene_id) %>% mutate(gene_coord_start = min(coord_start), gene_coord_end = max(coord_start)) %>% ungroup() %>%
    group_by(exon_id, chr, strand) %>%  summarise(gene_coord_start = min(gene_coord_start), gene_coord_end = max(gene_coord_start)) %>% ungroup() %>%
    arrange(exon_id)

assert_that(!any(is.na(exon_data$exon_id)))
assert_that(!any(duplicated(exon_data$exon_id)))


####
writeLines("####\n 60 - Generate lncRNA target")

lncRNA_target_exon <- Matrix(F, nrow = length(lncRNA_data$gene_id) ,
                             ncol = length(exon_data$exon_id),
                             dimnames = list(lncRNA_data$gene_id, exon_data$exon_id), sparse = T)

for (i in seq_len(nrow(lncRNA_target_exon))) {
    target_ids <- which(exon_data$chr == lncRNA_data$chr[i] & exon_data$strand == lncRNA_data$strand[i] &
                            exon_data$gene_coord_start - lncRNA_data$coord_end[i] <= opt$NUM_NUCLEOTIDES &
                                lncRNA_data$coord_start[i] - exon_data$gene_coord_end <= opt$NUM_NUCLEOTIDES &
                                pmin(lncRNA_data$coord_end[i], exon_data$gene_coord_end) - pmax(lncRNA_data$coord_start[i], exon_data$gene_coord_start) <= 0)
    if (length(target_ids) > 0) {
        lncRNA_target_exon[i, target_ids] <- T
    }

    if (i %% 500 == 499) {print(i)}
}

writeLines("lncRNA target gene, number lncRNA target a gene")
table(colSums(lncRNA_target_exon))

####
writeLines("####\n 70 - Write lncRNA target gene to local")
write_rds(lncRNA_target_exon, opt$LNCRNA_TARGET_MATRIX_EXON, compress = "gz")

# 562680 exons
# 0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17
# 10799 27681 45402 54826 57938 55220 52408 46892 38762 36293 28063 24407 17099 13766 10944  9704  7691  6397
# 18    19    20    21    22    23    24    25    26    27    28    29    30    31    32    33    34    35
# 5275  2681  2826  1911  1300   850   788   640   340   425   260   244    74   107   177   194   270    26
