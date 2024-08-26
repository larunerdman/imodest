
####
# Genreate methylation target matrix for genes in one chromosome
#   Output format: methylation probe x gene symbol / exon id

library(methods)
library(readr)
library(dplyr)
library(optparse)
library(stringr)
library(Matrix)
library(assertthat)

option_list = list(
    make_option(c("--NUM_NUCLEOTIDES"),
                type = "integer",
                default = 1000000,
                help = "Exon speicifc or gene specific", metavar = "character"),
    make_option(c("--EXON_SPECIFIC"),
                type = "integer",
                default = 0,
                help = "Exon speicifc or gene specific", metavar = "character"),
    make_option(c("--c"),
                type = "integer",
                default = 1,
                help = "Exon speicifc or gene specific", metavar = "character"),
    make_option(c("--GENE_GTF"),
                type = "character",
                default = "/Volumes/USB128/ref/gencode.v19/gencode.v19.annotation.gtf",
                help = "Gene symbol to isoform conversion table",
                metavar = "character"),
    make_option(c("--METHYLATION_PROBE_INFO_450"),
                type = "character",
                default = "/Users/MJ/Desktop/methylation_preprocess/probeInfo_450.rds",
                help = "",
                metavar = "character"),
    make_option(c("--TEMP_DIR"),
                type = "character",
                default = "/Users/MJ/Desktop/methylation_preprocess/temp/",
                help = "", metavar = "character")
)


print(opt <- parse_args(OptionParser(option_list = option_list)))


####################################################
writeLines("####\n preprocess methylation probe")

methylation_data <- read_rds(opt$METHYLATION_PROBE_INFO_450) %>%
    select(probe_id = 1, chr = 12, coord = 13) %>%
    na.omit() %>%
    mutate(chr = str_c("chr", chr)) %>%
    arrange(probe_id)

if (opt$EXON_SPECIFIC == 0) {

    ####################################################
    writeLines("####\n preprocess gene annotation")
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

    Chromosomes <- intersect(methylation_data$chr, gene_data$chr)
    methylation_data <- methylation_data %>% filter(chr %in% Chromosomes) %>% mutate(id = seq_len(n()))
    gene_data <- gene_data %>% filter(chr %in% Chromosomes) %>% mutate(id = seq_len(n()))

    methylation_target_gene <- Matrix(F, nrow = length(methylation_data$probe_id),
                                      ncol = length(gene_data$gene_id),
                                      dimnames = list(methylation_data$probe_id,
                                                      gene_data$gene_id),
                                      sparse = T)

    methylation_data <- methylation_data %>% filter(chr == Chromosomes[opt$c])
    gene_data <- gene_data %>% filter(chr == Chromosomes[opt$c])

    for (i in seq_len(nrow(gene_data))) {
        target_i <- gene_data[i,]
        methylation_ids <- methylation_data$id[target_i$coord_start - methylation_data$coord <= opt$NUM_NUCLEOTIDES &
                                                   methylation_data$coord - target_i$coord_end <= opt$NUM_NUCLEOTIDES]
        if (length(methylation_ids) > 0) {
            methylation_target_gene[methylation_ids, target_i$id] <- T
        }
        if (i %% 1000 == 999) {print(i)}
    }


    assert_that(!is.unsorted(rownames(methylation_target_gene)))
    assert_that(!is.unsorted(colnames(methylation_target_gene)))

    ####################################################
    writeLines("####\n Write to local")
    print(Chromosomes[opt$c])
    write_rds(methylation_target_gene, str_c(opt$TEMP_DIR, "40_methylation_target_matrix_gene_", Chromosomes[opt$c],".rds"), compress = "gz")

} else {

    ####################################################
    writeLines("####\n preprocess exon annotation")
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

    Chromosomes <- intersect(methylation_data$chr, exon_data$chr)
    methylation_data <- methylation_data %>% filter(chr %in% Chromosomes) %>% mutate(id = seq_len(n()))
    exon_data <- exon_data %>% filter(chr %in% Chromosomes) %>% mutate(id = seq_len(n()))

    methylation_target_exon <- Matrix(F, nrow = length(methylation_data$probe_id) ,
                                      ncol = length(exon_data$exon_id),
                                      dimnames = list(methylation_data$probe_id, exon_data$exon_id),
                                      sparse = T)

    methylation_data <- methylation_data %>% filter(chr == Chromosomes[opt$c])
    exon_data <- exon_data %>% filter(chr == Chromosomes[opt$c])

    for (i in seq_len(nrow(exon_data))) {
        target_i <- exon_data[i,]
        methylation_ids <- methylation_data$id[target_i$gene_coord_start - methylation_data$coord <= opt$NUM_NUCLEOTIDES &
                                                   methylation_data$coord - target_i$gene_coord_end <= opt$NUM_NUCLEOTIDES]
        if (length(methylation_ids) > 0) {
            methylation_target_exon[methylation_ids, target_i$id] <- T
        }
        if (i %% 1000 == 999) {print(i)}
    }


    assert_that(!is.unsorted(rownames(methylation_target_exon)))
    assert_that(!is.unsorted(colnames(methylation_target_exon)))

    ####################################################
    writeLines("####\n Write to local")
    print(Chromosomes[opt$c])
    write_rds(methylation_target_exon, str_c(opt$TEMP_DIR, "40_methylation_target_matrix_exon_", Chromosomes[opt$c],".rds"), compress = "gz")
}
