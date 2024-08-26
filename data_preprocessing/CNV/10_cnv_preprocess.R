
library(methods)
library(readr)
library(stringr)
library(dplyr, warn.conflicts = F)
library(optparse)
library(tidyr)
library(assertthat)
####
#   Generate 2D cnv matrix
#       row x col: gene cnv x patient
#
#       Note 1: chr23 in cnv data

option_list = list(
    make_option(c("--EXON_SPECIFIC"),
                type = "integer",
                default=0,
                help="Exon speicifc or gene specific", metavar = "character"),
    make_option(c("--CNV_DATA"),
                type = "character",
                default="/Users/MJ/Desktop/LIHC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt",
                help="CNV", metavar = "character"),
    make_option(c("--GENE_GTF"),
                type = "character",
                default="/Volumes/USB128/ref/gencode.v19/gencode.v19.annotation.gtf",
                help="Gene symbol to isoform conversion table",
                metavar = "character"),
    make_option(c("--TEMP_FOLDER"),
                type = "character",
                default="/Users/MJ/Desktop/cnv_preprocess/temp/",
                help="", metavar = "character"),
    make_option(c("--CNV_MATRIX"),
                type = "character",
                default="/Volumes/USB128/cnv_preprocess_output/LIHC_cnv_matrix.rds",
                help="", metavar = "character")
)

print(opt <- parse_args(OptionParser(option_list=option_list)))

################################################
writeLines("\n #### 10 - Process exon annotation")
cnv_data <- read_tsv(opt$CNV_DATA,
                     col_types = cols(
                         Sample = col_character(),
                         Chromosome = col_integer(),
                         Start = col_integer(),
                         End = col_integer(),
                         Num_Probes = col_character(),
                         Segment_Mean = col_double()
                     )) %>%
    mutate(sample_id = str_extract(Sample, "^TCGA-..-....-...") %>% str_replace_all("-", "."), chr = str_c("chr", Chromosome)) %>%
    select(sample_id, chr, segment_start = Start, segment_end = End, segment_mean = Segment_Mean) %>%
    arrange(sample_id, chr, segment_start, segment_end)


if (opt$EXON_SPECIFIC == 0) {
    ################################################
    writeLines("\n #### 20 - Process gene annotation")
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

    Chromosomes <- sort(intersect(gene_data$chr, cnv_data$chr))

    # calculate gene cnv for each chromosome
    for (i in seq_along(Chromosomes)) {
        print(Chromosomes[i])
        cnv_data_gene <- inner_join(cnv_data %>% filter(chr == Chromosomes[i]),
                                    gene_data %>% filter(chr == Chromosomes[i]), by = c("chr")) %>% select(-chr) %>%
            filter(segment_end >= coord_start & coord_end >= segment_start) %>%
            group_by(sample_id, gene_id) %>%
            mutate(overlap = pmin(segment_end, coord_end) - pmax(segment_start, coord_start) + 1) %>%
            summarise(cnv = sum(overlap * segment_mean) / sum(overlap)) %>%
            spread(sample_id, cnv, fill = NA)

        write_rds(cnv_data_gene, str_c(opt$TEMP_FOLDER, "cnv_data_gene_", Chromosomes[i], ".rds"), compress = "gz")
    }

    # Collect gene cnv from temp dir
    for (i in seq_along(Chromosomes)) {
        print(Chromosomes[i])
        cnv_data_gene_i <- read_rds(str_c(opt$TEMP_FOLDER, "cnv_data_gene_", Chromosomes[i], ".rds"))
        if (i == 1) {
            num_samples <- ncol(cnv_data_gene_i) - 1
            cnv_data_gene <- cnv_data_gene_i
        } else {
            assert_that(num_samples == ncol(cnv_data_gene_i) - 1)
            assert_that(all(colnames(cnv_data_gene_i) == colnames(cnv_data_gene)))
            cnv_data_gene <- bind_rows(cnv_data_gene, cnv_data_gene_i)
        }
    }

    # Convert to 2d matrix
    cnv_data_gene_2d <- cnv_data_gene %>% select(-gene_id) %>% as.matrix()
    rownames(cnv_data_gene_2d) <- cnv_data_gene$gene_id

    # Write to local
    cnv_data_gene_2d <- cnv_data_gene_2d[order(rownames(cnv_data_gene_2d)),]
    cnv_data_gene_2d <- cnv_data_gene_2d[,order(colnames(cnv_data_gene_2d))]
    cnv_data_gene_2d <- cnv_data_gene_2d[complete.cases(cnv_data_gene_2d),]
    write_rds(t(cnv_data_gene_2d), opt$CNV_MATRIX, compress = "gz")

} else {
    writeLines("\n #### 20 - Process exon annotation")
    ################################################
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
               gene_type = description %>% str_extract("gene_type \"[\\w|\\.\\-]+") %>% str_replace("gene_type \"", "")) %>%
        filter(gene_type %in% c("protein_coding", "polymorphic_pseudogene", "IG_V_gene",
                                "TR_V_gene", "IG_C_gene", "IG_J_gene", "TR_J_gene", "TR_C_gene",
                                "IG_D_gene", "TR_D_gene", "lincRNA", "miRNA")) %>%
        select(exon_id, chr, coord_start, coord_end) %>%
        distinct() %>%
        arrange(exon_id)

    Chromosomes <- sort(intersect(exon_data$chr, cnv_data$chr))

    # calculate exon cnv for each chromosome
    for (i in seq_along(Chromosomes)) {
        print(Chromosomes[i])
        cnv_data_exon <- inner_join(cnv_data %>% filter(chr == Chromosomes[i]),
                                    exon_data %>% filter(chr == Chromosomes[i]), by = c("chr")) %>% select(-chr) %>%
            filter(segment_end >= coord_start & coord_end >= segment_start) %>%
            group_by(sample_id, exon_id) %>%
            mutate(overlap = pmin(segment_end, coord_end) - pmax(segment_start, coord_start) + 1) %>%
            summarise(cnv = sum(overlap * segment_mean) / sum(overlap)) %>%
            spread(sample_id, cnv, fill = NA)

        write_rds(cnv_data_exon, str_c(opt$TEMP_FOLDER, "cnv_data_exon_", Chromosomes[i], ".rds"), compress = "gz")
    }

    # Collect exon cnv from temp dir
    for (i in seq_along(Chromosomes)) {
        print(Chromosomes[i])
        cnv_data_exon_i <- read_rds(str_c(opt$TEMP_FOLDER, "cnv_data_exon_", Chromosomes[i], ".rds"))
        if (i == 1) {
            num_samples <- ncol(cnv_data_exon_i) - 1
            cnv_data_exon <- cnv_data_exon_i
        } else {
            assert_that(num_samples == ncol(cnv_data_exon_i) - 1)
            assert_that(all(colnames(cnv_data_exon_i) == colnames(cnv_data_exon)))
            cnv_data_exon <- bind_rows(cnv_data_exon, cnv_data_exon_i)
        }
    }

    # Convert to 2d matrix
    cnv_data_exon_2d <- cnv_data_exon %>% select(-exon_id) %>% as.matrix()
    rownames(cnv_data_exon_2d) <- cnv_data_exon$exon_id

    # Write to local
    cnv_data_exon_2d <- cnv_data_exon_2d[order(rownames(cnv_data_exon_2d)),]
    cnv_data_exon_2d <- cnv_data_exon_2d[,order(colnames(cnv_data_exon_2d))]
    cnv_data_exon_2d <- cnv_data_exon_2d[complete.cases(cnv_data_exon_2d),]
    write_rds(t(cnv_data_exon_2d), opt$CNV_MATRIX, compress = "gz")
}
