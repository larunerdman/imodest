
#####
# Generate lncRNA target matrix
#
# Note:
#   1) cis only: one the same strand of the chromsome
#   2) Use gencode.v19.long_noncoding_RNAs.gtf for lncRNA coordinate
#   3) Use TANRIC for lncRNA epression

library(methods)
library(readr)
library(tidyr)
library(dplyr, warn.conflicts = F)
library(stringr)
library(tibble)
library(assertthat)
library(Matrix)
library(optparse)

option_list = list(
    make_option(c("--NUM_NUCLEOTIDES"),
                type = "integer",
                default = "1000000",
                help = "Exon speicifc or gene specific", metavar = "character"),
    make_option(c("--EXON_SPECIFIC"),
                type = "integer",
                default = "1",
                help = "Exon speicifc: 1, or gene specific: 0", metavar = "character"),
    make_option(c("--UCSC_EXON_INFO"),
                type = "character",
                default = "/Users/MJ/Desktop/RNA_preprocess/Ucsc_Exon_Info.rds",
                help = "Gene symbol to isoform conversion table",
                metavar = "character"),
    make_option(c("--LNCRNA_DATA"),
                type = "character",
                default = "/Users/MJ/Desktop/lncRNA_preprocess/gencode.v19.long_noncoding_RNAs.gtf",
                help = "",
                metavar = "character"),
    make_option(c("--c"),
                type = "integer",
                default = 1,
                help = "Specify the chromosome for exon-specific model",
                metavar = "character"),
    make_option(c("--TEMP_DIR"),
                type = "character",
                default = "/hpf/largeprojects/agoldenb/mingjie/lncRNA_preprocess/temp/",
                help = "", metavar = "character"),
    make_option(c("--OUTPUT_DIR"),
                type = "character",
                default = "/hpf/largeprojects/agoldenb/mingjie/lncRNA_preprocess/",
                help = "", metavar = "character")
)


opt <- parse_args(OptionParser(option_list = option_list))
print(opt)

Num_nucleotides <- opt$NUM_NUCLEOTIDES

#### 10 ) Read in gencode lncRNA reference data
lncRNA_data <- read_tsv(opt$LNCRNA_DATA,
                        skip = 5, col_names = F) %>% select(chr=1, type = 3, coord_start = 4, coord_end = 5, strand = 7, id = 9) %>%
    filter(type == "gene") %>%
    mutate(gene_id = id %>% str_extract("gene_id \"ENSGR?\\d+") %>% str_extract("ENSGR?\\d+")) %>%
    select(-id) %>% arrange(gene_id)

assert_that(!any(is.na(lncRNA_data$gene_id)))
assert_that(!any(duplicated(lncRNA_data$gene_id)))


#### 20 - Create lncRNA target matrix for isoform specific model
Ucsc_Gene_Info <- read_rds(opt$UCSC_EXON_INFO) %>%
    select(-exon_id, -isoform_id) %>%
    group_by(ucsc_gene_symbol, chr, strand) %>%
    summarise(coord_start = min(coord_start), coord_end = max(coord_end)) %>% ungroup() %>%
    group_by(ucsc_gene_symbol) %>%
    filter(n() == 1) %>% ungroup() %>%
    arrange(ucsc_gene_symbol) %>%
    mutate(id = seq_len(n()))

if (opt$EXON_SPECIFIC == 0) {
    LncRNA_Target <- Matrix(F, nrow = length(lncRNA_data$gene_id) ,
                            ncol = length(Ucsc_Gene_Info$ucsc_gene_symbol),
                            dimnames = list(lncRNA_data$gene_id, Ucsc_Gene_Info$ucsc_gene_symbol), sparse = T)


    for (i in seq_along(lncRNA_data$gene_id)) {
        lncRNA_data_i <- lncRNA_data[i,]
        target_i <- Ucsc_Gene_Info %>% filter((chr == lncRNA_data_i$chr) & (strand == lncRNA_data_i$strand) &
                                                  (coord_start - lncRNA_data_i$coord_end <= Num_nucleotides) &
                                                  (lncRNA_data_i$coord_start  - coord_end <= Num_nucleotides) &
                                                  (pmin(lncRNA_data_i$coord_end, coord_end) - pmax(lncRNA_data_i$coord_start, coord_start) <= 0))
        if (length(target_i$id) > 0) {
            LncRNA_Target[i, target_i$id] <- T
        }

        if (i %% 500 == 499) {print(i)}
    }

    assert_that(!is.unsorted(rownames(LncRNA_Target)))
    assert_that(!is.unsorted(colnames(LncRNA_Target)))

    #hist(colSums(LncRNA_Target))

    # write to local
    write_rds(LncRNA_Target, str_c(opt$OUTPUT_DIR, "40_lncRNA_target_score_gene.rds"), compress = "gz")
}


######
## 30 - Generate exon-specific lncRNA target score by proximity
if (opt$EXON_SPECIFIC == 1) {
    Ucsc_Exon_Info <- read_rds(opt$UCSC_EXON_INFO) %>%
        select(-isoform_id) %>%
        inner_join(Ucsc_Gene_Info %>%
                       rename(gene_coord_start = coord_start, gene_coord_end = coord_end),
                   by = c("chr", "strand", "ucsc_gene_symbol")) %>%
        distinct() %>%
        group_by(exon_id, chr, strand) %>%
        summarise(gene_coord_start = min(gene_coord_start), gene_coord_end = max(gene_coord_end)) %>% ungroup() %>%
        arrange(exon_id)

    assert_that(!any(duplicated(Ucsc_Exon_Info$exon_id)))

    #Exon_ids <- unique(Ucsc_Exon_Info$exon_id)
    #Ucsc_Exon_Info$group_id <- Ucsc_Exon_Info %>% group_indices(exon_id)

    Chromosomes <- intersect(Ucsc_Exon_Info$chr, lncRNA_data$chr)

    lncRNA_data_c <- lncRNA_data %>% filter(chr == Chromosomes[opt$c])
    Ucsc_Exon_Info_c <- Ucsc_Exon_Info %>% filter(chr == Chromosomes[opt$c]) %>% mutate(id = seq_len(n()))
    print(Chromosomes[opt$c])

    LncRNA_Target_Exon <- Matrix(F, nrow = length(lncRNA_data_c$gene_id) ,
                                 ncol = length(Ucsc_Exon_Info_c$exon_id),
                                 dimnames = list(lncRNA_data_c$gene_id, Ucsc_Exon_Info_c$exon_id), sparse = T)


    for (i in seq_along(lncRNA_data_c$gene_id)) {
        lncRNA_data_i <- lncRNA_data_c[i,]
        target_i <- Ucsc_Exon_Info_c %>% filter((gene_coord_start - lncRNA_data_i$coord_end <= Num_nucleotides) &
                                                    (lncRNA_data_i$coord_start  - gene_coord_start <= Num_nucleotides) &
                                                    (pmin(lncRNA_data_i$coord_end, gene_coord_start) - pmax(lncRNA_data_i$coord_start, gene_coord_start) <= 0))

        if (length(target_i$id) > 0) {
            LncRNA_Target_Exon[i, target_i$id] <- T
        }

        if (i %% 500 == 499) {print(i)}

    }

    assert_that(!is.unsorted(rownames(LncRNA_Target_Exon)))
    assert_that(!is.unsorted(colnames(LncRNA_Target_Exon)))
    #hist(colSums(LncRNA_Target_Exon))

    # write to local
    write_rds(LncRNA_Target_Exon, str_c(opt$TEMP_DIR, "40_lncRNA_target_score_exon_",  Chromosomes[opt$c], ".rds"), compress = "gz")

}



# lncRNA_data_i <- lncRNA_data[i,]
# target_i <- Ucsc_Gene_Info %>% filter((chr == lncRNA_data_i$chr) & (strand == lncRNA_data_i$strand) &
#                                           (coord_start - lncRNA_data_i$coord_end <= Num_nucleotides) &
#                                           (lncRNA_data_i$coord_start  - coord_end <= Num_nucleotides) &
#                                           (abs(lncRNA_data_i$coord_start - coord_start) >= 100 |abs(lncRNA_data_i$coord_end - coord_end) >= 100)
# )
#
#
# target_j <- Ucsc_Gene_Info %>% filter((chr == lncRNA_data_i$chr) & (strand == lncRNA_data_i$strand) &
#                                           (coord_start - lncRNA_data_i$coord_end <= Num_nucleotides) &
#                                           (lncRNA_data_i$coord_start  - coord_end <= Num_nucleotides) &
#                                           (pmin(lncRNA_data_i$coord_end, coord_end) - pmax(lncRNA_data_i$coord_start, coord_start) <= 0))
# )
#
# target_k <- Ucsc_Gene_Info %>% filter((chr == lncRNA_data_i$chr) & (strand == lncRNA_data_i$strand) &
#                                           (coord_start - lncRNA_data_i$coord_end <= Num_nucleotides) &
#                                           (lncRNA_data_i$coord_start  - coord_end <= Num_nucleotides)
#
# )
