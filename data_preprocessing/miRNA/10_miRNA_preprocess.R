
### Generaete miRNA target matrix

library(tibble)
library(Matrix)
library(stringr)
library(stringr)
library(tidyr)
library(dplyr)
library(readr)
library(optparse)

option_list = list(
    make_option(c("--MIRNA_TARGET"),
                type = "character",
                default = "/Volumes/USB128/ref/targetScan.v7.1/Summary_Counts.all_transcripts.txt",
                help = "Exon speicifc: 1, or gene specific: 0", metavar = "character"),
    make_option(c("--GENE_GTF"),
                type = "character",
                default = "/Users/MJ/Desktop/ref/gencode.v19/gencode.v19.annotation.gtf",
                help = "",
                metavar = "character"),
    make_option(c("--MIRNA_FAMILY"),
                type = "character",
                default = "/Users/MJ/Desktop/miRNA_preprocess/miR_Family_Info.txt",
                help = "",
                metavar = "character"),
    make_option(c("--MIRNA_TARGET_MATRIX_GENE"),
                type = "character",
                default = "/Users/MJ/Desktop/lncRNA_preprocess/10_miRNA_target_matrix_gene.rds",
                help = "", metavar = "character"),
    make_option(c("--MIRNA_TARGET_MATRIX_EXON"),
                type = "character",
                default = "/Users/MJ/Desktop/lncRNA_preprocess/10_miRNA_target_matrix_exon.rds",
                help = "", metavar = "character")
)

print(opt <- parse_args(OptionParser(option_list = option_list)))


# Read TargetScan miRNA target prediction
HAS_ID <- 9606
miRNA_target <- read_tsv(opt$MIRNA_TARGET, na="NULL", col_types = cols(
    `Transcript ID` = col_character(),
    `Gene Symbol` = col_character(),
    `miRNA family` = col_character(),
    `Species ID` = col_integer(),
    `Total num conserved sites` = col_integer(),
    `Number of conserved 8mer sites` = col_integer(),
    `Number of conserved 7mer-m8 sites` = col_integer(),
    `Number of conserved 7mer-1a sites` = col_integer(),
    `Total num nonconserved sites` = col_integer(),
    `Number of nonconserved 8mer sites` = col_integer(),
    `Number of nonconserved 7mer-m8 sites` = col_integer(),
    `Number of nonconserved 7mer-1a sites` = col_integer(),
    `Number of 6mer sites` = col_integer(),
    `Representative miRNA` = col_character(),
    `Total context++ score` = col_double(),
    `Cumulative weighted context++ score` = col_double(),
    `Aggregate PCT` = col_double()
)) %>% filter(`Species ID` == HAS_ID) %>% select(transcript_id=1,miRNA_family=3,score=16)
miRNA_target <- miRNA_target[complete.cases(miRNA_target),]

# Read in targetScan miRNA family info
miRNA_family <- read_tsv(opt$MIRNA_FAMILY,
                         col_types = cols(`miR family` = col_character(),
                                          `Seed+m8` = col_character(),
                                          `Species ID` = col_integer(),
                                          `MiRBase ID` = col_character(),
                                          `Mature sequence` = col_character(),
                                          `Family Conservation?` = col_integer(),
                                          `MiRBase Accession` = col_character()
                         )) %>%
    filter(`Species ID` == HAS_ID) %>%
    select(seed=2, miRNA_id=7)
miRNA_family <- unique(miRNA_family[complete.cases(miRNA_family),])

# Assign score for each miRNA family member
miRNA_target_complete <- miRNA_target %>% inner_join(miRNA_family, by = c("miRNA_family" = "seed")) %>% select(-miRNA_family)


# transcript Gtf
transcript_data <- read_tsv(opt$GENE_GTF,
                            skip=5, col_names = F) %>% select(type=3, descrption=9) %>%
    filter(type == "transcript") %>%
    mutate(gene_id = descrption %>% str_extract("gene_id \"[\\w|\\.]+") %>% str_replace("gene_id \"", ""),
           transcript_id = descrption %>% str_extract("transcript_id \"[\\w|\\.]+") %>% str_replace("transcript_id \"", "")) %>%
    select(-descrption, - type)

# Generate miRNA target to gene
miRNA_target_gene <- miRNA_target_complete %>%
    inner_join(transcript_data, by = c("transcript_id")) %>% select(-transcript_id) %>%
    group_by(miRNA_id, gene_id) %>% summarise(score = min(score)) %>% ungroup() %>%
    spread(gene_id, score, fill = 0)

miRNA_target_gene_2d <- miRNA_target_gene %>%
    select(-miRNA_id) %>% as.matrix() %>% Matrix(sparse = T)

rownames(miRNA_target_gene_2d) <- miRNA_target_gene$miRNA_id
write_rds(miRNA_target_gene_2d, opt$MIRNA_TARGET_MATRIX_GENE, compress = "gz")

# Generate miRNA target to exon
exon_data <- read_tsv(opt$GENE_GTF,
                      skip=5, col_names = F) %>%
    select(chr=1, type=3, coord_start = 4, coord_end = 5, strand = 7, descrption = 9) %>%
    filter(type == "exon") %>%
    mutate(exon_id = str_c(chr, strand, coord_start, coord_end, sep = "_"),
           transcript_id = descrption %>% str_extract("transcript_id \"[\\w|\\.]+") %>% str_replace("transcript_id \"", "")) %>%
    select(exon_id, transcript_id)

miRNA_target_exon <- miRNA_target_complete %>%
    inner_join(exon_data, by = c("transcript_id")) %>% select(-transcript_id) %>%
    group_by(miRNA_id, exon_id) %>% summarise(score = min(score)) %>% ungroup() %>%
    spread(exon_id, score, fill = 0)


miRNA_target_exon_2d <- miRNA_target_exon %>%
    select(-miRNA_id) %>% as.matrix() %>% Matrix(sparse = T)

rownames(miRNA_target_exon_2d) <- miRNA_target_exon$exon_id
write_rds(miRNA_target_exon_2d, opt$MIRNA_TARGET_MATRIX_EXON, compress = "gz")

