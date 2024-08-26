# This script was designed to correct the spreadsheet
#    HumanMethylation450_15017482_v1-2.csv due to outdated HMM-islands

library(readr)
library(dplyr, warn.conflicts = F)
library(stringr)
probeInfo_450 <- read_csv("/Users/MJ/Desktop/methylation_preprocess/HumanMethylation450_15017482_v1-2.csv", skip = 7,
                          col_types = cols(
                              .default = col_character(),
                              AddressA_ID = col_integer(),
                              AddressB_ID = col_integer(),
                              Genome_Build = col_integer(),
                              MAPINFO = col_integer(),
                              Coordinate_36 = col_integer(),
                              Random_Loci = col_logical(),
                              Methyl27_Loci = col_logical(),
                              Enhancer = col_logical(),
                              DHS = col_logical()
                          )) %>%
    select(-15, -16)

HMM_islands <- read_tsv("/Users/MJ/Desktop/methylation_preprocess/model-based-cpg-islands-hg19.txt") %>% select(chr=1, start=2, end=3) %>% arrange(chr)
HMM_islands$chr <- HMM_islands$chr %>% str_replace("chr", "")
HMM_islands <- HMM_islands %>% mutate(island_id=str_c(chr,":", start, "-", end))

# Remap HMM islands
chromosomes <- unique(intersect(HMM_islands$chr, probeInfo_450$CHR))

probeInfo_450 <- probeInfo_450 %>% mutate(id=seq_len(nrow(probeInfo_450)))
probeInfo_450$HMM_Island <- NA

for (i in seq_len(length(chromosomes))) {
    HMM_islands_i <- HMM_islands %>% filter(chr==chromosomes[i])
    probeInfo_450_i <- probeInfo_450 %>% filter(CHR==chromosomes[i])
    print(i)
    for (j in seq_len(nrow(HMM_islands_i))) {

        id_j <- (probeInfo_450_i %>%
                     filter(MAPINFO >= HMM_islands_i$start[j]) %>%
                     filter(MAPINFO <= HMM_islands_i$end[j]))$id

        if (length(id_j != 0)) {
            probeInfo_450$HMM_Island[id_j] <- HMM_islands_i$island_id[j]
        }
    }
}

write_rds(probeInfo_450, "/Users/MJ/Desktop/methylation_preprocess/probeInfo_450.rds")