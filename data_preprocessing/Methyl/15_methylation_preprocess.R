# library(readr)
# library(dplyr, warn.conflicts = F)


## This script perform CpG probes pruning

#########
# Method 1: intersection between methlatuon 27k porbes and 450k probes
#########

# probeInfo_450 <- read_tsv("/Users/MJ/Desktop/methylation_preprocess/jhu-usc.edu_LGG.HumanMethylation450.14.lvl-3.TCGA-TQ-A7RS-01A-12D-A33U-05.txt",
#                skip = 1,
#                col_types = cols(
#                    `Composite Element REF` = col_character(),
#                    Beta_value = col_double(),
#                    Gene_Symbol = col_character(),
#                    Chromosome = col_character(),
#                    Genomic_Coordinate = col_integer()
#                )) %>%
#     select(probe_id=1, chr=4, coord=5) %>%
#     arrange(chr, coord)
# probeInfo_450 <- probeInfo_450[complete.cases(probeInfo_450),]

# probeInfo_27 <- read_tsv("/Users/MJ/Desktop/methylation_preprocess/jhu-usc.edu_BRCA.HumanMethylation27.2.lvl-3.TCGA-B6-A0IM-01A-11D-A032-05.txt",
#                       skip = 1,
#                       col_types = cols(
#                           `Composite Element REF` = col_character(),
#                           Beta_value = col_double(),
#                           Gene_Symbol = col_character(),
#                           Chromosome = col_character(),
#                           Genomic_Coordinate = col_integer()
#                       )) %>%
#     select(probe_id=1, Chro=4, coord=5) %>%
#     arrange(Chro, coord)
# probeInfo_27 <- probeInfo_27[complete.cases(probeInfo_27),]

# probe_id_27_450_merged <- intersect(probeInfo_27$probe_id, probeInfo_450$probe_id)
# write_rds(probe_id_27_450_merged, "/Users/MJ/Desktop/methylation_preprocess/probe_id_27_450_merged.rds")
# write_rds(probeInfo_27, "/Users/MJ/Desktop/methylation_preprocess/probeInfo_27.rds")
# write_rds(probeInfo_450, "/Users/MJ/Desktop/methylation_preprocess/probeInfo_450.rds")



##############
# Method 2:
#  Collapase CpG probes in to CpG islands
##############
# library(readr)
# library(dplyr, warn.conflicts = F)
# probeInfo_450 <- read_csv("/Users/MJ/Desktop/methylation_preprocess/HumanMethylation450_15017482_v1-2.csv", skip = 7,
#                           col_types = cols(
#                               .default = col_character(),
#                               AddressA_ID = col_integer(),
#                               AddressB_ID = col_integer(),
#                               Genome_Build = col_integer(),
#                               MAPINFO = col_integer(),
#                               Coordinate_36 = col_integer(),
#                               Random_Loci = col_logical(),
#                               Methyl27_Loci = col_logical(),
#                               Enhancer = col_logical(),
#                               DHS = col_logical()
#                           )) %>%
#     select(probe_id=1, chr=12, coord=13)
#
# probeInfo_450 <- probeInfo_450[complete.cases(probeInfo_450),] %>% arrange(chr, coord)
#
# # CpG islands initiation
# islands <- cbind(start=probeInfo_450$coord, end=probeInfo_450$coord+1, size=1, width=2, group_id=1:nrow(probeInfo_450))
# rownames(islands) <- probeInfo_450$probe_id
#
# # Extension of cpgCollapase from 'minfi' package
# # Assume probeInfor are sort by chr, coord
# cpgCollapase <- function(probeInfo, probes, islands, maxGap, maxIslandWidth) {
#     probeInfo <- probeInfo %>% filter(probe_id %in% probes)
#
#     lastHead <- which(islands[,3]==1)[1]
#     lastChr <- probeInfo$chr[lastHead]
#
#     if (is.na(lastHead)) {
#         return(islands)
#     }
#
#     for (i in (lastHead:nrow(probeInfo))[-1]) {
#         if (islands[i,3] != 1) {
#             lastHead <- NA
#             lastChr <- NA
#         } else if (is.na(lastHead) || lastChr != probeInfo$chr[i]) {
#             lastHead <- i
#             lastChr <- probeInfo$chr[i]
#         } else if (islands[i, 2] - islands[i-1, 1] > maxGap) {
#             lastHead <- i
#         } else if (islands[i, 2] - islands[lastHead, 1] > maxIslandWidth) {
#             lastHead <- i
#         } else {
#             islands[lastHead, 2:3] <- c(islands[i, 2], islands[lastHead, 3] + 1)
#             islands[i, 3] <- 0
#             islands[i, 5] <- lastHead
#         }
#     }
#
#     islands[,4] <- (islands[,2]-islands[,1] + 1) * (islands[,3] != 0)
#     return(islands)
# }
#
# ts_500_1500 <- cpgCollapase(probeInfo_450, probeInfo_450$probe_id, islands, 500, 1500)
# ts_2000_10000 <- cpgCollapase(probeInfo_450, probeInfo_450$probe_id, islands, 2000, 10000)
# ts_500_1500_2000_10000 <- cpgCollapase(probeInfo_450, probeInfo_450$probe_id, ts_500_1500, 2000, 10000)
# write_rds(ts_500_1500, "/Users/MJ/Desktop/analysis/test_3/methylation_grouping_500_1500.rds")
# write_rds(ts_2000_10000, "/Users/MJ/Desktop/analysis/test_3/methylation_grouping_2000_10000.rds")
# write_rds(ts_500_1500_2000_10000, "/Users/MJ/Desktop/analysis/test_3/methylation_grouping_500_1500_2000_10000.rds")
#
#
#
# ts_2000_7500 <- cpgCollapase(probeInfo_450, probeInfo_450$probe_id, islands, 2000, 7500)
# ts_4000_15000 <- cpgCollapase(probeInfo_450, probeInfo_450$probe_id, islands, 4000, 15000)
# ts_2000_7500_4000_15000 <- cpgCollapase(probeInfo_450, probeInfo_450$probe_id, ts_2000_7500, 4000, 15000)
# write_rds(ts_2000_7500, "/Users/MJ/Desktop/analysis/test_3/methylation_grouping_2000_7500.rds")
# write_rds(ts_4000_15000, "/Users/MJ/Desktop/analysis/test_3/methylation_grouping_4000_15000.rds")
# write_rds(ts_2000_7500_4000_15000, "/Users/MJ/Desktop/analysis/test_3/methylation_grouping_2000_7500_4000_15000.rds")
#
# ts_5000_10000 <- cpgCollapase(probeInfo_450, probeInfo_450$probe_id, islands, 5000, 10000)
# ts_10000_20000 <- cpgCollapase(probeInfo_450, probeInfo_450$probe_id, islands, 10000, 20000)
# ts_5000_10000_10000_20000 <- cpgCollapase(probeInfo_450, probeInfo_450$probe_id, ts_5000_10000, 10000, 20000)
# write_rds(ts_5000_10000, "/Users/MJ/Desktop/analysis/test_3/methylation_grouping_5000_10000.rds")
# write_rds(ts_10000_20000, "/Users/MJ/Desktop/analysis/test_3/methylation_grouping_10000_20000.rds")
# write_rds(ts_5000_10000_10000_20000, "/Users/MJ/Desktop/analysis/test_3/methylation_grouping_5000_10000_10000_20000.rds")


##############
# Method 3:
#  Collapase CpG probes in to CpG islands predicted by Wu H, et al
#   https://www.ncbi.nlm.nih.gov/pubmed/?term=Wu%20H%5BAuthor%5D&cauthor=true&cauthor_uid=20212320
##############

# library(readr)
# library(dplyr, warn.conflicts = F)
# library(stringr)
#
# probeInfo_450 <- read_csv("/Users/MJ/Desktop/methylation_preprocess/HumanMethylation450_15017482_v1-2.csv", skip = 7,
#                           col_types = cols(
#                               .default = col_character(),
#                               AddressA_ID = col_integer(),
#                               AddressB_ID = col_integer(),
#                               Genome_Build = col_integer(),
#                               MAPINFO = col_integer(),
#                               Coordinate_36 = col_integer(),
#                               Random_Loci = col_logical(),
#                               Methyl27_Loci = col_logical(),
#                               Enhancer = col_logical(),
#                               DHS = col_logical()
#                           )) %>%
#     select(probe_id=1, chr=12, coord=13)
#
# probeInfo_450 <- probeInfo_450[complete.cases(probeInfo_450),] %>% arrange(chr, coord) %>% mutate(island_id=NA)
# probeInfo_450 <- probeInfo_450 %>% mutate(id=1:nrow(probeInfo_450))
#
# # model-based-cpg-islands-hg19.txt
# predicted_cpg_islands <- read_tsv("/Users/MJ/Desktop/CGI-Hsapiens-90.txt") %>% select(1, 2, 3)
#
# predicted_cpg_islands$chr <- predicted_cpg_islands$chr %>% str_replace("chr", "")
# predicted_cpg_islands <- predicted_cpg_islands %>% arrange(chr, start) %>% mutate(island_id= 1:nrow(predicted_cpg_islands))
#
# chromosomes <- unique(intersect(predicted_cpg_islands$chr, probeInfo_450$chr))
#
# for (i in seq_len(length(chromosomes))) {
#     predicted_cpg_islands_i <- predicted_cpg_islands %>% filter(chr==chromosomes[i])
#     probeInfo_450_i <- probeInfo_450 %>% filter(chr==chromosomes[i])
#     print(i)
#     for (j in seq_len(nrow(predicted_cpg_islands_i))) {
#
#         id_j <- (probeInfo_450_i %>%
#                      filter(coord >= predicted_cpg_islands_i$start[j]) %>%
#                      filter(coord <= predicted_cpg_islands_i$end[j]))$id
#
#         if (length(id_j != 0)) {
#             probeInfo_450$island_id[id_j] <- predicted_cpg_islands_i$island_id[j]
#         }
#     }
# }
#
#


########
# Method 4: collapase based on methylation probes 450k design feature
########

library(readr)
library(dplyr, warn.conflicts = F)
library(stringr)
probeInfo_450 <- read_rds("/Users/MJ/Desktop/methylation_preprocess/probeInfo_450.rds") %>%
    select(1, 12, 13, 20, 22, 23, 24, 27, 28) %>% arrange(CHR, MAPINFO)


probeInfo_450 <- probeInfo_450[!is.na(probeInfo_450$HMM_Island) |
                                   !is.na(probeInfo_450$UCSC_RefGene_Name) |
                                   !is.na(probeInfo_450$UCSC_CpG_Islands_Name) |
                                   !is.na(probeInfo_450$Enhancer),]

#### Group by UCSC_CpG_Islands
# Simplify multiple gene group
probeInfo_450$UCSC_RefGene_Name <- probeInfo_450$UCSC_RefGene_Name %>%
    str_split(";") %>% sapply(function(x) x[[1]])

probeInfo_450$UCSC_RefGene_Group <- probeInfo_450$UCSC_RefGene_Group %>%
    str_split(";") %>% sapply(function(x) x[[1]])

# Merge 1stExon with gene body
probeInfo_450$UCSC_RefGene_Group <- probeInfo_450$UCSC_RefGene_Group %>%
    str_replace("1stExon", "Body")

# Merge TSS200 and TSS1500 into TSS
probeInfo_450$UCSC_RefGene_Group <- probeInfo_450$UCSC_RefGene_Group %>%
    str_replace("TSS200", "TSS") %>% str_replace("TSS1500", "TSS")

UCSC_gene <- probeInfo_450 %>%
    select(UCSC_RefGene_Name=4, UCSC_RefGene_Group=5) %>%
    arrange(UCSC_RefGene_Name, UCSC_RefGene_Group) %>%
    unique() %>%
    filter(!is.na(UCSC_RefGene_Name) & !is.na(UCSC_RefGene_Group))

UCSC_gene <- UCSC_gene %>% mutate(group_id=seq_len(nrow(UCSC_gene)))

UCSC_gene_group <- UCSC_gene %>% inner_join(probeInfo_450 %>% select(1, 4, 5), by = c("UCSC_RefGene_Name", "UCSC_RefGene_Group"))


#### Group by UCSC island
probeInfo_450 <- probeInfo_450 %>% filter(!IlmnID %in% UCSC_gene_group$IlmnID)

# Merge island shelf with island shore
probeInfo_450$Relation_to_UCSC_CpG_Island <- probeInfo_450$Relation_to_UCSC_CpG_Island %>%
    str_replace("_Shore", "") %>%
    str_replace("_Shelf", "")

UCSC_island <- probeInfo_450 %>%
    select(UCSC_CpG_Islands_Name=6, Relation_to_UCSC_CpG_Island=7) %>%
    arrange(UCSC_CpG_Islands_Name, Relation_to_UCSC_CpG_Island) %>%
    unique() %>%
    filter(!is.na(UCSC_CpG_Islands_Name) & !is.na(Relation_to_UCSC_CpG_Island))

UCSC_island <- UCSC_island %>% mutate(group_id=seq_len(nrow(UCSC_island)) + max(UCSC_gene_group$group_id))

UCSC_island_group <- UCSC_island %>% inner_join(probeInfo_450 %>% select(1, 6, 7), by = c("UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island"))


#### Group by HMM island
probeInfo_450 <- probeInfo_450 %>% filter(!IlmnID %in% UCSC_island_group$IlmnID)

HMM_island <- probeInfo_450 %>% select(HMM_Island=9) %>%
    arrange(HMM_Island) %>%
    unique() %>%
    filter(!is.na(HMM_Island))

HMM_island <- HMM_island %>% mutate(group_id=seq_len(nrow(HMM_island))+max(UCSC_island_group$group_id))

HMM_island_group <- HMM_island %>% inner_join(probeInfo_450 %>% select(1, 9), by = c("HMM_Island"))


#### Group by Enchancer
probeInfo_450 <- probeInfo_450 %>% filter(!IlmnID %in% HMM_island_group$IlmnID) %>% arrange(CHR, MAPINFO)

enchancer_group <- probeInfo_450 %>% select(IlmnID=1) %>% mutate(group_id=seq_len(nrow(probeInfo_450)) + max(HMM_island_group$group_id))


#### Collect grouping results
probeGroup_450 <- bind_rows(UCSC_gene_group %>% select(4, 3),
                            UCSC_island_group %>% select(4, 3),
                            HMM_island_group %>% select(3, 2),
                            enchancer_group)

#### Tesing for intercrossing grouping

probeInfo_450 <- read_rds("/Users/MJ/Desktop/methylation_preprocess/probeInfo_450.rds") %>%
    select(1, 12, 13, 20, 22, 23, 24, 27, 28) %>% arrange(CHR, MAPINFO)


probeGroup_450_ts <- probeGroup_450 %>% inner_join(probeInfo_450, by=c("IlmnID")) %>% arrange(CHR, MAPINFO)


group_ts <- probeGroup_450_ts$group_id

# error_group <- c()
#
# for (i in seq_along(group_ts)) {
#     j <- i + 1
#     if (j <= length(group_ts) && group_ts[j] != group_ts[i]) {
#         while (j <= length(group_ts) && j <= i + 100) {
#             if (group_ts[j] == group_ts[i]) {
#                 error_group <- c(error_group, i)
#             }
#             j <- j + 1
#         }
#     }
# }

new_group_id <- max(group_ts) + 1
for (i in seq_along(group_ts)) {
    group_i <- group_ts[i]
    j <- i + 1
    if (j <= length(group_ts) && group_ts[j] != group_i) {
        while (j <= length(group_ts) && j <= i + 100) {
            if (group_ts[j] == group_i) {
                group_ts[j] <- new_group_id
            }
            j <- j + 1
        }
    }
    new_group_id <- new_group_id + 1
}


probeGroup_450_ts$group_id <- group_ts

probeGroup_450_v <- probeGroup_450_ts$group_id
names(probeGroup_450_v) <- probeGroup_450_ts$IlmnID

write_rds(probeGroup_450_v, "/Users/MJ/Desktop/methylation_preprocess/probeGroup_450.rds")
