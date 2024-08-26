#LIB_LOCATION <- NULL
LIB_LOCATION <- "/hpf/tools/centos6/R/3.2.3/lib64/R/library"
library(readr, lib.loc = LIB_LOCATION)
library(dplyr, lib.loc = LIB_LOCATION)
library(stringr, lib.loc = LIB_LOCATION)
require(optparse, lib.loc = LIB_LOCATION) # make_option


option_list = list(
    make_option(c("-A", "--cancer"), type = "character", default=NULL,
                help="Cancer set abbreviation name", metavar = "character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

source("/hpf/largeprojects/agoldenb/mingjie/TCGA-Assembler.1.0.3/Module_B.r")
#source("/Volumes/USB128/TCGA-Assembler.1.0.3/Module_B.r")

#CANCER_TYPE <- "OV"
CANCER_TYPE <- opt[[1]]

#setwd(str_c("/Volumes/USB128/", CANCER_TYPE, "_raw"))

if (file.exists(str_c(CANCER_TYPE, "_methyl_27.rda"))) {
    # Load READ methylation27 and methylation450 data.
    load(str_c(CANCER_TYPE, "_methyl_27.rda"))
    Methylation27Data <- list(Des = Des, Data = Data);
    load(str_c(CANCER_TYPE, "_methyl_450.rda"))
    Methylation450Data <- list(Des = Des, Data = Data);
    
    # Merge READ methylation27 and methylation450 data.
    Methylation27_450_Merged = MergeMethylationData(input1 = Methylation27Data, 
                                                    input2 = Methylation450Data,
                                                    outputFileName = str_c(CANCER_TYPE, ".methyl_27_450_merge"),
                                                    outputFileFolder = "./")
}
