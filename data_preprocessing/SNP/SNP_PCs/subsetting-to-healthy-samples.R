## 
##    TESTING SCRIPT 
##

# setwd("C:/Users/Owner/Desktop/Goldenberg Lab/Modes of regulation project/Genotyping-stuff/Making-fam-files/")
# 
# fam.file <- read.table('LIHC-no-.fam',header=FALSE,as.is=TRUE)
# head(fam.file)
# 
# dim(fam.file[as.numeric(substr(x = fam.file$V2,start = 14,stop = 15)) > 9,])
# 
# healthy.inds <- fam.file[as.numeric(substr(x = fam.file$V2,start = 14,stop = 15)) > 9,c('V1','V2')]
# head(healthy.inds)
# dim(healthy.inds)

## 
## 

library("optparse")

# Read in options from command line
option_list = list(
  make_option(c("-f", "--famfile"), type = "character",
              help="upper case TCGA cancer name", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list)

opt = parse_args(opt_parser)

filename = opt[[1]] # creat cancer identifier variable

fam.file <- read.table(paste0(filename,'.fam'),header=FALSE,as.is=TRUE)

healthy.inds <- fam.file[as.numeric(substr(x = fam.file$V2,start = 14,stop = 15)) > 9,c('V1','V2')]

write.table(healthy.inds,file = paste0(filename,'-healthy-inds.txt'),quote = FALSE,row.names = FALSE,col.name = FALSE)