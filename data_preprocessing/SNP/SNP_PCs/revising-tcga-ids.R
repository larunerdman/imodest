library("optparse")

# Read in options from command line
option_list = list(
  make_option(c("-i", "--infile"), type = "character",
              help="infile identifier", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list)

opt = parse_args(opt_parser)

infile = opt[[1]] # creat infile identifier variable

fam.file <- read.table(paste0(infile,'.fam'),header=FALSE,as.is=TRUE)
fam.file$revised.iid <- substr(x = fam.file$V1,start = 1,stop = 16)
fam.file$revised.iid <- gsub("-",".",fam.file$revised.iid)

write.table(fam.file[,c('revised.iid','revised.iid','V3','V4','V5','V6')],file = paste0(infile,'-revised-ids.fam'),quote=FALSE,row.names=FALSE,col.names=FALSE)

rmvl.ids <- fam.file$revised.iid[substr(fam.file$revised.iid,1,4) != 'TCGA']

write.table(data.frame(cbind(rmvl.ids,rmvl.ids)),file = paste0(infile,'-rmvl-ids.txt'),quote=FALSE,row.names=FALSE,col.names=FALSE)
