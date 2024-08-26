library("optparse")

# Read in options from command line
option_list = list(
  make_option(c("-i", "--infile"), type = "character",
              help="infile identifier", metavar = "character")
)

opt_parser = OptionParser(option_list=option_list)

opt = parse_args(opt_parser)

infile = opt[[1]] # creat infile identifier variable

outlier <- read.table(infile,header=FALSE,as.is=TRUE)
outlier$id <- sub(":.*","",outlier$V3)

outliers.to.rmv.df <- outlier[,c('id','id')]

write.table(outliers.to.rmv.df,file = paste0('prepped-for-plink-',infile),quote=FALSE,row.names=FALSE,col.names=FALSE)
