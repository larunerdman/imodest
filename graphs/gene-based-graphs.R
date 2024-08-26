######################################################
######################################################
#
#              Gene-based graphs 
#
######################################################
######################################################

## user-defined functions
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
gg_color_hue <- function(n) {
  ## ggplot default rainbow colors from: https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette 
  # hues = seq(15, 375, length = n + 1)
  hues = seq(25, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## required libraries
library("ggplot2")
library("reshape")
library("RColorBrewer")

## setting graph colors 
tum.cols <- gg_color_hue(21)
reg.cols <- brewer.pal(n = 6,name = "Dark2")

###///////////////////////////////////////////////////##
###
###     Reading in data files and getting simple stats 
###         
###///////////////////////////////////////////////////##


cans <- "BLCA, BRCA, CESC, ESCA, HNSC, KIRC, KIRP, LGG, LIHC, LUAD, LUSC, PAAD, PCPG, PRAD, READ-COAD, SARC, SKCM, STAD, TGCT, THCA, UCEC"
tums <- unlist(strsplit(x = cans,split = ",\ "))
length(tums)
names(tum.cols) <- tums
tum.vec <- tums
cancers.of.int <- tums

  #####
  ###   READING IN TUMOUR R2 FILES
  #####

  ## SET TO YOUR GIVEN PATH TO GET R2 AND SUMMARY FILES
tumour.file.list <- list.files(path = 'C:/Users/larun/Desktop/Modes/Data//tumour pred model performance/r2')
# tumour.file.list <- tumour.file.list[unlist(sapply(tum.vec, function(cancer){grep(pattern = cancer,x = tumour.file.list)}))]
tumour.file.list <- tumour.file.list[na.omit(unlist(sapply(tums,function(tum){grep(tum,tumour.file.list)})))]
tumour.file.list <- tumour.file.list[grep("_tumor_r_square_1.rds",tumour.file.list)]

tumour.r2.dfs <- list()
gene.tum.details.df <- list()
clinic.list <- list()

k=1
for(tumour.file.in in tumour.file.list){
  
  ## Reading in PRESS R2 results
  tumour.r2.dfs[[k]] <- readRDS(paste0('C:/Users/larun/Desktop/Modes/Data//tumour pred model performance/r2/',tumour.file.in)) 
  cancer <- gsub(pattern = "_.*$",replacement = "",x = tumour.file.in)
  names(tumour.r2.dfs)[k] <- cancer
  
  
  ## Reading in cancer model details
  gene.tum.details.df[[k]] <- readRDS(paste0('C:/Users/larun/Desktop/Modes/Data/Info for each gene model each tissue/',list.files('C:/Users/larun/Desktop/Modes/Data/Info for each gene model each tissue/')[grep(pattern = cancer,x = list.files('C:/Users/larun/Desktop/Modes/Data/Info for each gene model each tissue/'))]))
  names(gene.tum.details.df)[k] <- cancer
  
  clinic.list[[k]] <-  readRDS(paste0('C:/Users/larun/Desktop/Modes/Data/clinical-model-input/',list.files('C:/Users/larun/Desktop/Modes/Data/clinical-model-input/')[grep(pattern = cancer,x = list.files('C:/Users/larun/Desktop/Modes/Data/clinical-model-input//'))]))
  names(clinic.list)[k] <- cancer
  
  k=k+1
}

# save(list = c("tumour.dfs"),file = "C:/Users/larun/Desktop/Modes/App/App prototype/Data/corr-tumour.dfs.RData")

length(tumour.r2.dfs)
names(tumour.r2.dfs)

# 
tum.r2.dfs  <- tumour.r2.dfs#[tum.vec]
tum.r2.dfs.sub  <- tumour.r2.dfs#[tum.vec]

names(gene.tum.details.df)
length(gene.tum.details.df)

### MAKING OVERALL GENE SUMMARY
num.reg.sum <- lapply(gene.tum.details.df,function(x){x$NUM_REGULATORS})
overall.summary <- Reduce(rbind,num.reg.sum)
overall.summary <- overall.summary[!duplicated(overall.summary$target_id),]


########################################################################
#####**************************************************************#####
##
##          PLOTTING INDIVIDUAL GENES 
##
#####**************************************************************#####
########################################################################

r2.added.list <- list()
tum.vec <- tums
for(cancer in tum.vec){
  r2.added.list[[cancer]] <- readRDS(paste0("C:/Users/larun/Desktop/Modes/Data/Clustering-output/",cancer,"-var-added-list.rds"))
}


## e.g. 
mygene = "TP53"

gene.list <- lapply(r2.added.list,function(x){
  if(mygene %in% names(x)){
    newdf <- x[[mygene]]
    rownames(newdf) <- 1:nrow(newdf)
    return(newdf)
  } else{
    return(NA)
  }
})

full.gene.df <- Reduce(rbind,gene.list[names(gene.list[is.na(gene.list) == FALSE])])
gene_spec_cancers = names(gene.list[is.na(gene.list) == FALSE])
full.gene.df$cancer <- rep(gene_spec_cancers,each = nrow(gene.list[is.na(gene.list) == FALSE][[1]]))

# full.gene.df.mlt <- melt(full.gene.df)
# full.gene.df.mlt <- full.gene.df.mlt[full.gene.df.mlt$variable == "r2.growth",]

ggplot(full.gene.df,aes(x = var.added.f, y = r2.growth,fill = cancer)) + 
  geom_boxplot() + xlab("Regulator") + ylab("R2 added") + ylim(c(-1,1)) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) + 
  scale_fill_manual(values = tum.cols)


### 
###     V2 -- PLOTTING CONDITIONAL GAIN IN PRESS R2
###

cond.var = c(1,2)
excl.var = c(6)

# sapply output: vectors of same length
# reduce: recurses over intersect, 
# I only want models will these cond vars in them. How much do these add to the model, when 
# all the base models have both 1 and 2 in them 

graph_df = full.gene.df[Reduce(intersect,data.frame(sapply(cond.var,function(reg)grep(reg,full.gene.df$base.model)))),]

# if excluded var is in any of the base models, excl this models
# double check empty lists bugs
graph_df = graph_df[grep(pattern = paste0(excl.var,collapse = "|"),x = graph_df$base.model,invert = TRUE),]

ggplot(graph_df,aes(x = var.added.f, y = r2.growth,fill = cancer)) + 
  geom_boxplot() + xlab("Regulator") + ylab("R2 added") + ylim(c(-1,1)) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) + 
  scale_fill_manual(values = tum.cols)
