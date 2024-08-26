library("reshape")

###
###   FUNCTIONS 
### 

get.r2.addition.ZERO <- function(df,row.num,ylim=c(-1,1)){
  vars.removed <- matrix(nrow=0,ncol=4)
  colnames(vars.removed) <- c("var.added","base.model","r2.growth","mod.r2")
  
  df[df < 0] <- 0
  
  for(k in rev(unique(nchar(colnames(df))))[1:(length(rev(unique(nchar(colnames(df)))))-1)]){
    mods <- colnames(df)[which(nchar(colnames(df)) == k)]
    mods.min.1 <- colnames(df)[which(nchar(colnames(df)) == (k-1))]
    
    for(i in 1:6){
      vars.removed <- rbind(vars.removed,
                            cbind(rep(i,length(mods.min.1[grep(pattern = i,x = mods.min.1,invert = TRUE)])),
                                  mods.min.1[grep(pattern = i,x = mods.min.1,invert = TRUE)],
                                  unlist(df[row.num,mods[grep(pattern = i,x = mods)]] - 
                                           df[row.num,mods.min.1[grep(pattern = i,x = mods.min.1, invert = TRUE)]]),
                                  unlist(df[row.num,mods[grep(pattern = i,x = mods)]] - 
                                           rep(0,length(df[row.num,mods.min.1[grep(pattern = i,x = mods.min.1, invert = TRUE)]])))))
    }
  }
  
  # dim(vars.removed)
  # head(vars.removed)
  
  vars.removed <- as.data.frame(vars.removed)
  vars.removed$r2.growth <- as.numeric(as.character(vars.removed$r2.growth))
  vars.removed$mod.r2 <- as.numeric(as.character(vars.removed$mod.r2))
  vars.removed$var.added.f <- factor(vars.removed$var.added,levels=1:6,labels=c("miRNA","TF","lncRNA","methyl","CNV","SNP"))
  vars.removed$base.model <- as.character(vars.removed$base.model)
  vars.removed$gene <- rownames(df)[row.num]
  # str(vars.removed)
  # head(vars.removed)
  
  return(vars.removed) 
  
  ## DARKER = VARIABLE ADDED TO MODEL WITH FEWER PREDICTORS
  # graph.list <- split(x = vars.removed$r2.growth,f = vars.removed$var.added.f)
  # alpha.list <- split(x = 1-nchar(x = vars.removed$base.model)/10,f = vars.removed$var.added.f)
  
  # par(mfrow=c(1,1))
  # stripchart(split(x = vars.removed$r2.growth,f = vars.removed$var.added.f),method = "jitter",vertical = T,pch=19,col = add.alpha(col = "darkorchid",alpha = 0.7))
  
  ## STRIPCHART
  # stripchart(split(x = vars.removed$r2.growth,f = vars.removed$var.added.f),method = "jitter",vertical = T,pch=19,col = 'white',ylim=ylim,
  #            xlab = "",ylab = "",cex=2)
  # for(k in 1:6){
  #   points(x = rep(k,length(graph.list[[k]])) + rnorm(n = length(graph.list[[k]]),0,0.1),y = graph.list[[k]],col = add.alpha(col = "blue",alpha.list[[k]]*2),pch=19,cex=2)
  #   points(x = k, y = df[row.num,paste0("glm_",k)],pch = 18,col = 'red',cex=2)
  # }
  
  
}

create_var_added_files <- function(r2.df.list = ol.genes.r2.dfs,gene.name.ref.list = gene.tum.details.df,
                                   cancer="UCEC",cluster=FALSE){
  if(length(grep(pattern = cancer,tumour.file.list)) == 1){
    # i = 1
    var.added.list <- list()
    var.added.df <- data.frame(matrix(nrow=nrow(t(r2.df.list[[cancer]])),ncol=nrow(get.r2.addition.ZERO(df = t(r2.df.list[[cancer]]),row.num = 1))))
    for(i in 1:nrow(t(r2.df.list[[cancer]]))){
      var.added.df[i,] <- get.r2.addition.ZERO(df = t(r2.df.list[[cancer]]),row.num = i)$r2.growth
      var.added.list[[i]] <- get.r2.addition.ZERO(df = t(r2.df.list[[cancer]]),row.num = i)
    }
    ## SHOULD BE THE SAME ACROSS CANCERS, HAD ISSUE WITH THCA SO JUST LEAVING IT WITH HNSC
    rownames(var.added.df) <- gene.name.ref.list[["HNSC"]][[1]]$gene_name[match(row.names(t(r2.df.list[[cancer]])),gene.name.ref.list[["HNSC"]][[1]]$target_id)]
    names(var.added.list)  <- gene.name.ref.list[["HNSC"]][[1]]$gene_name[match(row.names(t(r2.df.list[[cancer]])),gene.name.ref.list[["HNSC"]][[1]]$target_id)]
    
    saveRDS(var.added.list,paste0("C:/Users/larun/Desktop/Modes/Data/Clustering-output/",cancer,"-var-added-list.rds"))
    saveRDS(var.added.df,paste0("C:/Users/larun/Desktop/Modes/Data/Clustering-output/",cancer,"-var-added-df.rds"))
    
    if(cluster){
      cat("Clustering started... \n")
      dist.out <- dist(var.added.df)
      hclust.out <- hclust(dist.out)
      #
      clusts <- sapply(2:100,function(i){cutree(hclust.out,i)})
      sil.means <- sapply(2:100,function(i){test <- cutree(tree = hclust.out,k = i) ; sil.out <- silhouette(x = test,dist = dist.out) ; return(mean(sil.out[,3]))})
      sil.sds <- sapply(2:100,function(i){test <- cutree(tree = hclust.out,k = i) ; sil.out <- silhouette(x = test,dist = dist.out) ; return(sd(sil.out[,3]))})
      
      saveRDS(dist.out,paste0("C:/Users/larun/Desktop/Modes/Data/Clustering-output/",cancer,"-dist-mat.rds"))
      saveRDS(hclust.out,paste0("C:/Users/larun/Desktop/Modes/Data/Clustering-output/",cancer,"-hclust.rds"))
      saveRDS(sil.means,paste0("C:/Users/larun/Desktop/Modes/Data/Clustering-output/",cancer,"-sil-means.rds"))
      saveRDS(sil.sds,paste0("C:/Users/larun/Desktop/Modes/Data/Clustering-output/",cancer,"-sil-sds.rds"))
      saveRDS(clusts,paste0("C:/Users/larun/Desktop/Modes/Data/Clustering-output/",cancer,"-clusts.rds"))
      cat("...clustering complete.\n")
    }
    
    cat(paste0("Cancer: ",cancer," complete. \n"))
  }  
}

###
###     Reading in data files and getting simple stats 
###         

cans <- "BLCA, BRCA, CESC, ESCA, HNSC, KIRC, KIRP, LGG, LIHC, LUAD, LUSC, PAAD, PCPG, PRAD, READ-COAD, SARC, SKCM, STAD, TGCT, THCA, UCEC"
tums <- unlist(strsplit(x = cans,split = ",\ "))
length(tums)
names(tum.cols) <- tums
tum.vec <- tums
cancers.of.int <- tums

#####
###   READING IN TUMOUR FILES -- R2
#####

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

## Confirm 21 PRESS R2 files read in 
length(tumour.r2.dfs)
## Confirm r2 df named for  
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

## 
##    CLUSTERING GENES -- only needs to be run once, output written to files
##

library(cluster)

cancers.of.int <- c("HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PAAD","PCPG","PRAD","UCEC","ESCA","SARC","SKCM","CESC")

gene.info.list = as.data.frame(gene.tum.details.df[[cancer]][[1]])
head(gene.info.list)

all.an.cancers <- c("HNSC","KIRC","KIRP","LIHC","LUAD","LUSC","PAAD","PCPG","PRAD","UCEC","ESCA","SARC","SKCM","CESC")
all.an.cancers <- tums

gene.names.list <- lapply(t(tum.r2.dfs),colnames)
overlapping.genes <- na.omit(Reduce(f = function(x,y)intersect(x,y),gene.names.list))
str(overlapping.genes) ## 10362
#
ol.genes.r2.dfs <- lapply(t(tum.r2.dfs),function(x){x[,overlapping.genes]})
names(ol.genes.r2.dfs) <- names(tum.r2.dfs)

## 
## LOOP TO GENERATE VARIANCE ADDED DF'S AND CLUSTERS
##

for(my.cancer in cancers.of.int){
  create_var_added_files(cancer=my.cancer)
}

