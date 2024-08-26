
library(reshape)
library(ggplot2)

### FUNCTIONS FOR THE SCRIPT

get.top.coefs = function(in_list,summary_df){
  out_list = list()
  
  for(gene in na.omit(names(in_list))){
    cat(gene)
    gene_list = in_list[[gene]]
    
    coef_names = na.omit(Reduce(intersect,lapply(gene_list,names)))
    my_df = data.frame(matrix(nrow=length(coef_names),ncol=length(names(gene_list))))
    rownames(my_df) = coef_names
    colnames(my_df) = names(gene_list)
    for(cancer in names(gene_list)){
      my_df[na.omit(names(gene_list[[cancer]])),cancer] <- gene_list[[cancer]][na.omit(names(gene_list[[cancer]]))]
      my_df[,cancer] = as.numeric(my_df[,cancer])
    }
    if(nrow(my_df)>0){
      my_df[abs(my_df) > 5] = NA
      my_means = apply(my_df,1,function(x){mean(na.omit(x))})
      my_means[my_means > 5] = NA
      max_mean_coef_names = names(my_means)[order(abs(na.omit(my_means)),decreasing = TRUE)][1:10]
      
      out_list[[gene]] = lapply(gene_list,function(x)na.omit(x[max_mean_coef_names]))
    }
    
  }
  
  return(out_list)
}


### Getting tumour colors
tum.cols <- gg_color_hue(21)
cans <- "BLCA, BRCA, CESC, ESCA, HNSC, KIRC, KIRP, LGG, LIHC, LUAD, LUSC, PAAD, PCPG, PRAD, READ-COAD, SARC, SKCM, STAD, TGCT, THCA, UCEC"
tums <- unlist(strsplit(x = cans,split = ",\ "))
length(tums)
names(tum.cols) <- tums

reg.vec <- c("miRNA","TF","lncRNA","methyl","CNV","SNP")

### READ IN SUMMARY FILES


summary.dir = "C:/Users/larun/Desktop/Modes/Data/Info for each gene model each tissue/"
summaries <- list()
for(tum in tums){
  if(paste0(tum,"_0_data_summary.rds") %in% list.files(summary.dir)){
    summaries[[tum]] <- readRDS(paste0(summary.dir,tum,"_0_data_summary.rds"))
  }
}

num.reg.sum <- lapply(summaries,function(x){x$NUM_REGULATORS})
str(num.reg.sum)
overall.summary <- Reduce(rbind,num.reg.sum)
head(overall.summary)

overall.summary <- overall.summary[!duplicated(overall.summary$target_id),]
dim(overall.summary)
head(overall.summary)

### 
###   COMPUTE MEANS FOR miRNA, TF, lncRNA
###

  ## STARTING WITH miRNA
reg = 1
reg.name = reg.vec[reg]

cat("Regulator: ",reg.name,"\n")

coef.dir = "C:/Users/larun/Desktop/Modes/Data/model-coefs/coef-summaries/"
coefs <- list()
for(tum in tums){
  if(paste0(tum,"_",reg.name,"_means.rds") %in% list.files(coef.dir)){
    coefs[[tum]] <- readRDS(paste0(coef.dir,tum,"_",reg.name,"_means.rds"))
  }
}
length(coefs)

# str(coefs)

  ## re-arrange list
names(coefs)
head(names(coefs$BLCA))

gene_mirna_list = list()
# names(gene_mirna_list) <- overall.summary$gene_name

for(gene in overall.summary$target_id){
  gene_mirna_list[[gene]] = lapply(coefs,function(x)x[[gene]][["wt"]])
}
# names(gene_mirna_list) = overall.summary$targ

head(names(gene_mirna_list[[1]]))

mirna_top_coefs = get.top.coefs(in_list = gene_mirna_list,summary_df = overall.summary)
length(mirna_top_coefs)

saveRDS(object = mirna_top_coefs,file = "C:/Users/larun/Desktop/Modes/Data/top10_mirna_coefs.rds")


## Next, TFs: 

reg = 2
reg.name = reg.vec[reg]

cat("Regulator: ",reg.name,"\n")

coef.dir = "C:/Users/larun/Desktop/Modes/Data/model-coefs/coef-summaries/"
coefs <- list()
for(tum in tums){
  if(paste0(tum,"_",reg.name,"_means.rds") %in% list.files(coef.dir)){
    coefs[[tum]] <- readRDS(paste0(coef.dir,tum,"_",reg.name,"_means.rds"))
  }
}

length(coefs)
names(coefs)

gene_tf_list = list()

for(gene in overall.summary$target_id){
  gene_tf_list[[gene]] = lapply(coefs,function(x)x[[gene]][["wt"]])
}

tf_top_coefs = get.top.coefs(in_list = gene_tf_list,summary_df = overall.summary)
saveRDS(object = tf_top_coefs,file = "C:/Users/larun/Desktop/Modes/Data/top10_tf_coefs.rds")


## Next, lncRNA: 

reg = 3
reg.name = reg.vec[reg]

cat("Regulator: ",reg.name,"\n")

coef.dir = "C:/Users/larun/Desktop/Modes/Data/model-coefs/coef-summaries/"
coefs <- list()
for(tum in tums){
  if(paste0(tum,"_",reg.name,"_means.rds") %in% list.files(coef.dir)){
    coefs[[tum]] <- readRDS(paste0(coef.dir,tum,"_",reg.name,"_means.rds"))
  }
}

length(names(coefs))

gene_lncrna_list = list()

for(gene in overall.summary$target_id){
  gene_lncrna_list[[gene]] = lapply(coefs,function(x)x[[gene]][["wt"]])
}

lncrna_top_coefs = get.top.coefs(in_list = gene_lncrna_list,summary_df = overall.summary)
saveRDS(object = lncrna_top_coefs,file = "C:/Users/larun/Desktop/Modes/Data/top10_lncrna_coefs.rds")
lncrna_top_coefs$ENSG00000004455.12
