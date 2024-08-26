
## CALLING OPTIONS FROM COMMANDLINE
library(optparse)

option_list <- list(
  make_option(c("--REGULATORS"), type = "numeric", default=1,
              help = "regulators",metavar = "character"),
  make_option(c("--CANCER"), type="character", default="UCEC",
              help = "cancer", metavar = "character")
)

opt <- parse_args(OptionParser(option_list=option_list))
cancer = opt$CANCER
reg=opt$REGULATORS

  # OUTPUT FOLDER CURRENTLY FIXED BUT CAN BE ADDED TO OPTIONS
out.folder="/hpf/largeprojects/agoldenb/lauren/Modes/model_coef_summary"

# MATRIX CONVERTING GENE IDs from ENSEMBLE TO GENE SYMBOLS
transmat.genes <- (readRDS(paste0("/hpf/largeprojects/agoldenb/mingjie/temp/data_summary/",cancer,"_0_data_summary.rds"))[["NUM_REGULATORS"]][,c("gene_name","target_id")])
names(transmat.genes) <- c("new","old")

## TRANSMAT CHANGES THE NAME OF THE REGULATOR, IF REG IS MIRNA (1) OR TF (2) WANT TO CONVERT NAMES TO COMMONLY USED SYMBOL
transmat=NULL
write.rds = TRUE

if(reg == 1){
  mirna.table <- read.table("/hpf/largeprojects/agoldenb/mingjie/temp/miR_Family_info.txt",sep = "\t",header = TRUE)
  transmat <- mirna.table[,c("MiRBase.ID","MiRBase.Accession")]
  names(transmat) <- c("new","old")
}

if(reg == 2){
  transmat <- (readRDS(paste0("/hpf/largeprojects/agoldenb/mingjie/temp/data_summary/",cancer,"_0_data_summary.rds"))[["NUM_REGULATORS"]][,c("gene_name","target_id")])
  names(transmat) <- c("new","old")
  
}

## GENERATE A LIST OF ALL THE MODEL NUMBERS
filenamenumbers <- unlist(sapply(1:6,function(x){
  comb.mat <- combn(x = 1:6,m = x)
  apply(comb.mat,2,function(y){paste(y,collapse = "")})
}))

## READ IN PRESS R2 STATS
cancer.model.acc = readRDS(paste0("/hpf/largeprojects/agoldenb/mingjie/temp/predict_on_same_tissue_cohort_Sep19_cached/",cancer,"_0_tumor_r_square_1.rds"))

## SET CANCER FOLDER NAME
cancer.folder.name = paste0("/hpf/largeprojects/agoldenb/mingjie/temp/model_coeff_Sep19_cached/", cancer, "_0")

## WANT TO RUN AVERAGE COMPUTATION IN PARALLEL
require(doParallel)
require(foreach)

## GENERATE LIST OF REGULATORS FOR EACH GENE: EACH LIST ELEMENT NAMED BY GENE AND CONTAINS VECTOR OF REGULATORS 
  ## FROM THE SPECIFIED GROUP WHICH WERE INCLUDED IN THE MODEL FOR THAT GENE
reg.names <- lapply(readRDS(paste0(cancer.folder.name,"/",cancer, "_0_model_coeff_glm_",reg,".rds")),names)

## CREATE VECTOR WITH NAMES OF GENES THAT HAVE HAD THEIR EXPRESSION MODELED
ol.genes <- names(reg.names)

cl <- makeCluster(detectCores() - 2)
registerDoParallel(cl)

## CREATE LIST OF COEFFICIENTS WHERE EACH ELEMENT IS A MODEL USING THE GIVEN REGULATOR
  ## WITHIN EACH MODEL ELEMENT, THERE IS A LIST OF GENES 
  ## WITHIN EACH GENE LIST, THERE IS A VECTOR OF NAMED COEFFICIENTS FOR THAT GENE IN THAT MODEL
genes.coef.list <- foreach(i = grep(reg,filenamenumbers)) %dopar% {
  cat(paste0("Reading file: ",paste0(cancer.folder.name,"/",cancer, "_0_model_coeff_glm_",filenamenumbers[i],".rds\n")))

  coef.file = readRDS(paste0(cancer.folder.name,"/",cancer, "_0_model_coeff_glm_",filenamenumbers[i],".rds"))
  
  out.list <- sapply(ol.genes,function(x){
    coef.file[[x]][reg.names[[x]]]
  })     
}

##  MAKE SURE TO ONLY INCLUDE GENES WHICH HAVE HAD EVERY MODEL RUN 
names.list <- lapply(genes.coef.list,names)
ol.genes <- Reduce(intersect,names.list)
ol.genes <- na.omit(ol.genes[ol.genes %in% colnames(cancer.model.acc)])

## CALCULATE WEIGHTED AND UNWEIGHTED COEFFICIENT MEANS ACROSS ALL MODELS THE REGULATOR WAS IN
gene.coef.means <- sapply(ol.genes,function(x){
  out.mat <- Reduce(cbind,lapply(genes.coef.list,function(y)y[[x]]))
  colnames(out.mat) <- paste0("glm_",filenamenumbers[grep(reg,filenamenumbers)])
  if(!is.null(transmat)){
    rownames(out.mat) <- transmat$new[match(rownames(out.mat),transmat$old)]
  }
  
  gene.wts <- cancer.model.acc[,x]
  gene.wts[gene.wts < 0] <- 0
  gene.wts <- gene.wts[colnames(out.mat)]
  
  gene.nowt.coef.means <- apply(out.mat,1,function(x){mean(na.omit(x))})
  gene.wt.coef.means <- apply(out.mat,1,function(x){weighted.mean(x = x,w = gene.wts,na.rm = TRUE)})
  
  gene.coef.list <- list(gene.nowt.coef.means,gene.wt.coef.means)
  names(gene.coef.list) <- c("nowt","wt")
  
  return(gene.coef.list)
  })

## FORMAT COEFFICIENT MEANS INTO A LIST WHERE EACH ELEMENT IS A GENE
  ## WITHIN EACH GENE ELEMENT, THERE ARE TWO ELEMENTS "nowt" and "wt" 
  ## WHERE "nowt" IS THE NON-WEIGHTED MEAN COEFFICIENTS 
  ## "wt" IS THE WEIGHTED MEAN COEFFICIENTS
gene.coef.means2 <- list()
for(k in 1:(length(ol.genes))){
  my.out.list <- list(gene.coef.means[[((k)*2)-1]],
                      gene.coef.means[[((k)*2)]])
  names(my.out.list) <- c("nowt","wt")
  gene.coef.means2[[ol.genes[k]]] <- my.out.list  
}

saveRDS(gene.coef.means2,file = paste0(out.folder,"/",cancer,"_",c("miRNA","TF","lncRNA","methyl","CNV","SNP")[as.numeric(as.character(reg))],"_means.rds"))  

