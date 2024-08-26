#!/bin/bash

CANCERS=(BLCA CESC ESCA KIRC LIHC LUSC PAAD PRAD SARC STAD THCA UCEC BRCA COAD GBM HNSC KIRP LGG LUAD OV PCPG READ SKCM TGCT THYM)

for i in "${CANCERS[@]}" ; do 

	CANCER=$i
	cp $CANCER.fam /hpf/largeprojects/agoldenb/mingjie/SNP_preprocess_input/$CANCER/genotype-calls-july-2016/

done
