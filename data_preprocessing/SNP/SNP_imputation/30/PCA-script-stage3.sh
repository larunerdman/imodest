#!/bin/bash

echo WORKING_DIR: $WORKING_DIR
echo DATA_DIR: $DATA_DIR
echo INPUT: $INPUT

cd $DATA_DIR

module load R 
module load plink
module load python/2.7.9

Rscript --vanilla /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/creating-pca-ol-file.R -i $INPUT-nohapols-pca.outlier

plink --bfile $INPUT-nohapols --remove prepped-for-plink-$INPUT-nohapols-pca.outlier --make-bed --out ${INPUT}2

plink --bfile ${INPUT}2 --maf 0.05 --indep-pairwise 1500 100 0.2 --out pruned3

plink --bfile ${INPUT}2 --extract pruned3.prune.in --recode --out ${INPUT}2-pruned

/hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/making-smartpca-files.sh ${INPUT}2-pruned ${INPUT}2-pca

rm ${INPUT}2-pca-parfile
python /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/making-parfile.py -p ${INPUT}2-pca -o ${INPUT}2-pca-parfile

module purge

qsub -v PARFILE=${INPUT}2-pca-parfile,LOGFILE=${INPUT}2-pca.log,DATA_DIR=$DATA_DIR \
	 -o $OE -e $OE \
	 -l vmem=10g,walltime=15:00:00 \
	 $WORKING_DIR/running-smartpca.sh