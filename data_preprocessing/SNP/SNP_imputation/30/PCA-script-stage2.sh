#!/bin/bash

echo WORKING_DIR: $WORKING_DIR
echo DATA_DIR: $DATA_DIR
echo INPUT: $INPUT
echo HAPMAP_KEEP: $HAPMAP_KEEP

module load plink
module load python/2.7.9


cd $DATA_DIR

plink --bfile $INPUT-ur --keep $HAPMAP_KEEP --make-bed --out $INPUT-nohapols

plink --bfile $INPUT-nohapols --maf 0.05 --indep-pairwise 1500 100 0.2 --out prune2

plink --bfile $INPUT-nohapols --extract prune2.prune.in --recode --out $INPUT-nohapols-pruned

/hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/making-smartpca-files.sh $INPUT-nohapols-pruned $INPUT-nohapols-pca

rm $INPUT-nohapols-pca-parfile
python /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/making-parfile.py -p $INPUT-nohapols-pca -o $INPUT-nohapols-pca-parfile

module purge

qsub -v PARFILE=$INPUT-nohapols-pca-parfile,LOGFILE=$INPUT-nohapols-pca1.log,DATA_DIR=$DATA_DIR \
	 -o $OE -e $OE \
	 -l vmem=10g,walltime=15:00:00 \
	 $WORKING_DIR/running-smartpca.sh

