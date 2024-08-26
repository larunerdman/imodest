#!/bin/bash

echo WORKING_DIR: $WORKING_DIR
echo DATA_DIR: $DATA_DIR
echo INPUT: $INPUT

module load plink
module load python/2.7.9
module load R 

cd $DATA_DIR

hapmap_dir=$WORKING_DIR/hapmap-data

Rscript --vanilla $WORKING_DIR/revising-tcga-ids.R -i $INPUT

cp $INPUT.bim $INPUT-revised-ids.bim
cp $INPUT.bed $INPUT-revised-ids.bed

plink --bfile $INPUT-revised-ids --remove $INPUT-rmvl-ids.txt --make-bed --out $INPUT-2

plink --bfile $INPUT-2 --maf 0.05 --indep-pairwise 1500 100 0.2 --out pruned 

plink --bfile $INPUT-2 --extract pruned.prune.in --make-bed --out $INPUT-pruned

plink --bfile $INPUT-pruned --genome --out $INPUT-relatedness

## Identify individuals who are first cousins or closer:
awk '{if($10 > 0.12) print $1,$2}' $INPUT-relatedness > $INPUT-related-inds-to-remove.txt

plink --bfile $INPUT-2 --remove $INPUT-related-inds-to-remove.txt --make-bed --out $INPUT-ur

plink --bfile $INPUT-ur --bmerge $hapmap_dir/hapmap-hg19-ceu-tsi-ur.bed $hapmap_dir/hapmap-hg19-ceu-tsi-ur.bim $hapmap_dir/hapmap-hg19-ceu-tsi-ur.fam --recode --out $INPUT-hap-ceu-tsi-merged

plink --bfile $INPUT-ur --exclude $INPUT-hap-ceu-tsi-merged.missnp --make-bed --out $INPUT-missnp-rmvd

plink --bfile $INPUT-missnp-rmvd --bmerge $hapmap_dir/hapmap-hg19-ceu-tsi-ur.bed $hapmap_dir/hapmap-hg19-ceu-tsi-ur.bim $hapmap_dir/hapmap-hg19-ceu-tsi-ur.fam --make-bed --out $INPUT-hap-ceu-tsi-merged

plink --bfile $INPUT-hap-ceu-tsi-merged --missing --out $INPUT-hap-ceu-tsi-merged-missing

awk '{if($5 < 0.05) print}' $INPUT-hap-ceu-tsi-merged-missing.lmiss > no-lc-snp-list.txt

plink --bfile $INPUT-hap-ceu-tsi-merged --extract no-lc-snp-list.txt --make-bed --out $INPUT-hap-ceu-tsi-merged-nolc

Rscript --vanilla $WORKING_DIR/case-ctrl-array-fam-file-creation.R -i $INPUT-hap-ceu-tsi-merged-nolc

cp $INPUT-hap-ceu-tsi-merged-nolc.bim $INPUT-hap-ceu-tsi-merged-nolc-arraycasectrl.bim
cp $INPUT-hap-ceu-tsi-merged-nolc.bed $INPUT-hap-ceu-tsi-merged-nolc-arraycasectrl.bed

plink --bfile $INPUT-hap-ceu-tsi-merged-nolc-arraycasectrl --assoc --out testing-array-against-hap

awk '{if($9 < 0.001 || $5 > 0.45 || $6 > 0.45) print $2}' testing-array-against-hap.assoc > assoc-snps-to-rmv.txt

plink --bfile $INPUT-hap-ceu-tsi-merged-nolc --exclude assoc-snps-to-rmv.txt --make-bed --out $INPUT-hap-ceu-tsi-merged-nolc-noassoc-snps

plink --bfile $INPUT-hap-ceu-tsi-merged-nolc-noassoc-snps --maf 0.05 --indep-pairwise 1500 100 0.2 --out pruned-hap-ceu-tsi-merged

plink --bfile $INPUT-hap-ceu-tsi-merged-nolc-noassoc-snps --extract pruned-hap-ceu-tsi-merged.prune.in --recode --out $INPUT-hap-ceu-tsi-merged-pruned

$WORKING_DIR/making-smartpca-files.sh $INPUT-hap-ceu-tsi-merged-pruned $INPUT-hap-ceu-tsi-merged-pca

rm $INPUT-tsi-ceu-hap-pca-parfile
python $WORKING_DIR/making-parfile.py -p $INPUT-hap-ceu-tsi-merged-pca -o $INPUT-tsi-ceu-hap-pca-parfile

module purge
qsub -v PARFILE=$INPUT-tsi-ceu-hap-pca-parfile,LOGFILE=$INPUT-hap-tsi-ceu-pca1.log,DATA_DIR=$DATA_DIR \
	 -o $OE -e $OE \
	 -l vmem=10g,walltime=15:00:00 \
	 $WORKING_DIR/running-smartpca.sh



