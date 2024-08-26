#!/bin/bash

WORKING_DIR=/hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files

echo WORKING_DIR: $WORKING_DIR
echo DATA_DIR: $DATA_DIR
echo INPUT: $INPUT

#DATA_DIR=$DATA_DIR
#INPUT=$INPUT

module load plink
module load python/2.7.9
module load R 

cd $DATA_DIR

hapmap_dir=/hpf/largeprojects/agoldenb/lauren/Reference-data-1000g-PCA/

#Rscript --vanilla /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/revising-tcga-ids.R -i $INPUT
Rscript --vanilla $WORKING_DIR/subsetting-to-healthy-samples.R -f $INPUT

#cp $INPUT.bim $INPUT-revised-ids.bim
#cp $INPUT.bed $INPUT-revised-ids.bed

plink --bfile $INPUT --keep $INPUT-healthy-inds.txt --make-bed --out $INPUT-2

plink --bfile $INPUT-2 --maf 0.05 --indep-pairwise 1500 100 0.2 --out pruned 

plink --bfile $INPUT-2 --extract pruned.prune.in --make-bed --out $INPUT-pruned

plink --bfile $INPUT-pruned --genome --out $INPUT-relatedness

## Identify individuals who are first cousins or closer:
awk '{if($10 > 0.45) print $1,$2}' $INPUT-relatedness.genome > $INPUT-related-inds-to-remove.txt

plink --bfile $INPUT-2 --remove $INPUT-related-inds-to-remove.txt --make-bed --out $INPUT-ur

plink --bfile $INPUT-ur --bmerge $hapmap_dir/1000genomes-pca-ref-indep-allpops.bed $hapmap_dir/1000genomes-pca-ref-ur-noambig-allpops.bim $hapmap_dir/1000genomes-pca-ref-ur-noambig-allpops.fam --recode --out $INPUT-hap-ceu-tsi-merged

plink --bfile $INPUT-ur --exclude $INPUT-hap-ceu-tsi-merged.missnp --make-bed --out $INPUT-missnp-rmvd

plink --bfile $INPUT-missnp-rmvd --bmerge $hapmap_dir/1000genomes-pca-ref-ur-noambig-cau.bed $hapmap_dir/1000genomes-pca-ref-ur-noambig-cau.bim $hapmap_dir/1000genomes-pca-ref-ur-noambig-cau.fam --allow-no-sex --make-bed --out $INPUT-hap-ceu-tsi-merged

plink --bfile $INPUT-hap-ceu-tsi-merged --missing --out $INPUT-hap-ceu-tsi-merged-missing

awk '{if($5 < 0.05) print}' $INPUT-hap-ceu-tsi-merged-missing.lmiss > no-lc-snp-list.txt

plink --bfile $INPUT-hap-ceu-tsi-merged --extract no-lc-snp-list.txt --make-bed --out $INPUT-hap-ceu-tsi-merged-nolc

Rscript --vanilla /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/case-ctrl-array-fam-file-creation.R -i $INPUT-hap-ceu-tsi-merged-nolc

cp $INPUT-hap-ceu-tsi-merged-nolc.bim $INPUT-hap-ceu-tsi-merged-nolc-arraycasectrl.bim
cp $INPUT-hap-ceu-tsi-merged-nolc.bed $INPUT-hap-ceu-tsi-merged-nolc-arraycasectrl.bed

plink --bfile $INPUT-hap-ceu-tsi-merged-nolc-arraycasectrl --assoc --out testing-array-against-hap

awk '{if($9 < 0.01 || $5 > 0.45 || $6 > 0.45 || $5 < 0.05 || $6 < 0.05) print $2}' testing-array-against-hap.assoc > assoc-snps-to-rmv.txt
# awk '{if($5 > 0.49 || $6 > 0.49) print $2}' testing-array-against-hap.assoc > assoc-snps-to-rmv.txt
#awk '{if($9 == "NA" || $5 > 0.4) print $2}' testing-array-against-hap.assoc > assoc-snps-to-rmv.txt

plink --bfile $INPUT-hap-ceu-tsi-merged-nolc --exclude assoc-snps-to-rmv.txt --make-bed --out $INPUT-hap-ceu-tsi-merged-nolc-noassoc-snps
#plink --bfile $INPUT-hap-ceu-tsi-merged-nolc --make-bed --out $INPUT-hap-ceu-tsi-merged-nolc-noassoc-snps

plink --bfile $INPUT-hap-ceu-tsi-merged-nolc-noassoc-snps --maf 0.05 --indep-pairwise 1500 100 0.2 --out pruned-hap-ceu-tsi-merged

plink --bfile $INPUT-hap-ceu-tsi-merged-nolc-noassoc-snps --extract pruned-hap-ceu-tsi-merged.prune.in --recode --out $INPUT-hap-ceu-tsi-merged-pruned

/hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/making-smartpca-files.sh $INPUT-hap-ceu-tsi-merged-pruned $INPUT-hap-ceu-tsi-merged-pca

rm $INPUT-tsi-ceu-hap-pca-parfile
python /hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files/making-parfile.py -p $INPUT-hap-ceu-tsi-merged-pca -o $INPUT-tsi-ceu-hap-pca-parfile

module purge
qsub $WORKING_DIR/running-smartpca.sh -v PARFILE=$INPUT-tsi-ceu-hap-pca-parfile,LOGFILE=$INPUT-hap-tsi-ceu-pca1.log,DATA_DIR=$DATA_DIR



