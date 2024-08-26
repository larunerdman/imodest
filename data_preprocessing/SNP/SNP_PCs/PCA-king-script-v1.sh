#!/bin/bash

# “i” and “input” have the required argument of the input file
# “r” and “reference” have no arguments, flags use of reference sample in input data

# set an initial value for the reference data flag
REF=0

# read the options
TEMP=`getopt -o d:e:i:r --long datadir:,exclude:,input:,reference -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
        -d|--datadir)
            case "$2" in
                *) DATA_DIR=$2 ; shift 2 ;;
            esac ;; 
		-e|--exclude)
            case "$2" in
                *) EXCLUDE=$2 ; shift 2 ;;
            esac ;;
		-i|--input)
            case "$2" in
                *) INPUT=$2 ; shift 2 ;;
            esac ;;
        -r|--reference) REF=1 ; shift ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

echo $EXCLUDE

WORKING_DIR=/hpf/largeprojects/grp1250/TCGA_ICGC/data/geno-calling-files
hapmap_dir=/hpf/largeprojects/agoldenb/lauren/Reference-data-1000g-PCA/

echo WORKING_DIR: $WORKING_DIR
echo DATA_DIR: $DATA_DIR
echo INPUT: $INPUT

if [ $REF -eq 1 ]
then 
	echo "Using 1000Genomes Caucasian reference samples"
fi

module load plink/1.90b3x
module load python/2.7.9
module load king
module load R

cd $DATA_DIR

if [ $REF -eq 1 ]
then
# Create list of healthy individuals to subset to. 
	Rscript --vanilla $WORKING_DIR/subsetting-to-healthy-samples.R -f $INPUT

	#cp $INPUT.bim $INPUT-revised-ids.bim
	#cp $INPUT.bed $INPUT-revised-ids.bed

# Subset data to healthy individuals only 
	plink --bfile $INPUT --keep $INPUT-healthy-inds.txt --make-bed --out $INPUT-2

# Prune data to check for relatedness
	plink --bfile $INPUT-2 --maf 0.05 --indep-pairwise 1500 100 0.2 --out pruned 

# Extract pruned SNPs
	plink --bfile $INPUT-2 --extract pruned.prune.in --make-bed --out $INPUT-pruned

# 
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

	plink --bfile $INPUT-hap-ceu-tsi-merged-nolc-noassoc-snps --extract pruned-hap-ceu-tsi-merged.prune.in --make-bed --out $INPUT-hap-ceu-tsi-merged-pruned

	king -b $INPUT-hap-ceu-tsi-merged-pruned.bed --kinship --prefix $INPUT-hap-ceu-tsi-merged-pruned
else
	if [ -z "$EXCLUDE" ]
	then
		plink --bfile $INPUT --maf 0.05 --indep-pairwise 1500 100 0.2 --out prune2

		plink --bfile $INPUT --extract prune2.prune.in --make-bed --out $INPUT-pruned

		king -b $INPUT-pruned.bed --kinship --prefix $INPUT-pruned
	else
		plink --bfile $INPUT --remove $EXCLUDE --make-bed --out $INPUT-ols-rmvd 
		
		plink --bfile $INPUT-ols-rmvd --maf 0.05 --indep-pairwise 1500 100 0.2 --out prune2

		plink --bfile $INPUT-ols-rmvd --extract prune2.prune.in --make-bed --out $INPUT-pruned

		king -b $INPUT-pruned.bed --kinship --prefix $INPUT-pruned
	fi
fi
