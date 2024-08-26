#!/bin/bash

module load plink/1.90b3x
	
cancer=$1

plink --bfile $cancer --missing --out $cancer-miss1 --noweb

# Low call SNPs
awk '{if($5 > 0.1) print $2}' $cancer-miss1.lmiss > low-call-snps-to-rmv.txt

plink --bfile $cancer --exclude low-call-snps-to-rmv.txt --make-bed --out $cancer-nolc-snps --noweb
plink --bfile $cancer-nolc-snps --missing --out $cancer-miss2 --noweb

# Low call individuals
awk '{if($6 > 0.1) print $1,$2}' $cancer-miss2.imiss > low-call-inds-to-rmv.txt

plink --bfile $cancer-nolc-snps --remove low-call-inds-to-rmv.txt --make-bed --out $cancer-nolc --noweb


# Hardy-wienberg
plink --bfile $cancer-nolc --hardy --out $cancer-hardy --noweb
awk '{if($9 < 0.00001) print $2}' $cancer-hardy.hwe > snps-out-of-hwe.txt

plink --bfile $cancer-nolc --exclude snps-out-of-hwe.txt --make-bed --out $cancer-postqc --noweb

