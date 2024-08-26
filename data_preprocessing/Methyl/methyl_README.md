# Methylation data pre-processing
## 10 - correct the spreadsheet HumanMethylation450_15017482_v1-2.csv due to outdated HMM-islands
## 15 - Perform CpG probes pruning
## 30 - Pre-process methylation text-based data from fire-browse and generate 2d methylation matrix
## 35 - Correct methylation raw matrix
## 36 - Collect corrected methylation matrix for normal sample and tumor sample
### prerequisite: 35
## 40 - Generate methylation target matrix for genes in one chromosome
## 50 - Collect methylation target matrix for genes/exons in one chromosome
### prerequisite: 40
## test - compare 27k and 450k probe
