import os
import pandas as pd
import sys, getopt
import subprocess as sp

def main(argv):
	cancer = ''
	try:
		opts, args = getopt.getopt(argv,"hc:",["cancer="])
	except getopt.GetoptError:
		print ' -c cancer'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print '-c cancer name (capitalized)'
			sys.exit()
		elif opt in ("-c", "--cancer"):
			cancer = arg		 

	print 'Cancer is "', cancer,'"'
	
	return cancer
	
cancer = main(sys.argv[1:])


#os.system('plink --tped ' + cancer + '-no-.tped --tfam ' + cancer + '.fam --make-bed --out ' + cancer + '-1 --noweb')

#os.system('plink --bfile ' + cancer + '-1 --keep ' + cancer + '.sample_id.tsv --make-bed --out ' + cancer + ' --noweb')

os.system('plink --tped ' + cancer + '-no-.tped --tfam ' + cancer + '-no-.fam --make-bed --out ' + cancer + ' --noweb')

os.system('plink --bfile %s --missing --out miss1 --noweb' % cancer)

# Low call SNPs
lmiss_file=pd.read_table('miss1.lmiss',header=0,sep=r'\s*')
lowcall_snps_df = lmiss_file.SNP[lmiss_file.F_MISS > 0.1]
n_lc_snps=len(lowcall_snps_df.tolist())

with open('qc.out','a') as outfile:
	outfile.write(str(n_lc_snps) + ' SNPs removed for having a call rate below 0.1' + '\n\n')

with open('low-call-snps-to-rmv.txt','a') as lcsnps:
	lcsnps.write('\n'.join(lowcall_snps_df.tolist()))	 
print 'Low call SNPs identified and written to \'low-call-snps-to-rmv.txt\''

os.system('plink --bfile %s --exclude low-call-snps-to-rmv.txt --make-bed --out nolc-snps --noweb' % cancer)
os.system('plink --bfile nolc-snps --missing --out miss2 --nowe b')

# Low call individuals
imiss_file=pd.read_table('miss2.imiss',header=0,sep=r'\s*')
lowcall_inds_df = imiss_file[imiss_file.F_MISS > 0.05]

n_lc_inds = lowcall_inds_df.FID.count()

with open('qc.out','a') as outfile:
	outfile.write(str(n_lc_inds) + ' individuals removed for having a call rate below 0.05' + '\n\n')

lowcall_inds_df.to_csv('low-call-inds-to-rmv.txt',sep='\t')

os.system('plink --bfile nolc-snps --remove low-call-inds-to-rmv.txt --make-bed --out nolc --noweb')
os.system('plink --bfile nolc-snps --hardy --out hardy --noweb')

# Hardy-wienberg
hardy_file = pd.read_table('hardy.hwe',header=0,sep=r'\s*')
snps_out_of_hwe = hardy_file.SNP[hardy_file.P < 1e-05]

n_snps_out_of_hwe = len(snps_out_of_hwe.tolist())

with open('qc.out','a') as outfile:
	outfile.write(str(n_snps_out_of_hwe) + ' SNPs removed for having a chi-sq p-value < 1e-05' + '\n\n')

with open('snps-out-of-hwe.txt','a') as hweoutfile:
	hweoutfile.write('\n'.join(snps_out_of_hwe.tolist()))

os.system('plink --bfile nolc --exclude snps-out-of-hwe.txt --make-bed --out nolc-nohwe --noweb')
os.system('plink --bfile nolc-nohwe --check-sex --out sex-check --noweb')

# Sex check
sex_file = pd.read_table('sex-check.sexcheck',header=0,sep=r'\s*')
prob_sex_inds = sex_file[sex_file.STATUS == 'PROBLEM']

n_prob_sex_inds = prob_sex_inds.FID.count()

with open('qc.out','a') as outfile:
	outfile.write(str(n_prob_sex_inds) + ' individuals removed for having ambiguous sex status'+ '\n\n')

prob_sex_inds.to_csv('ambig-sex-inds.txt',sep='\t')

os.system('plink --bfile nolc-nohwe --remove ambig-sex-inds.txt --make-bed --out ' + cancer + '-postqc')


with open(cancer+'-postqc.log') as logfile:
	for line in logfile:
		if 'Total genotyping rate' in line:
			special_line1=line			
		elif 'pass filters and QC' in line:
			special_line2=line
		else:
			continue

with open('qc.out','a') as outfile:
	if 'special_line1' in locals() or 'special_line1' in globals():
		outfile.write(''.join(special_line1) + '\n')
	if 'special_line2' in locals() or 'special_line2' in globals():
		outfile.write(''.join(special_line2) + '\n')
	
os.system('plink --bfile ' + cancer + '-postqc --maf 0.05 --indep-pairwise 1500 100 0.2 --out pruned --noweb')
os.system('plink --bfile ' + cancer + '-postqc --extract pruned.prune.in --recodeAD --out ' + cancer + ' --noweb')

rawfile = pd.read_csv(cancer+'.raw',header=0,sep=r'\s*')
rawfile2 = rawfile.iloc[:,range(0,6) + range(6,len(rawfile.columns),2)]

rawfile2.to_csv(cancer+'-pruned-snps.txt',sep='\t',index=False)
