import sys, getopt
import pandas as pd


def main(argv):
	samplefile = ''
	inputfile = ''
	outputfile = '' 
	try:
		opts, args = getopt.getopt(argv,"hc:a:o:",["callfile=","annofile=","ofile="])
	except getopt.GetoptError:
		print 'test.py -c <callfile> -a <annofile> -o <outputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'test.py -c <callfile> -a <annofile> -o <outputfile>'
			sys.exit()
		elif opt in ("-c", "--callfile"):
			callfile = arg		 
		elif opt in ("-a", "--annofile"):
			annofile = arg
		elif opt in ("-o", "--outfile"):
			outputfile = arg

	print 'Sample file is "', callfile,'"'
	print 'Input file is "', annofile,'"'
	print 'Output file is "', outputfile,'"'
	
	return callfile,annofile,outputfile
	
callfile,annofile,outputfile = main(sys.argv[1:])

anno_infile = pd.read_csv(annofile,skiprows=18)
revised_annofile = anno_infile.rename(index=str, columns={"Probe Set ID": "affy_snp","dbSNP RS ID":"rs_id","Allele A":"allele_a","Allele B":"allele_b"})


zero_add = '0\t0'
one_add = '0\t0'
two_add = '0\t0'
neg_one_add = '0\t0'

def naming_alleles(x,z_a = zero_add,o_a = one_add,t_a = two_add,n_o_a = neg_one_add):
	if x == '0': 
		return zero_add
	if x == '1': 
		return one_add
	if x == '2': 
		return two_add
	if x == '-1': 
		return neg_one_add

split = str.split

with open(callfile,'r') as birdseed_file: 
	for line in birdseed_file:
		if line[0] == "#":
			pass
		else:
			snp = split(line)[0]
			nums1 = split(line)[1:]
			nums = [s.strip('\n') for s in nums1]
			
			if snp in revised_annofile.affy_snp.values: 
				snp_slice = revised_annofile[revised_annofile.affy_snp == snp]
				
				alleles2 = [snp_slice.get_value(snp_slice.index[0],'allele_a'),snp_slice.get_value(snp_slice.index[0],'allele_b')]

				zero_add = "\t".join([alleles2[0], alleles2[0]])
				one_add = "\t".join(alleles2)
				two_add = "\t".join([alleles2[1], alleles2[1]])
				neg_one_add = "\t".join(['0','0'])
				

				all_alleles = [naming_alleles(x = x) for x in nums]    
				 
				rs_snp = snp_slice.get_value(snp_slice.index[0],'rs_id')
				chromosome = snp_slice.get_value(snp_slice.index[0],'Chromosome')
				bp = snp_slice.get_value(snp_slice.index[0],'Physical Position')
				tped_start = [chromosome, rs_snp, '0', bp]

						
				all_alleles2 = tped_start + all_alleles
				with open(outputfile,'a') as file: 
					file.write("\t".join(all_alleles2)+'\n')

