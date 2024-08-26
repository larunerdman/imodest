import numpy as np
import pandas as pd
import sys, getopt
import scipy as sp
from scipy import stats

def main(argv):
	inputfile = ''
	outputfile = '' 
	try:
		opts, args = getopt.getopt(argv,"hg:o:",["gfile=","ofile="])
	except getopt.GetoptError:
		print 'recoding-gen-files.py -g <genfile> -o <outputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'test.py -s <samplefile> -g <genfile> -o <outputfile> \n samplefile must have PHENO column'
			sys.exit()
		elif opt in ("-g", "--gfile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputfile = arg

	print 'Input file is "', inputfile,'"'
	print 'Output file is "', outputfile,'"'

	return inputfile,outputfile
	
	
	
	
if __name__ == "__main__":
   inputfile,outputfile = main(sys.argv[1:]) 

f = open(inputfile,'r')
i=1
for line in f.xreadlines():
	numbers_str = line.split()
	snp = numbers_str[2:3]
	numbers_str6 = numbers_str[6:]
	numbers_float = [float(x) for x in numbers_str6]
	# print("numbers converted to float values")
	linearray = np.array(numbers_float)
	# print("np array created")	
	
	newarray_size = linearray.size/3
	newarray=np.ndarray(shape=(1,newarray_size),dtype=float)
	
	startnum = range(0,linearray.size,3)
	midnum = range(1,linearray.size,3)	
	endnum = range(2,linearray.size,3)	
	
	for k in range(newarray_size):
		newarray[0,k] = 0*linearray[startnum[k]]+1*linearray[midnum[k]]+2*linearray[endnum[k]]	
	
	newarray_nona = newarray[0,]
    
	with open(outputfile,"a") as outfile:	
		np.savetxt(outfile,newarray_nona[None],fmt="%s")
	# print("line added to file")
	# print("line",i,"read in") 
	i=i+1
f.close()
	