
'''
Convert the output from sinvict into vcf file format, and filter out all mutations that were also found in the normal bam file.
'''

import sys
import numpy as np
import os
import math
from Bio import SeqIO

if len(sys.argv) <= 1:
	print("\nConvert the output from sinvict into vcf file format.\n")
	print("Usage: python sinvict2somaticVCF.py [tumor_calls_level1.sinvict] [tumor_calls_level2.sinvict] [tumor_calls_level3.sinvict] [tumor_calls_level4.sinvict] [tumor_calls_level5.sinvict] [tumor_calls_level6.sinvict] [normal_calls_level1.sinvict] [out_somatic.vcf] [reference.fa]\n")
	#print("Optional arguments, which are off by default:")
	#print("\t--max-adj-pvalue <float>")
	sys.exit(1)
if "-h" in sys.argv[1]:
	print("\nConvert the output from sinvict into vcf file format.\n")
	print("Usage: python sinvict2somaticVCF.py [tumor_calls_level1.sinvict] [tumor_calls_level2.sinvict] [tumor_calls_level3.sinvict] [tumor_calls_level4.sinvict] [tumor_calls_level5.sinvict] [tumor_calls_level6.sinvict] [normal_calls_level1.sinvict] [out_somatic.vcf] [reference.fa]\n")
	#print("Optional arguments, which are off by default:")
	#print("\t--max-adj-pvalue <float>")
	sys.exit(1)

# functions:
# this function gets the ref base of a genomic position
def getRefBase(chromosome, pos, records):
	startPosBed=pos-1	# 0-based
	endPosBed=pos	 	# 1-based
	for record in records:
		if record.id == chromosome:
			seqArr=str(record.seq[startPosBed:endPosBed])
	seqArr=convertToUpperCase(seqArr)
	#print(seqArr)
	return seqArr

def getNeighbors(chromosome, pos, records, num_neighbors):
	startPosBed=pos-(num_neighbors+1)   # 0-based
	endPosBed=pos+num_neighbors         # 1-based
	seqArr=[]
	for record in records:
		if record.id == chromosome:
			seqArr=list(record.seq[startPosBed:endPosBed])
	#print(seqArr)
	return seqArr

def convertToUpperCase(base):
	upperCaseBase=base
	if base == "t":
		upperCaseBase="T"
	elif base == "c":
		upperCaseBase="C"
	elif base == "a":
		upperCaseBase="A"
	elif base == "g":
		upperCaseBase="G"
	elif base == "n":
		upperCaseBase="N"
	return upperCaseBase

tumor_sinvict_list=[]
for i in range(1,7):  # argument 1 through 6
	tumor_sinvict_list.append(sys.argv[i])

normal_sinvict = sys.argv[7]
outVCF = sys.argv[8]
referenceFasta = sys.argv[9]
#if len(sys.argv)>9:
#	for i in range(10,len(sys.argv)):
#	if sys.argv[i]=="--max-adj-pvalue":
#		maxPval=float(sys.argv[i+1])

infileTU_list = []
for i in range(0,6): # 6 times: 0-5
	infileTU_list.append(open(tumor_sinvict_list[i],'r'))

# open the reference fasta file to be able to check the ref bases at a position (important to write out indels)
handle = open(referenceFasta, "rU")
records = list(SeqIO.parse(handle, "fasta"))
handle.close()

infileNO = open(normal_sinvict,'r')
if os.path.exists(outVCF):
	print("Outfile already exists: %s")%(outVCF)
	sys.exit(1)
outfile = open(outVCF,'w')





# print header lines already
outfile.write("##fileformat=VCFv4.0\n")
outfile.write("##source=SiNVICT\n")
outfile.write("##reference=%s\n" % referenceFasta)
outfile.write("##INFO=<ID=TV,Number=1,Type=Integer,Description=\"Number of reads supporting variant in tumor\">\n")
outfile.write("##INFO=<ID=CT,Number=1,Type=Integer,Description=\"Coverage in tumor\">\n")
outfile.write("##INFO=<ID=VF,Number=.,Type=Float,Description=\"Tumor variant allele frequency of the SNV\">\n")
outfile.write("##INFO=<ID=CL,Number=.,Type=Float,Description=\"Confidence level of variant (one of 1,2,3,4,5,6; the bigger the better)\">\n")
outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")


# create a list with all germline mutations (all mutation found in the matched normal of lowest confidence level, i.e. including all vars)
delim="__"
germlineVars = []
cntNO=0
for line in infileNO:
	lineSplit = line.strip().split("\t")
	chr=lineSplit[0]
	pos=lineSplit[1]
	ref=lineSplit[3]
	var=lineSplit[5]
	varID=chr +delim+ pos +delim+ ref +delim+ var
	cntNO+=1
	germlineVars.append(varID)

print("Found %d variants in the matched normal file." % cntNO)

cntTUtotal=0
cntSomatic=0
alreadyAdded=[]
# loop over lists with variants
for i in range(5,-1,-1): # 5 4 3 2 1 0 # we want to start with the best ones, because they are subsets of the worse lists, and we want the correct qual score for each variant
	currentFile=tumor_sinvict_list[i]
	# check that file is not empty:
	if os.stat(currentFile).st_size == 0:
		print("File is empty: %s." % currentFile)
		continue
	print("Currently working on: %s." % currentFile)
	# otherwise read in the mutations and write them into the somatic vcf if they are not in the matched normal, and if they were not already written out in a higher confidence level
	for line in infileTU_list[i]:
		if i==0: # if it is the lowest confidence file, we count all mutations to know how many there were in total
			cntTUtotal+=1
		lineSplit = line.strip().split("\t")
		chr=lineSplit[0]
		posRaw=lineSplit[1]
		refRaw=lineSplit[3]
		varRaw=lineSplit[5]
		varID=chr +delim+ posRaw +delim+ refRaw +delim+ varRaw
		
		# only go ahead if not in matched normal:
		if varID in germlineVars:
			continue
		# also only go ahead if not already written out:
		if varID in alreadyAdded:
			continue
		
		TV=int(lineSplit[6])
		CT=int(lineSplit[4])
		VAFperc=float(lineSplit[7])
		VF=(float(VAFperc)/100)
		CL=int(i)+1
		
		# now check whether it is an indel or not
		if varRaw.startswith("-"): # deletion
			# get the position and ref allele one base before
			pos=int(posRaw)-1
			refAlleleBefore=getRefBase(chr, pos, records)
			ref=refAlleleBefore + varRaw.replace("-","")
			# make sure that all alleles are in capital letters:
			refIndBases=list(ref)
			refIndBasesCapital=[]
			for thisBase in refIndBases:
				capitalI=convertToUpperCase(thisBase)
				refIndBasesCapital.append(capitalI)
			# join it back together to a string
			ref="".join(refIndBasesCapital)	
			var=refAlleleBefore
		elif varRaw.startswith("+"): # insertion
			pos=int(posRaw)
			ref=refRaw
			var=refRaw+varRaw.replace("+","")
			# make sure that all alleles are in capital letters:
			varIndBases=list(var)
			varIndBasesCapital=[]
			for thisBase in varIndBases:
				capitalI=convertToUpperCase(thisBase)
				varIndBasesCapital.append(capitalI)
			# join it back together to a string
			var="".join(varIndBasesCapital)
		else:
			pos=int(posRaw)
			ref=refRaw
			var=varRaw
		# write it out
		alreadyAdded.append(varID)
		outfile.write(
			'%s\t%d\t.\t%s\t%s\t%d\t.\tTV=%d;CT=%d;VF=%.5f;CL=%d\t.\n' % (chr,pos,ref,var,CL,TV,CT,VF,CL))
		cntSomatic+=1
	infileTU_list[i].close()
outfile.close()
print("Found in total %d variants in the tumor file." % cntTUtotal)
percentageFiltered=float((float(cntSomatic)/float(cntTUtotal))*100.0)
print("Wrote %d variants to the somatic VCF file. Filtered out %d variants as germline. Percentage of somatic mutations from total: %.2f%%.\n" %(cntSomatic,cntTUtotal-cntSomatic,percentageFiltered))
