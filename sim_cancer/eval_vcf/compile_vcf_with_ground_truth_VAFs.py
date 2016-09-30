import sys
import math
import random
import string
import os

if len(sys.argv) < 3:
	print("\nFrom all the mutations in the subclones and their abundancies in the tumor, create the entire list with all mutations and their ground truth variant allele frequencies.")
	print("Usage: python compile_vcf_with_ground_truth_VAFs.py [subclone_var_list_dir] [abundance_file.txt] [out_dir]")
	print("This program needs to be run twice:")
	print("1) create bed file with all variants, which is to be intersected with all subclone's bedfiles: all_variants_TEMP.bed")
	print("2) continue with reading in for each variant, in which subclone's regions it would fall into. Compute VAF.")
	sys.exit(1)

var_dir = sys.argv[1]
abundance_file = sys.argv[2]
out_dir = sys.argv[3]


## functions
# this function gets a variant and a bed file, and checks whether the variant is in the bed region
def check_var_region(chr, pos, ref, var, bed_file):
	infile=open(bed_file, 'r')
	isIn=False
	for line in infile:
		lineSplit=line.strip().split("\t")
		thisChr=lineSplit[0]
		if thisChr == chr:
			thisStart=int(lineSplit[1])+1  # start position in bed file is 0-based
			thisEnd=int(lineSplit[2])
			var_length=max([ len(ref), len(var) ])  # variant may be indel, then we want it to be included fully in the region
			var_end=int(pos)+int(var_length)-1
			if thisStart <= pos and var_end <= thisEnd: # then the variant is within the bed regions!
				isIn=True
				break
	return isIn


random.seed(44)

# define outfile and make sure it does not exist yet
outName = out_dir + "/all_tumor_eval.vcf"
if os.path.exists(outName):
	print("Outfile already exists: %s")%(outName)
	sys.exit(1)

# read in the abundancies of the subclones
infile = open(abundance_file, 'r')
abundance_reads=[]
sum=0
firstline = infile.readline() # just the first line with the number of reads assigned to each subclone copy 0, 1
line_split = firstline.strip().split()  # splits on whitespace
numClones = len(line_split)
for i in range(0, numClones):
	abundance_reads.append(float(line_split[i]))
	sum += int(line_split[i])
abundance_perc=[]
for i in range(0, numClones):
	abundance_perc.append(float(float(abundance_reads[i])/float(sum)))
	print("Abundance of subclone %d: %.5f" % (i+7,abundance_perc[i]))

all_vars = {} # create a dictionary that will contain all variants and the list in which subclones it is (for the numerator of the VAF)
delim='__'
for i in range(7,15): # loop over all subclones
	for j in range(0,2):  # loop over all copies (diploid)
		thisClone=str(i) + "_" + str(j) + "_TU"
		for currSubcloneCopy in [thisClone, thisClone + '_gain']:
			currSubcloneCopyVcf = var_dir +  "/" + currSubcloneCopy + ".vcf"
			if os.path.exists(currSubcloneCopyVcf):
				print('Found vcf file %s. Read in the variants and add them to \'all_vars\'.\n' % currSubcloneCopyVcf)
				infile = open(currSubcloneCopyVcf, 'r') # read in the variants from the current file
				for line in infile:
					line_split = line.strip().split("\t")
					this_chr=line_split[0]
					if this_chr == "#CHROM":  # header line
						continue 
					this_pos=line_split[1]
					this_ref=line_split[3]
					this_var=line_split[4]
					varkey = this_chr + delim + this_pos + delim + this_ref + delim + this_var
					if varkey in all_vars:
						all_vars[varkey].append(currSubcloneCopy)
					else:
						all_vars[varkey]=[currSubcloneCopy]
# print one entry in all_vars:
print(all_vars[varkey])

# write out all the variants into a file to intersect with each subclone bed file
outNameVars = out_dir + "/all_variants_TEMP.bed"
if os.path.exists(outNameVars):
	print("Outfile already exists: %s")%(outNameVars)
else:
	outfile = open(outNameVars, 'w')
	for var_key in all_vars:
		var_split = var_key.split(delim)
		this_chr = var_split[0]
		this_pos = int(var_split[1])-1 # bed format has a 0-based start position
		this_ref = var_split[2]
		this_var = var_split[3]
		var_length = max([ len(this_ref), len(this_var) ])  # variant may be indel, then we want it to be included fully in the region
		var_end=this_pos+int(var_length)
		outfile.write("%s\t%d\t%d\t%s\t%s\n" % (this_chr, this_pos, var_end, this_ref, this_var))
	sys.exit(1)


# Now check for each variant whether the region is in a gain/loss
all_vars_regions = {} # create another dictionary that will contain for each variant, the bed files where it would fall into (for the denominator of the VAF)
for i in range(7,15): # loop over all subclones
	for j in range(0,2):  # loop over all copies (diploid)
		thisClone=str(i) + "_" + str(j) + "_TU"
		for currSubcloneCopy in [thisClone, thisClone + '_gain']:
			currSubcloneCopyBed = var_dir +  "/" + currSubcloneCopy + "_all_variants_regions.bed"
			if os.path.exists(currSubcloneCopyBed):
				print('Found bed file %s. Read in the variants to know which ground truth variants of \'all_vars\' are in these regions.' % currSubcloneCopyBed)
				infile = open(currSubcloneCopyBed, 'r')
				for line in infile:
					line_split = line.strip().split("\t")
					this_chr=line_split[0]
					this_pos=int(line_split[1])+1  # reading from 0-based bed format
					this_ref=line_split[3]
					this_var=line_split[4]
					varkey = this_chr + delim + str(this_pos) + delim + this_ref + delim + this_var
					if varkey in all_vars_regions:
						all_vars_regions[varkey].append(currSubcloneCopy)
					else:
						all_vars_regions[varkey]=[currSubcloneCopy]

# print one enrty from all_vars_regions:
print(all_vars_regions[varkey])


# # check whether all variants from all_vars are also in all_vars_regions:
# for var_key in all_vars:
# 	if not var_key in all_vars_regions:
# 		print("var_key in all_vars but not in all_vars_regions: %s" % var_key)
# 
# sys.exit(1)
# Note: Not all vars in all_vars are also in all_vars_regions, because the normal mutations were not restricted to the bed regions yet


# Now it is possible to compute the VAF for each variant and print it to outfile
outfile = open(outName,'w')
outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")
for var_key in all_vars_regions:
	denom=0.0
	numer=0.0
	thisVAF=-1.0
	for i in range(0,8): # loop over all subclones(i+7)
		numCopiesRegion=0
		numCopiesVar=0
		for thisClone in [str(i+7) + "_0_TU", str(i+7) + "_1_TU", str(i+7) + "_0_TU_gain", str(i+7) + "_1_TU_gain"]:
			#print(thisClone)
			#print(all_vars_regions[var_key])
			if thisClone in all_vars_regions[var_key]:
				numCopiesRegion+=1
			if thisClone in all_vars[var_key]:
				numCopiesVar+=1
		denom+=float(numCopiesRegion)*float(abundance_perc[i]/float(2.0))
		numer+=float(numCopiesVar)*float(abundance_perc[i]/float(2.0))
	thisVAF=float(float(numer)/float(denom))
	var_split = var_key.split(delim)
	this_chr = var_split[0]
	this_pos = var_split[1]
	this_ref = var_split[2]
	this_var = var_split[3]
	outfile.write("%s\t%s\t.\t%s\t%s\t%.3f\t.\tSRF=1,SRR=1,SAF=1,SAR=1,FREQ=%f\t.\n" % (this_chr, this_pos, this_ref, this_var, thisVAF, thisVAF))

