import sys
import math
import random
import string


## assumption: none of the cnv regions overlap each other! (should be the case anyway)

cnv_file = sys.argv[1]
out_dir = sys.argv[2]

random.seed(44)


infile=open(cnv_file, 'r')

myList_names=["7891011121314_0_TU","7_0_TU_8_0_TU_9_0_TU_10_0_TU","11_0_TU_12_0_TU_13_0_TU_14_0_TU","7_0_TU_8_0_TU","9_0_TU_10_0_TU","11_0_TU_12_0_TU","13_0_TU_14_0_TU","7_0_TU","8_0_TU","9_0_TU","10_0_TU","11_0_TU","12_0_TU","13_0_TU","14_0_TU"]
myList_cnvs =[[]             ,[]     ,[]        ,[]  ,[]   ,[]    ,[]    ,[] ,[] ,[] ,[]  ,[]  ,[]  ,[]  ,[]]

for line in infile:
	line_split = line.strip().split("\t")
	chr=line_split[0]
	posStart=line_split[1]
	posEnd=line_split[2]
	copyNumber=line_split[3] # one of 0(hom del), 1(het del), 3(het ampl), 4 (hom ampl)

	cnkey = chr + '_' + posStart + '_' + posEnd + '_' + copyNumber

	treeLevel=random.random() # random number between 0 and 1
	if treeLevel > 0.9:  # 7891011121314
		myList_cnvs[0].append(cnkey)
	elif treeLevel > 0.6: # 78910 11121314
		nodeInTree=random.random()
		if nodeInTree < 0.5:
			myList_cnvs[1].append(cnkey)
		else:
			myList_cnvs[2].append(cnkey)
	elif treeLevel > 0.3: # 78 910 1112 1314
		nodeInTree=random.random()
		if nodeInTree < 0.25:
			myList_cnvs[3].append(cnkey)
		elif nodeInTree < 0.5:
			myList_cnvs[4].append(cnkey)
		elif nodeInTree < 0.75:
			myList_cnvs[5].append(cnkey)
		else:
			myList_cnvs[6].append(cnkey)
	else: # all leaves of tree
		nodeInTree=random.random()
		if nodeInTree < 0.125:
			myList_cnvs[7].append(cnkey)
		elif nodeInTree < 0.250:
			myList_cnvs[8].append(cnkey)
		elif nodeInTree < 0.375:
			myList_cnvs[9].append(cnkey)
		elif nodeInTree < 0.500:
			myList_cnvs[10].append(cnkey)
		elif nodeInTree < 0.625:
			myList_cnvs[11].append(cnkey)
		elif nodeInTree < 0.750:
			myList_cnvs[12].append(cnkey)
		elif nodeInTree < 0.875:
			myList_cnvs[13].append(cnkey)
		else:
			myList_cnvs[14].append(cnkey)



## Now write the CNVs to files - one for each clone
cnt=0
for subcloneList in myList_cnvs:
	name0=myList_names[cnt]
	outName0gain = out_dir + '/cnv_gains_' + name0 +'.bed'
	outName0loss = out_dir + '/cnv_losses_' + name0 +'.bed'
	outFile0gain = open(outName0gain, 'w')
	outFile0loss = open(outName0loss, 'w')
	name1=name0.replace("_0_","_1_")
	outName1gain = out_dir + '/cnv_gains_' + name1 +'.bed'
	outName1loss = out_dir + '/cnv_losses_' + name1 +'.bed'
	outFile1gain = open(outName1gain, 'w')
	outFile1loss = open(outName1loss, 'w')
	for thisCNV in subcloneList:
		cnInfos = thisCNV.strip().split("_") ## chr + '_' + posStart + '_' + posEnd + '_' + copyNumber
		startPos=int(cnInfos[1])-1 # bed file has 0-based start
		if cnInfos[3] == "1": # loss
			copyStatus="loss"
		elif cnInfos[3] == "3": # gain
			copyStatus="gain"
		isHomozygous=random.random()
		if isHomozygous < 0.5: # heterozygous
			whichCopy=random.random()
			if whichCopy < 0.5:
				if copyStatus=="loss":
					outFile0loss.write('chr%s\t%d\t%s\t%s\n' % (cnInfos[0],startPos,cnInfos[2],copyStatus))
				else:
					outFile0gain.write('chr%s\t%d\t%s\t%s\n' % (cnInfos[0],startPos,cnInfos[2],copyStatus))
			else:
				if copyStatus=="loss":
					outFile1loss.write('chr%s\t%d\t%s\t%s\n' % (cnInfos[0],startPos,cnInfos[2],copyStatus))
				else:
					outFile1gain.write('chr%s\t%d\t%s\t%s\n' % (cnInfos[0],startPos,cnInfos[2],copyStatus))
		else: # homozygous
			if copyStatus=="loss":
				outFile0loss.write('chr%s\t%d\t%s\t%s\n' % (cnInfos[0],startPos,cnInfos[2],copyStatus))
				outFile1loss.write('chr%s\t%d\t%s\t%s\n' % (cnInfos[0],startPos,cnInfos[2],copyStatus))
			else:
				outFile0gain.write('chr%s\t%d\t%s\t%s\n' % (cnInfos[0],startPos,cnInfos[2],copyStatus))
				outFile1gain.write('chr%s\t%d\t%s\t%s\n' % (cnInfos[0],startPos,cnInfos[2],copyStatus))
	cnt+=1

