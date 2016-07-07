import sys
import math
import random
import string

# funtion to go up in a binary tree
def goUp(pos):
    return int((pos + 1) / 2) - 1

# funtion to go to th left node in a binary tree
def goLeft(pos):
    return 2 * (pos + 1) - 1

# funtion to go to the right node in a binary tree
def goRight(pos):
    return 2 * (pos + 1)

def sample_freq():
    return(random.random())


# check if two overlapping variants support each other
def checkOverlapCollision(pos, line):
    leavePos=pos[0]

    for i in range(1,3):
        if (pos[1] == 0) or (pos[1] == i): # if homo two checks, otherwise one
            parent = i - 1;

            lineChrom = line.split("\t")[0]
            linePos = int(line.split("\t")[1])
            lineRefLen = len(line.split("\t")[3])
            lineAlt = line.split("\t")[4]
            lineAltLen = len(lineAlt)

            if (len(vcfs[leavePos][parent]) > 0): # do not insert a variant at the same position


                vcfLineLast=vcfs[leavePos][parent][-1]
                vcfChromLast = vcfLineLast.split("\t")[0]
                vcfPosLast = int(vcfLineLast.split("\t")[1])

                if (lineChrom == vcfChromLast) and (linePos == vcfPosLast):
                    #print line
                    #print vcfLineLast
                    return True

            if lineChrom == overlapCheck[leavePos][parent][0]:
                vcfLine=vcfs[leavePos][parent][overlapCheck[leavePos][parent][1]]
                vcfPos = int(vcfLine.split("\t")[1])
                vcfRefLen = len(vcfLine.split("\t")[3])
                vcfAlt = vcfLine.split("\t")[4]
                vcfAltLen = len(vcfAlt)

                if linePos < overlapCheck[leavePos][parent][2]:
                    return True # if this is commented out supporting mutations are allowed
                    vcfSubStart = linePos - vcfPos
                    #check if overlap or contained
                    if (vcfRefLen - vcfSubStart + 1 < lineRefLen):
                        overlapEnd = vcfRefLen - vcfSubStart
                        if vcfAlt[vcfSubStart:] != lineAlt[:overlapEnd]:
                            return True
                    else:
                        if vcfAlt[vcfSubStart:vcfSubStart+lineAltLen] != lineAlt:
                            return True
    return False

def checkOverlapCollisionDownstream(pos, line):
    if goLeft(pos[0]) >= len(vcfs):
        return checkOverlapCollision(pos, line)
    else:
        return checkOverlapCollisionDownstream((goLeft(pos[0]), pos[1]), line) or checkOverlapCollisionDownstream((goRight(pos[0]), pos[1]), line)

#check in the lower levels of the tree if newPos and posToCheck are in the same sub-tree
def checkDownstream(newPos, posToCheck, line):
    if newPos[0] == posToCheck[0]:
        if ((posToCheck[1] == 1) and (newPos[1] == 2)) or ((posToCheck[1] == 2) and (newPos[1] == 1)):
            return False
        return True
    elif newPos[0] > len(modProbs) -1:
        return False
    return checkDownstream((goLeft(newPos[0]), newPos[1]), posToCheck, line) or checkDownstream((goRight(newPos[0]), newPos[1]), posToCheck, line)

# place a mutation into all nodes of the subree with node pos as root
def placeMutationInSubtree(line, pos):
    if insertNormal and pos[0] == 0: # the "and pos[0]" is necessary because of the recursive call
        if pos[1] == 0:
            vcfsNormals[0].append(line)
            vcfsNormals[1].append(line)
        elif pos[1] == 1:
            vcfsNormals[0].append(line)
        else:
            vcfsNormals[1].append(line)

    if float(len(probs) + 1) > 2.0 * float(pos[0] + 1):
        placeMutationInSubtree(line, (goLeft(pos[0]), pos[1]))
        placeMutationInSubtree(line, (goRight(pos[0]), pos[1]))
        return

    if pos[0] < len(probs):
        lineChrom = line.split("\t")[0]
        linePos = int(line.split("\t")[1])
        lineLen = len(line.split("\t")[3])

        for i in range(1, 3):

            if (pos[1] == 0) or (pos[1] == i):
                parent = i - 1
                vcfs[pos[0]][parent].append(line)

                # keep track of variant to check for overlap collisions
                if (lineChrom != overlapCheck[pos[0]][parent][0]) or (linePos + lineLen > overlapCheck[pos[0]][parent][2]):
                    if parent == 0:
                        newTuple = ((lineChrom, len(vcfs[pos[0]][parent]) - 1, linePos + lineLen), overlapCheck[pos[0]][1])
                    else :
                        newTuple = (overlapCheck[pos[0]][0], (lineChrom, len(vcfs[pos[0]][parent]) - 1, linePos + lineLen))
                    overlapCheck[pos[0]] = newTuple

# determine in which node in the tree a certain mutation is placed
def placeMutation(line, tryMutationNo, placeSuccessNo):
    if insertNormal:
        freq = 1.0
    else:
        freq = sample_freq()
    for i in range(0, len(treeLevelRange)):
        if freq >= treeLevelRange[i]:
            rand=random.random() # node in the tree level
            parent=random.random() # homo or hetero

            for j in range(0, len(cygocity[insertNormal])):
                if parent >= cygocity[insertNormal][j]:
                    parent=j
                    break
            for j in range(pow(2, i)-1, pow(2, i+1) - 1):
                if rand <= modProbs[j]:
                    # just for statistics
                    if insertNormal:
                        #continue
                        tryMutationNo = tryMutationNo + 1
                    else:
                        tryMutationTu[j] = tryMutationTu[j] + 1

                    if not checkOverlapCollisionDownstream((j, parent), line):
                        # just for statistics
                        if insertNormal:
                            placeSuccessNo = placeSuccessNo + 1
                        else:
                            placeSuccessTu[j] = placeSuccessTu[j] + 1

                        placeMutationInSubtree(line, (j, parent))
                        return (tryMutationNo, placeSuccessNo)
                    else :
                        vcfFails.append(line)
                    break
            return (tryMutationNo, placeSuccessNo)

'''
Here starts the program
'''

# get the command line arguments [1] = tumor vcf, [2] = normal vcf, [3] = configuration file
vcfInTumorName=sys.argv[1]
vcfInNormalName=sys.argv[2]
configInName=sys.argv[3]

random.seed(44)

#arrays for statistics
tryMutationNo = 0
tryMutationTu = []
placeSuccessNo = 0
placeSuccessTu = []

for i in range(15):
    tryMutationTu.append(0)
    placeSuccessTu.append(0)

# read config file
entry=-1
with open(configInName, 'r') as configIn:
    # treeLevelRange stores the frequency cut off that determine in which level of the tree a vcf record is placed
    treeLevelRange = []
    # prob stores the probability to put a vcf record into a certain node with in a tree level
    probs = []
    # homo, father, mother
    cygocity = ([], [])
    for line in configIn:
        if not line.startswith('#'):
            if line.startswith('>treeLevelRange'):
                entry=0
            elif line.startswith('>treeStructure'):
                entry=1
            elif line.startswith('>cygocityT'):
                entry=2
            elif line.startswith('>cygocityN'):
                entry=3
            elif entry == 0:
                treeLevelRange.extend([float(i) for i in line.strip().split(" ")])
            elif entry == 1:
                line=[float(i) for i in line.strip().split(" ")]
                for i in line:
                    j = float(i) / sum(line)
                    probs.append(j)
            elif entry == 2:
                cygocity[0].extend([float(i) for i in line.strip().split(" ")])
            elif entry == 3:
                cygocity[1].extend([float(i) for i in line.strip().split(" ")])

# modProbes are the modified probs that gurantee that the sum of probs in each tree level sum to 1
modProbs=probs
tmp=0;
for i in range(0, len(modProbs)):
    if i+1 == pow(2, int(math.log(i+1, 2))):
        tmp=modProbs[i]
    else:
        tmp+=modProbs[i]
        modProbs[i]=tmp

print "THE CHROMOSOMES OF THE TWO VCF FILES HAVE TO BE SORTED LEXICOGRAPHICALLY!"

print modProbs

for i in range(0, len(treeLevelRange)):
    treeLevelRange[i] = float(len(treeLevelRange) - 1 - i)/float(len(treeLevelRange))

print "TreeLevelRange: ", treeLevelRange

# initialize the vcf files that we create
vcfs = []
overlapCheck = []
for i in range(0, len(probs)):
    vcfs.append(([], []))
    overlapCheck.append((("",-1,-1), ("",-1,-1)))

vcfsNormals = ([], [])

vcfFails = []

vcfInTumor = open(vcfInTumorName, 'r')
vcfInNormal = open(vcfInNormalName, 'r')

# parse the vcf records and call the functions above
header = []
counter = 0

# ignore the header of the tumor vcf file
lineTumor = vcfInTumor.readline();
while lineTumor.startswith('#'):
    lineTumor = vcfInTumor.readline();

# read the header of the normal vcf file
lineNormal = vcfInNormal.readline();
while lineNormal.startswith('#'):
    header.append(lineNormal)
    lineNormal = vcfInNormal.readline();

insertNormal = False
doublePos=False
doubleCounter = 0
while (lineNormal != "") or (lineTumor != ""):
    doublePos = False

    # skip half of the variants randomly
    #if (random.random()>0.5) and (lineNormal != ""):
    #    lineNormal = vcfInNormal.readline()

    #if (random.random()>0.5) and (lineTumor != ""):
    #    lineTumor = vcfInTumor.readline()

    if (lineNormal == "") and (lineTumor == ""):
        break
    elif lineNormal == "":
        line=lineTumor
        insertNormal = False
        lineTumor = vcfInTumor.readline()
    elif lineTumor == "":
        line=lineNormal
        insertNormal = True
        lineNormal = vcfInNormal.readline()
    else:
        tumorChrom = lineTumor.split("\t")[0]
        normalChrom = lineNormal.split("\t")[0]
        tumorPos = int(lineTumor.split("\t")[1])
        normalPos = int(lineNormal.split("\t")[1])
        if tumorChrom == normalChrom:
            if tumorPos == normalPos:
                doublePos = True
                line = lineTumor
                insertNormal = False
                lineTumor = vcfInTumor.readline()
                doubleCounter += 1
                print "Double position should not have happend!"

            elif tumorPos < normalPos:
                line = lineTumor
                insertNormal = False
                lineTumor = vcfInTumor.readline()
            else:
                line = lineNormal
                insertNormal = True
                lineNormal = vcfInNormal.readline()
        elif tumorChrom < normalChrom:
            line = lineTumor
            insertNormal = False
            lineTumor = vcfInTumor.readline()
        else:
            line = lineNormal
            insertNormal = True
            lineNormal = vcfInNormal.readline()

    if not doublePos:
        counter = counter + 1
        if counter % 100000 == 0:
            print counter
        a0=line.split("\t")[7].split(";")[5].split("=")[1]
        numAlt=len(a0.split(','))
        if numAlt > 1:
            positions=[]
            for altCount in range(0, numAlt):
                newVar=""
                tab=line.split("\t")
                for i in [0, 1, 2, 3]:
                    newVar=newVar+tab[i]+"\t"
                newVar=newVar+tab[4].split(',')[altCount]+"\t"
                for i in [5, 6]:
                    newVar=newVar+tab[i]+"\t"
                sem=tab[7].split(";")
                for i in range(0, len(sem)-1):
                    if sem[i].find(",") == -1:
                        newVar=newVar+sem[i]+';'
                    else:
                        eq=sem[i].split('=')
                        com=eq[1].split(',')
                        newVar=newVar+eq[0]+'='+com[altCount]+';'

                newVar=newVar+sem[-1]+"\t"
                newVar=newVar+tab[8]+"\t"
                col=tab[9].split(':')
                for i in [0, 1, 2, 3]:
                    newVar=newVar+col[i]+':'
                newVar=newVar+col[4].split(',')[altCount]+':'+col[5].split(',')[altCount]+':'
                com=col[6].split(',')
                for i in range(0,len(com)/numAlt-1):
                    newVar=newVar+com[i*numAlt+altCount]+','
                newVar=newVar+com[(len(com)/numAlt-1)*numAlt+altCount]
                if newVar[len(newVar) - 1] != '\n':
                    newVar=newVar + '\n'

                (tryMutationNo, placeSuccessNo) = placeMutation(newVar, tryMutationNo, placeSuccessNo)
        else:
            (tryMutationNo, placeSuccessNo) = placeMutation(line, tryMutationNo, placeSuccessNo)

print "tryMutation"
tmp=0;
tryMutationTuSum=[]
for i in range(0, len(tryMutationTu)):
    if i+1 == pow(2, int(math.log(i+1, 2))):
        tryMutationTuSum.append(tmp)
        tmp=tryMutationTu[i]
    else:
        tmp+=tryMutationTu[i]
tryMutationTuSum.append(tmp)
print tryMutationNo
print tryMutationTuSum
print tryMutationTu

print "placeSuccess"
tmp=0;
placeSuccessTuSum=[]
for i in range(0, len(placeSuccessTu)):
    if i+1 == pow(2, int(math.log(i+1, 2))):
        placeSuccessTuSum.append(tmp)
        tmp=placeSuccessTu[i]
    else:
        tmp+=placeSuccessTu[i]
placeSuccessTuSum.append(tmp)
print placeSuccessNo
print placeSuccessTuSum
print placeSuccessTu

print "doubleCounter"
print doubleCounter

# write the result to disk
for i in range(int(len(vcfs)/2), len(vcfs)):
    for j in range(0,2):
        outName=str(i)+"_"+str(j)+"_"+vcfInTumorName
        with open(outName, 'w') as out:
            for entry in header:
                out.write(entry);
            for entry in vcfs[i][j]:
                out.write(entry);
            out.close()

for j in range(0,2):
    outName=str(j)+"_"+vcfInNormalName
    with open(outName, 'w') as out:
        for entry in header:
            out.write(entry);
        for entry in vcfsNormals[j]:
            out.write(entry);
        out.close()

outName="Fails.vcf"
with open(outName, 'w') as out:
    for entry in header:
        out.write(entry);
    for entry in vcfFails:
        out.write(entry);
    out.close()
