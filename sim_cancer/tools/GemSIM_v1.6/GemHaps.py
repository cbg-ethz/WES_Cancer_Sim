#!/usr/bin/python

# Copyright (C) 2011, Kerensa McElroy.
# kerensa@unsw.edu.au

# This file is part of the sequence simulator GemSIM. 
# It is used to generate haplotypes for input into 
# the read generator GemReads.py. 
# Haplotypes are based on a reference sequence, with 
# randomly introduced SNPs.  

# GemSIM is free software; it may be redistributed and 
# modified under the terms of the GNU General Public 
# License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option)
# any later version.

# GemSIM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, without even the implied 
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU General Public License for more 
# details.

# You should have recieved a copy of the GNU General Public
# License along with GemSIM. If not, see 
# http://www.gnu.org/licenses/. 

import sys
import getopt
import random
import bisect

def bisect_choiceTUP(items):
    """Returns a function that makes a weighted random choice from a dictionary of tuples."""
    added_weights = []
    last_sum = 0.0
    for chrom,ref,lnth in items:
        weight=float(lnth)
        last_sum += weight
        added_weights.append(last_sum)
    def choice(rnd=random.random, bis=bisect.bisect):
        return items[bis(added_weights, rnd() * last_sum)]
    return choice

def gParse(gT,ref):
    gL=[]
    chr=bisect_choiceTUP(ref)
    for g,s in gT:
        x=0
        gD={}
        while x<int(s):
            chrom,ref,lnth=chr() 
            pos=random.randint(0,lnth-1)
            nuc=ref[pos]
            if nuc=='A' or nuc=='a':
                newNuc=random.choice('TGC')
            elif nuc=='T' or nuc=='t':
                newNuc=random.choice('AGC')
            elif nuc=='C' or nuc=='c':
                newNuc=random.choice('ATG')
            elif nuc=='G' or nuc=='g':
                newNuc=random.choice('ATC')
            if chrom in gD:
                if (pos+1) not in gD[chrom]:
                    gD[chrom][pos+1]=newNuc
                    x+=1
            else:
                gD[chrom]={pos+1:newNuc}
                x+=1
        gL.append((float(g),gD)) 
    return gL

def getRef(refFile):
    """Returns a list of chromosomes."""
    ref=''
    refList=[]
    f=open(refFile)
    head=f.readline().rstrip()[1:]
    line=f.readline()
    while line:
        if line[0]!='>':
            ref+=line.rstrip()
        else:
            ref=ref.upper()
            for l in 'RYLMKSWHBVD':
                ref=ref.replace(l,'N')
            refList.append((head,ref,len(ref)))
            ref=''
            head=line.rstrip()[1:]
        line=f.readline()
    refList.append((head,ref,len(ref)))
    return refList 

def usage():
    print '\n########################################################################'
    print '# GemSIM - Generic Error Model based SIMulator of N.G. sequencing data #'
    print '########################################################################\n'
    print '\nGemHaps.py:\n'
    print 'Uses a reference genome to create a set of related haplotypes for input into'
    print 'GemReads.py. Alternatively, users may manually create their own haplotype input'
    print 'file (see manual). Haplotype frequency, and the number of SNPs in each'
    print 'haplotype are specificed by the user. For instance, to specify that you want'
    print 'to include two haplotypes, one identical to the reference with frequency 80%,'
    print 'and one with 15 SNPs compared to the reference and frequency 20%, type '
    print "'.80,0 .20,15' after the option -g. SNPs are then randomly placed along the"
    print 'length of the genome.'
    print '\nNOTE: haplotypes MUST sum to 1!'
    print '\nOptions:'
    print '      -h prints these instructions.'
    print '      -r reference genome, in fasta format.'
    print "      -g haplotype list. Format '.80,0 .20,15' (see above, and manual).\n\n" 

def main(argv):
    ref=''
    gens=''
    try:
        opts, args = getopt.getopt(argv, "hr:g:o:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt =='-h':
            usage()
            sys.exit()
        elif opt =='-r':
            ref=arg
        elif opt =='-g':
            gens=arg
    if ref=='' or gens=='':
        usage()
        sys.exit(2)

    #checks haplotypes sum to 1
    Gs=gens.split(' ')
    gT=[]
    sum=0
    for g in Gs:
        g,s=g.split(',')
        gT.append((g,s))
        sum+=float(g)
    if str(sum)!='1.0':
        print '\nERROR!'
        print 'Genotypes do not sum to 1! Remember to add a genotype with 0 snps'
        print 'to include a genotype identical to the reference.\n'
        sys.exit()
    stem=ref.split('/')[-1]
    stem=stem.split('.')[0]
    ref = getRef(ref)
    gL=gParse(gT,ref)
    o=open(stem+'.txt', 'w')
    for g,d in gL:
        line=''
        line+=str(g)    
        keys=d.keys()
        keys.sort()
        for k in keys:
            snps=d[k].keys()
            snps.sort()
            for s in snps:
                line+=('\t'+k+'\t'+str(s)+'\t'+d[k][s])
        o.write(line+'\n')
    o.close()

if __name__=="__main__":
    main(sys.argv[1:])
