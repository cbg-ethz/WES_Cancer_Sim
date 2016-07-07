======================
= GemSIM version 1.2 =
======================

by Kerensa McElroy.

Copyright (c) 2011, Kerensa McElroy
kerensa@unsw.edu.au

LICENCE
=======

GemSIM is free software; it may be redistributed and modified 
under the terms of the GNU General Public License as published 
by the Free Software Foundation, either version 3 of the 
License, or (at your option) any later version.

GemSIM is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY, without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have recieved a copy of the GNU General Public
License along with GemSIM. If not, see 
http://www.gnu.org/licenses/. 


INTRODUCTION
============

GemSIM is a software package for generating realistic simulated 
next-generation sequencing reads with quality score values. Both 
Illumina (single or paired end) and Roche/454 reads  can be 
simulated using appropriate empirical error models. 


DESCRIPTION
===========

GemErr.py:

Takes a sam file and catalogues all the mismatches, insertions, and deletions
to create an error model for a particular sequencing run. Known true SNP
positions may be excluded.

Options:
      -h prints these instructions.
      -r read length. Set to LONGEST read in dataset.
      -f reference genome in fasta format
      -s input file in sam format.
      -n desired output filename prefix.
      -c T for circular reference, F for linear reference genome.
      -i use only every ith read for model (optional, must be odd).
      -m maximum indel size (optional, default=4).
      -p use only if your data contains paired end reads.
      -e comma separated list of reference positions to exclude e.g. '293, 342'

GemHaps.py:

Uses a reference genome to create a set of related haplotypes for 
input into GemReads.py. Alternatively, users may manually create 
their own haplotype input file (see manual). Haplotype frequency, 
and the number of SNPs in each haplotype are specificed by the user. 
For instance, to specify that you want to include two haplotypes, 
one identical to the reference with frequency 80%, and one with 15 
SNPs compared to the reference and frequency 20%, type '.80,0 .20,15' 
after the option -g. SNPs are then randomly placed along thelength 
of the genome.

NOTE: haplotypes MUST sum to 1!

Options:
      -h prints these instructions.
      -r reference genome, in fasta format.
      -g haplotype list. Format '.80,0 .20,15' (see above, and manual).
      -o output filename.


GemReads.py:

Takes a reference genome, an empirical error model, and a haplotype file
listing SNP locations and frequencies, and creates a simulated data set
of random reads, as would be produced by a next-gen sequencing run.
Output is in fastq format, suitable for input into popular alignment
software.

Options:
      -h prints these instructions.
      -r reference genome, in fasta format.
      -d Only for metagenome projects. Directory containing references.
      -a Only for metagenome projects. Species-abundance file.
      -n number of reads to produce. For paired end reads, number of pairs.
      -g haplotype file, specifying location and frequency of snps.
      -l length of reads. Integer value, or -l d for empirical distribution.
      -m error model file *_single.gzip or *_paired.gzip.
      -c T for circular reference, F for linear reference genome.
      -q quality score offset. Usually 33 or 64 (see manual).
      -o output file name prefix.
      -p use only to create paired end reads.

GemStats.py:

Takes error model files produce by GemErr.py, and generates statistics
for a particular error model. Output saved as .txt file.

Options:
      -h prints these instructions.
      -m error model file *_single.gzip or *_paired.gzip.
      -p use if model is for paired end reads.
      -n prefix for output files.


BUGS
====
please email kerensa@unsw.edu.au if you find one!


CHANGE LOG
==========

Changes since Version 1.0:

1.6:
* improved handling of multi-chromosome references
* added normal distribution option for fragment lengths
* added fragment length to read header for paired end reads
* added support for genotype directory in metagenomics mode
* changed command line option for metagenomic reference directory to -R
* changed minimum k-mer default option to 0 (required for speed for large genomes)

1.5:
* corrected but in GemReads.py that caused problems when reference genome
  featured lower case letters

1.4:
* corrected bug in GemReads.py that caused problems parsing input genomes
  with ambiguous letters

1.3:
* users now just specify a circular genome by supplying -c; no true/false
  argument is required.
* the minimum number of times a k-mer is required to be present in the 
  reference genome when calculating error models with GemErr now a 
  user-specified input parameter. 

1.2:
* fixed bug with extracting reads from linear genomes.
* ambiguous characters in reference are now replaced by 'N'.

1.1:
* added support for metagenomic projects.
* added support for linear genomes.
* fixed several bugs.
