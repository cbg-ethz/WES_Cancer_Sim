#CXX=clang++
CXX=g++


# Change this path to your path to samtools:
PATHTOSAMTOOLS=/usr/local/beerenwinkel/samtools-0.1.19/


CXXFLAGS+=-g -Wall -W -Wno-long-long -pedantic -Wno-variadic-macros
CXXFLAGS+=-O3 
CXXFLAGS += -std=c++0x

BCF_INCL = ${PATHTOSAMTOOLS}


BCF_O = $(BCF_INCL)/kstring.o $(BCF_INCL)/bgzf.o $(BCF_INCL)/knetfile.o $(BCF_INCL)/bedidx.o 

LIBS += pthread bam z curses bcf m 


INCLUDE_DIRS = ./imported_scripts/
LD_DIRS = ./imported_scripts/ ${PATHTOSAMTOOLS} ${PATHTOSAMTOOLS}/bcftools/


CXXFLAGS += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
LDFLAGS += $(foreach librarydir,$(LD_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(LIBS),-l$(library))

default: all
all: evaluate_vcf select_true_alignment filter_somatic_binom_test_JointSNVMix2 

##all: evaluate_vcf

evaluate_vcf.o: evaluate_vcf.cpp evaluate_vcf.h

evaluate_vcf: evaluate_vcf.o evaluate_vcf.h evaluate_vcf.cpp
	$(CXX) -o evaluate_vcf evaluate_vcf.o ${BCF_O} $(LDFLAGS) $(CXXFLAGS)

evaluate_vcf_pair: evaluate_vcf_pair.o
	$(CXX) -o evaluate_vcf_pair evaluate_vcf_pair.o ${BCF_O} $(LDFLAGS)

select_true_alignment: select_true_alignment.o
	$(CXX) -o select_true_alignment select_true_alignment.o  $(CPPFLAGS) $(LDFLAGS) $(CXXFLAGS)

filter_somatic_binom_test_JointSNVMix2: filter_somatic_binom_test_JointSNVMix2.cpp
	$(CXX) -o filter_somatic_binom_test_JointSNVMix2 filter_somatic_binom_test_JointSNVMix2.cpp

clean:
	@- $(RM) evaluate_vcf evaluate_vcf.o evaluate_vcf_pair.o select_true_alignment.o select_true_alignment evaluate_vcf_pair

print_vars:
	echo $(CXXFLAGS)

.PHONY: default all clean
