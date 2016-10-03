
importedScriptsDir=../../eval_vcf/imported_scripts/

INCLUDE_DIRS=${importedScriptsDir}

CXXFLAGS += -g -std=c++0x

SAMDIR = /usr/local/beerenwinkel/samtools-0.1.19
BCF_INCL = /usr/local/beerenwinkel/samtools-0.1.19/
BCF_LIB = /usr/local/beerenwinkel/samtools-0.1.19/bcftools/


INCLUDE_DIRS += ${SAMDIR} ${importedScriptsDir}/hdf5/hdf5-1.8.11/hdf5//include/ ${importedScriptsDir}/ngs/ ${BCF_LIB}
LD_DIRS += ${SAMDIR} ${BCF_LIB} ${importedScriptsDir} ${importedScriptsDir}/hdf5/hdf5-1.8.11/hdf5//lib/
LIBS += bam z curses bcf m 

# bam_region
LIBS += hdf5 hdf5_cpp hdf5 hdf5_hl

CPPFLAGS += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))
LDFLAGS += $(foreach librarydir,$(LD_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(LIBS),-l$(library))

.PHONY: all clean

all: mmr_for_variant_calling

mmr_for_variant_calling: mmr_for_variant_calling.o
	$(CXX) -o mmr_for_variant_calling mmr_for_variant_calling.o  $(CPPFLAGS) $(LDFLAGS)
