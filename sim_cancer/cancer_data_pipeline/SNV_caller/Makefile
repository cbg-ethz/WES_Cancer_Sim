all: deepSNVToVcf JointSNVMixToVcf


CXXFLAGS += -g -O3 #-std=c++0x
#CXXFLAGS += -g 

LIBS += bam z curses m pthread 
LDFLAGS += $(foreach library,$(LIBS),-l$(library))
LDFLAGS += -L${SAMDIR} -L${SAMDIR}/bcftools/
CPPFLAGS += $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir))

all: deepSNVToVcf JointSNVMixToVcf

deepSNVToVcf:deepSNVToVcf.cpp
	g++ -g deepSNVToVcf.cpp -o deepSNVToVcf ${LDFLAGS} ${CPPFLAGS}

JointSNVMixToVcf:JointSNVMixToVcf.cpp
	g++ JointSNVMixToVcf.cpp -o JointSNVMixToVcf
