#ifndef _TOOLS_BAM__
#define _TOOLS_BAM__

#include "bam_region.h"
#include "tools.h"
#include "bam.h"

vector<vector<int> > region_overlap(vector<Bam_Region*> regions1, vector<Bam_Region*> regions2);

vector<Bam_Region*> parse_bam_regions(char* fn_regions);
vector<Bam_Region*> parse_bam_regions(char* fn_regions, char* chr_name);

void set_chr_num(Region* reg, bam_header_t* header);

#endif
