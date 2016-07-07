#ifndef _GTF_TOOLS_H__
#define _GTF_TOOLS_H__
#include <assert.h>
#include "region.h"
#include <map>
	using std::map;

char* get_attribute(char* line, const char* tag);

vector<char*> get_fields(char* line);

vector<Region*> parse_gtf(char* gtf_file);

vector<Region*> parse_gtf(char* gtf_file, char* gene_name);

vector<Region*> parse_gff(char* gtf_file);

vector<Region*> parse_gff(char* gtf_file, const char* link_tag);


void write_gtf(FILE* fd, Region* region, const char* source);

void write_gff(FILE* fd, Region* region, const char* source);

vector<Region*> regions_from_map(map<string, Region*> transcripts);

vector<Region*> merge_overlapping_regions(vector<Region*> regions);

bool compare_second(segment s1, segment s2);

const char* determine_format(char* filename);

const char* determine_format(char* filename, const char* gff_link_tag);

#endif
