#ifndef __GENOME_H__
#define __GENOME_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


class Genome
{
	public:
		Genome();
		~Genome();
	int num_contigs ;
	char **contig_fnames;
	char **contig_fasta_fnames;
	char **contig_names;
	FILE **contig_files ;
	char *alph ;
	
	int num_est_files ;
	char **est_fnames ;
	
	int num_prot_files ;
	char **prot_fnames ;
	
	int num_cdna_files ;
	char **cdna_fnames ;
	
	char *basedir ;
	
	int num_annotation_files ;
	char **annotation_fnames ;
	
	int init_genome(char *fname) ;
	//void cleanup_genome() ;
	int get_contig_num(char *contig_name) ;
	char* get_contig_name(int chr_num);
	char* path_with_basedir(char* b1, char* basedir);
	char* read_line(FILE*infile, char* line, size_t max_len);
	char* read_flat_file(int chr_num);
	char* read_flat_file(int chr_num, int start, int stop);
	int contig_len(int chr_num);
	void complement(char* str, const int len);
	
	
	 int MAX_CDNA_DELETION ;
	 int MAX_CDNA_INSERTION ;
	 int MIN_INTRON_LEN     ;
	 int MAX_INTRON_LEN     ;
	 int MAX_DOUBLE_QGAP    ;
	 int MAX_DOUBLE_TGAP    ;
	 int TERMINAL_EXON_END_TOL  ;
	 int MERGE_EST_TRANSCRIPTS ;
	 int MAX_PERFECT_EST_MATCH_LEN ;
	 int BLAT_BEST_HIT_ONLY ;
	 int MIN_SSQUALITY_SCORE ;
	 int MIN_EXONQUALITY_SCORE ;
	
	 double MIN_EST_COVER_FRAC  ;
	 double MIN_PROT_COVER_FRAC  ;
	 double MIN_CDNA_COVER_FRAC  ;
	 double MAX_END_MISSING     ;
	 double MAX_START_MISSING   ;
	 double MAX_END_MISSING2    ;
	 double MAX_START_MISSING2  ;
	 double BLAT_BEST_HIT_MARGIN ;

};
#endif
