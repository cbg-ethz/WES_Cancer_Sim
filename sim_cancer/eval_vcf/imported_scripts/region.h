
#ifndef _REGIONS_H__
#define _REGIONS_H__

//#define READ_DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
//#include "genome.h"
#include <utility>
	using std::pair;
#include <vector>
	using std::vector;
#include <fstream>
#include <string>
	using std::string;

//typedef pair<int, int> segment;
class segment{
	public:
		segment(){first = 0; second = 0; flag = -1;};
		segment(int f, int s){first = f; second = s; flag=-1;};
		segment(int f, int s, int fl){first = f; second = s; flag=fl;};
		int first;
		int second;
		int flag;//4:CDS, 3: 3'UTR 5:5'UTR

		friend bool operator== (const segment& lhs, const segment& rhs)
		{return lhs.first==rhs.first && lhs.second==rhs.second;}

		friend bool operator<  (const segment& lhs, const segment& rhs)
		{return lhs.first<rhs.first || (!(rhs.first<lhs.first) && lhs.second<rhs.second);}

};
#define NO_CONNECTION -2
#define NEIGHBOR -1
#define CONNECTION 0

//#define INNER=0
//#define INITIAL=1
//#define TERMINAL=2

class Region
{
	public:
		int id;
		int start; 
		int stop;
		int chr_num;
		char* chr;
		char strand; 
		vector<segment> intron_list;
		vector<segment> unique_introns;
		vector<segment> segments;
		vector<vector<segment> > transcripts;
		vector<string> transcript_names;
		vector<string> gene_names;
		vector<vector<int> > transcript_paths;
		FILE* fd_out;

		//Genome* gio; 
			
		/** default constructor*/	
		Region();

		/** constructor*/
		Region(int pstart, int pstop, int pchr_num, char pstrand, const char* gio_fname);
		Region(int pstart, int pstop, int pchr_num, char pstrand);
		Region(int pstart, int pstop, char* chr, char pstrand);
		Region(Region* reg);

		/** destructor*/
		~Region();


		//char* get_sequence()
		//{
		//	if (!seq)
		//		load_genomic_sequence();
		//	return seq;
		//};

		//bool check_region()
		//{
		//	if (!(start>=0 && stop<gio->contig_len(chr_num)))
		//	{
		//		fprintf(stderr, "chr_num: %i, chr_len: %i\n", chr_num, gio->contig_len(chr_num));
		//	}
		//	return (start>=0 && stop<gio->contig_len(chr_num));	
		//};


		///** get the triplet codon for pos in transcript i */
		//char* get_triplet(int i, int pos, int* offset);

		///** find translation initiation site for transcript i*/
		//int get_TIS(int i); 
		///** check translation initiation consensus for transcript i*/
		//int check_TIS(int i); 

		///** find translation termination site for transcript i */
		//int get_TTS(int i); 
		///** check translation termination consensus for transcript i */
		//int check_TTS(int i); 

		//void set_gio(Genome* pgio);
	
		//void load_genomic_sequence(); 			

		char* get_region_str();

		virtual void print(_IO_FILE*& fd);

		void write_bed(FILE* fd)
		{
			fprintf(fd, "%s\t%i\t%i\n", chr, start, stop); 
		}

	private: 
		char* seq;

};

Region::Region()
{
	start = -1; 
	stop = -1;
	strand = '\0';
	chr_num = -1;
	chr = NULL;
	seq = NULL;
	//gio = NULL;
}


///** constructor*/
//Region::Region(int pstart, int pstop, int pchr_num, char pstrand, const char* gio_fname)
//{
//	start = pstart; 
//	stop = pstop;
//	strand = pstrand;
//	chr_num = pchr_num;
//	chr = NULL;
//	seq = NULL;
//	fd_out = stdout;
//
//	// initialize genome information object
//	gio = new Genome(); 
//	int ret = gio->init_genome((char*) gio_fname);
//	if (ret<0)
//	{   
//		fprintf(stderr, "error reading genome info file: %s\n", gio_fname);
//		return;
//	}
//}

/** constructor*/
Region::Region(int pstart, int pstop, int pchr_num, char pstrand)
{
	start = pstart; 
	stop = pstop;
	strand = pstrand;
	chr_num = pchr_num;
	chr = NULL;
	seq = NULL;
	fd_out = stdout;
	//gio = NULL;
}

/** constructor*/
Region::Region(int pstart, int pstop, char* pchr, char pstrand)
{
	start = pstart; 
	stop = pstop;
	strand = pstrand;
	chr_num = -1;
	chr = new char[strlen(pchr)+1];
	strcpy(chr, pchr);
	seq = NULL;
	fd_out = stdout;
	//gio = NULL;
}

/** constructor*/
Region::Region(Region* reg)
{
	start = reg->start; 
	stop = reg->stop;
	strand = reg->strand;
	chr_num = reg->chr_num;
	if (reg->chr)
	{
		chr = new char[strlen(reg->chr)+1];
		strcpy(chr, reg->chr);
	}
	else
	{
		chr = NULL;
	}
	seq = NULL;
	fd_out = reg->fd_out;
	//gio = NULL;
	segments = reg->segments;
	transcripts = reg->transcripts;
	transcript_names = reg->transcript_names;
	transcript_paths = reg->transcript_paths;
}

Region::~Region()
{
	delete[] seq;	
	delete[] chr;
	intron_list.clear();
	unique_introns.clear();
}

//void Region::set_gio(Genome* pgio)
//{
//	gio = pgio;
//}

//char* Region::get_triplet(int i, int pos, int* offset)
//{
//	assert(i<(int)transcripts.size());
//	assert(i>=0);
//
//	// TODO: this is just for debugging
//	//pos = get_TIS(i);
//
//	//if (pos==0)
//	//	return NULL;
//
//	// the sequence in not reverse complemented
//	char* seq = get_sequence();
//	int seq_len=stop-start+1;
//
//	bool exonic=false;
//	int len=0;
//	int cds_pos=0;
//	for (unsigned int j=0; j<transcripts[i].size(); j++)	
//	{
//
//		int flag = transcripts[i][j].flag; 
//		int exon_start = transcripts[i][j].first;
//		int exon_stop = transcripts[i][j].second;
//
//		if (flag==4)// CDS
//		{
//			len += exon_stop-exon_start+1; 
//			if (exon_stop<pos)
//			{
//				cds_pos += exon_stop-exon_start+1;
//			}
//			else if (exon_start<=pos && exon_stop>=pos)
//			{
//				cds_pos += pos-exon_start;
//				exonic = true;
//
//				*offset = cds_pos%3;
//				int local=pos-start-*offset; 
//
//				char* ret = new char[4]; 
//				ret[3] = 0;
//
//				if (local<0)
//				{
//					printf("local: %i<0; pos:%i, start:%i stop:%i offset:%i\n", local, pos, start, stop, *offset);
//					return NULL;
//				}
//				assert(local<seq_len);
//
//				if (strand=='+')
//				{
//					// check if codon is split between two exons
//					if (pos-exon_start<*offset)
//					{
//						// get the rest of the codon from previous exon
//						assert(j>0);
//						int pflag = transcripts[i][j-1].flag;
//						int pexon_stop = transcripts[i][j-1].second-start;
//						assert(pflag==4); 
//						
//						if (pos-exon_start<*offset-1)
//						{
//							ret[0] = seq[pexon_stop-1];
//							ret[1] = seq[pexon_stop];
//							ret[2] = seq[local+2];
//						}
//						else
//						{
//							ret[0] = seq[pexon_stop];
//							ret[1] = seq[local+1];
//							ret[2] = seq[local+2];
//						}
//					}
//					else if (exon_stop-(pos-*offset)<2)
//					{
//						// get the rest of the codon from next exon
//						assert(j+1<transcripts[i].size());
//						int nflag = transcripts[i][j+1].flag;
//						int nexon_start = transcripts[i][j+1].first-start;
//						if(nflag!=4)
//						{
//							//printf("chr:%s pos:%i trans_num:%i exon:%i start:%i, end:%i\n", chr, pos, i, j, transcripts[i][j].first, transcripts[i][j].second);
//							//if (transcript_names.size()>i)
//							//	printf("transcript_name: %s\n", transcript_names[i].c_str()); 
//							//printf("Checked a few casses and found annotation bugs\n");  
//							printf("[%s] Warning, this is a split codon, but the next exon of this transcript is not annotated as a coding exon\n", __func__); 
//							return NULL;
//						}
//
//						if (exon_stop-(pos-*offset)<1)
//						{
//							ret[0] = seq[local];
//							ret[1] = seq[nexon_start];
//							ret[2] = seq[nexon_start+1];
//						}
//						else
//						{
//							ret[0] = seq[local];
//							ret[1] = seq[local+1];
//							ret[2] = seq[nexon_start];
//						}
//
//					}
//					else
//					{
//						ret[0] = seq[local];
//						ret[1] = seq[local+1];
//						ret[2] = seq[local+2];
//					}
//					//printf("(+)triplet: %s, offset:%i char:%c\n", ret, *offset, seq[local+*offset]);
//					return ret;
//				}
//				else
//				{
//					// reverse
//					ret[0] = seq[local+2];
//					ret[1] = seq[local+1];
//					ret[2] = seq[local];
//					// complement
//					gio->complement(ret, 4);
//					//printf("(-)triplet: %s, offset:%i char:%c\n", ret, *offset, seq[local+*offset]);
//					// because of the reverse
//					*offset = 2-*offset;
//					return ret;
//				}
//			}
//			else 
//			{
//				assert(exon_start>pos);
//			}
//		}
//	}
//	return NULL;
//}
//
//int Region::get_TTS(int i)
//{
//	assert(i<(int)transcripts.size());
//	assert(i>=0);
//
//	int pos=0;
//	for (unsigned int j=0; j<transcripts[i].size(); j++)
//	{
//		int flag = transcripts[i][j].flag;
//		if (flag==4)//CDS
//		{
//			if (strand=='+')
//			{
//				// tis is the start of the first coding exon
//				pos = transcripts[i][j].second;
//			}
//			else
//			{
//				// tis is the end of the last coding exon
//				pos = transcripts[i][j].first;
//				break;
//			}
//		}
//	}	
//	return pos;
//}
//
//int Region::get_TIS(int i)
//{
//	assert(i<(int)transcripts.size());
//	assert(i>=0);
//
//	int pos=0;
//	for (unsigned int j=0; j<transcripts[i].size(); j++)
//	{
//		int flag = transcripts[i][j].flag;
//		if (flag==4)//CDS
//		{
//			if (strand=='+')
//			{
//				// tis is the start of the first coding exon
//				pos = transcripts[i][j].first;
//				break;
//			}
//			else
//			{
//				// tis is the end of the last coding exon
//				pos = transcripts[i][j].second;
//			}
//		}
//	}
//	return pos;
//}
//
//int Region::check_TIS(int i)
//{
//	assert(i<(int)transcripts.size());
//	assert(i>=0);
//
//	int pos = get_TIS(i);
//	if (pos==0)
//		return -1;
//
//	// transform to local coordinates
//	pos -= start;
//
//	char* seq = get_sequence();
//	//int seq_len = strlen(seq);
//	int seq_len=stop-start+1;
//	if (pos<0)
//	{
//		printf("check_TIS: Warning: found pos: %i<0\n", pos);
//		return -2;
//	}
//	assert(pos<seq_len);
//
//
//	if (strand=='+')
//	{
//		//printf("TIS(%c): %c%c%c %c%c%c %c%c%c\n", strand, seq[pos-3], seq[pos-2], seq[pos-1], seq[pos], seq[pos+1], seq[pos+2], seq[pos+3], seq[pos+4], seq[pos+5]);
//		return (seq[pos]=='A' || seq[pos]=='a') && (seq[pos+1]=='T' || seq[pos+1]=='t') && (seq[pos+2]=='G' || seq[pos+2]=='g');
//	}
//	else
//	{
//		//printf("TIS(%c): %c%c%c %c%c%c %c%c%c\n", strand, seq[pos-5], seq[pos-4], seq[pos-3], seq[pos-2], seq[pos-1], seq[pos], seq[pos+1], seq[pos+2], seq[pos+3]);
//		return (seq[pos]=='T' || seq[pos]=='t') && (seq[pos-1]=='A' || seq[pos-1]=='a') && (seq[pos-2]=='C' || seq[pos-2]=='c');
//	}
//}
//
//int Region::check_TTS(int i)
//{
//	int pos = get_TTS(i);
//	if (pos==0)
//		return -1;
//
//	// transform to local coordinates
//	pos -= start;
//
//	char* seq = get_sequence();
//	//int seq_len = strlen(seq);
//	int seq_len=stop-start+1;
//	if (pos<0)
//	{
//		printf("check_TTS: Warning: found pos: %i<0\n", pos);
//		return -2;
//	}
//
//	assert(pos<seq_len);
//
//
//	if (strand=='+')
//	{
//		//printf("TTS(%c): %c%c%c %c%c%c %c%c%c\n", strand, seq[pos-3], seq[pos-1], seq[pos], seq[pos+1], seq[pos+2], seq[pos+3], seq[pos+4], seq[pos+5], seq[pos+6]);
//
//		bool found = false;
//		found = found || ((seq[pos+1]=='T' || seq[pos+1]=='t') && (seq[pos+2]=='A' || seq[pos+2]=='a') && (seq[pos+3]=='A' || seq[pos+3]=='a'));
//		found = found || ((seq[pos+1]=='T' || seq[pos+1]=='t') && (seq[pos+2]=='G' || seq[pos+2]=='g') && (seq[pos+3]=='A' || seq[pos+3]=='a'));
//		found = found || ((seq[pos+1]=='T' || seq[pos+1]=='t') && (seq[pos+2]=='A' || seq[pos+2]=='a') && (seq[pos+3]=='G' || seq[pos+3]=='g'));
//		return found;
//	}
//	else
//	{
//		//printf("TTS(%c): %c%c%c %c%c%c %c%c%c\n", strand, seq[pos-6], seq[pos-5], seq[pos-4], seq[pos-3], seq[pos-2], seq[pos-1], seq[pos], seq[pos+1], seq[pos+2]);
//		bool found = false;
//		found = found || ((seq[pos-3]=='T' || seq[pos-3]=='t') && (seq[pos-2]=='T' || seq[pos-2]=='t') && (seq[pos-1]=='A' || seq[pos-1]=='a'));
//		found = found || ((seq[pos-3]=='T' || seq[pos-3]=='t') && (seq[pos-2]=='C' || seq[pos-2]=='c') && (seq[pos-1]=='A' || seq[pos-1]=='a'));
//		found = found || ((seq[pos-3]=='C' || seq[pos-3]=='c') && (seq[pos-2]=='T' || seq[pos-2]=='t') && (seq[pos-1]=='A' || seq[pos-1]=='a'));
//		return found;
//
//	}
//	return 0;
//}
//
//void Region::load_genomic_sequence()
//{
//	if (!check_region())
//	{
//		printf("load_genomic_sequence: check_region failed\n");
//		print(stderr);
//		exit(-1);
//	}
//	seq = gio->read_flat_file(chr_num, start-1, stop-1);
//}

void Region::print(_IO_FILE*& fd)
{
	fprintf(fd, "region %s\n", get_region_str());
	fprintf(fd, "region start:\t%i\n", start);
	fprintf(fd, "region stop:\t%i\n", stop);
	fprintf(fd, "region strand:\t%c\n", strand);
	if (chr_num>0)
		fprintf(fd, "region chr_num:\t%i\n", chr_num);
	//if (gio && chr_num>0)
	//	fprintf(fd, "region chr:\t%s\n", gio->get_contig_name(chr_num));
	else if (chr)
		fprintf(fd, "region chr:\t%s\n", chr);
}

char* Region::get_region_str()
{
	char* reg_str = new char[1000];
	if (chr)
		sprintf(reg_str, "%s:%i-%i", chr, start, stop);
	//else if (gio && chr_num>0)
	//	sprintf(reg_str, "%s:%i-%i", gio->get_contig_name(chr_num), start, stop);
	else
	{
		fprintf(stderr, "genome information object not set");
		exit(-1);
	}
	return reg_str;
}


#endif
