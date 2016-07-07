
#ifndef __BCF_UTILS_H__
#define __BCF_UTILS_H__

#include <string>
	using std::string;
#include <string.h>
#include <stdio.h>
#include <vector>
	using std::vector;
#include <map>
#include <unordered_set>
#include <unordered_map>
	using std::map;
	using std::unordered_multimap;
//#include "gtf_tools.h"
#include <algorithm>
#include "region.h"
#include <limits> 
#include <zlib.h>
#include <math.h>
#include <errno.h>

#include <ctype.h>
#include "variant.h"
#include <assert.h>


// read bcf file
#include "bcftools/bcf.h"
//#include "bcf.h"  // 2016-03-29 
#include "kstring.h"
#include "khash.h"
#include "kseq.h"
#include "bcf.h"
//#include "id_mapping.h"

KSTREAM_INIT(gzFile, gzread, 4096)

typedef struct {
    int d[4];
} anno16_t;

//static int test16(bcf1_t *b, anno16_t *a);
int test16(bcf1_t *b, anno16_t *a);
float get_AF(bcf1_t *b);
bool is_indel(bcf1_t *b);
int get_DP4(bcf1_t *b, anno16_t *a);
int read_bcf(char* fn_bcf, vector<variant>* variants);
int read_bcf(char* fn_bcf, vector<variant>* variants, float cutoff, float freq, int conf_cnt, const char* normal_tag);
int read_bcf_case_and_control(vector<char*>* control_bcf, vector<char*>* cases_bcf,  vector<vector<variant>* >* all_var, float bcf_qual_cutoff, float conf_frac);


//static int test16(bcf1_t *b, anno16_t *a)
int test16(bcf1_t *b, anno16_t *a)
{
    char *p;
    int i, anno[16];
    a->d[0] = a->d[1] = a->d[2] = a->d[3] = 0.;
    if ((p = strstr(b->info, "I16=")) == 0) return -1;
    p += 4;
    for (i = 0; i < 16; ++i) {
        errno = 0; anno[i] = strtol(p, &p, 10);
        if (anno[i] == 0 && (errno == EINVAL || errno == ERANGE)) return -2;
        ++p;
    }
    memcpy(a->d, anno, 4 * sizeof(int));
    return 1;
}

char gt_int2char(int y, int sample)
{
	if (y>>7&1)
	{
		return '.';
	}
	else if (sample==0)
	{
		return '0' + (y>>3&7); 
	}
	else 
	{
		return '0' + (y&7); 
	}
}


int get_DP4(bcf1_t *b, anno16_t *a)
{
	char *p;
	if ((p = strstr(b->info, "DP4=")) == 0) return -1;
	p += 4;
	for (int i=0; i<4; i++)
	{
		a->d[i] = strtol(p, &p, 10);
		if (i<3)
			assert(*p==',');
		p++; 
	}
	return 1; 
}

static void *locate_field(const bcf1_t *b, const char *fmt, int l)
{
	int i;
	uint32_t tmp;
	tmp = bcf_str2int(fmt, l);
	for (i = 0; i < b->n_gi; ++i)
		if (b->gi[i].fmt == tmp) break;
	return i == b->n_gi? 0 : b->gi[i].data;
}

bool is_indel(bcf1_t *b)
{
	//char *p;
	//if ((p = strstr(b->info, "INDEL")) == 0) return false;
	return bcf_is_indel(b); 
}

bool get_from_info(const char* tag, char* info, float* res)
{
	if (!info)
	{
		printf("[%s] Warning: info fields is NULL (tag: %s)\n", __func__, tag); 
		return false; 
	}

	char *p;
	if ((p = strstr(info, tag)) == 0) 
	{
		return false;
	}
    p += strlen(tag);
    *res = strtod(p, NULL);
	if ((errno != 0)) 
	{
		printf("Error parsing double from string: tag:%s, p:%s info:%s, %f\n", tag, p, info, *res); 
		printf("[%s] Error: %s\n", __func__, strerror(errno)); 
		errno=0; 
		char* p2 = strstr(p, "E"); 
		char* p3 = strstr(p, "e"); 
		if ((p2 && atoi(++p2)<-300) || (p3 && atoi(++p3)<-300))
		{
			printf("[%s] val too small, return 0.0 exp: %s %s\n", __func__, p2, p3); 
			*res = 0.0; 
			return true; 
		}
		return false;
	}
    return true;
}

bool get_from_info(const char* tag, char* info, int* res)
{
	char *p;
	if ((p = strstr(info, tag)) == 0) 
	{
		return false;
	}
    p += strlen(tag);
    *res = strtol(p, NULL, 0);
	if ((errno != 0)) 
	{
		printf("[%s] Error: %s\n", __func__, strerror(errno)); 
		return -2.0;
	}
    return true;
}

float get_AF(bcf1_t *b)
{
    char *p;
    if ((p = strstr(b->info, "AF1=")) == 0) 
	{
		return -1.0;
	}
    p += 4;
    float ret = strtof(p, NULL);
	if ((errno != 0)) 
	{
		printf("Error: %s\n", strerror(errno)); 
		return -2.0;
	}
    return ret;
}

typedef struct {
	gzFile fp;
	FILE *fpout;
	kstream_t *ks;
	void *refhash;
	kstring_t line;
	int max_ref;
} vcf_t;

int my_bcf_sync(bcf1_t *b)
{   
    char *p, *tmp[5];
    int n = b->n_smpl;
    // set ref, alt, flt, info, fmt
    b->ref = b->alt = b->flt = b->info = b->fmt = 0;
    for (p = b->str, n = 0; p < b->str + b->l_str; ++p) {
        if (*p == 0 && p+1 != b->str + b->l_str) {
            if (n == 5) {
                ++n;
                break;
            } else tmp[n++] = p + 1;
        }
    }
    if (n != 5) {
		if (n==4)
		{
    		b->ref = tmp[0]; 
			b->alt = tmp[1]; 
			b->flt = tmp[2]; 
			b->info = tmp[3];
		}
		else
		{
       		fprintf(stderr, "[%s] incorrect number of fields (%d != 5) at %d:%d\n", __func__, n, b->tid, b->pos);
		}
        return -1;
    }
    return 0;
}


int my_vcf_read(bcf_t *bp, bcf_hdr_t *h, bcf1_t *b)
{

	bool sync_done=false; 
	int dret, k, i, sync = 0;
	vcf_t *v = (vcf_t*)bp->v;
	char *p, *q;
	kstring_t str, rn;
	ks_tokaux_t aux, a2;
	//if (!bp->is_vcf) return bcf_read(bp, h, b);
	v->line.l = 0;
	str.l = 0; str.m = b->m_str; str.s = b->str;
	rn.l = rn.m = h->l_nm; rn.s = h->name;
	if (ks_getuntil(v->ks, '\n', &v->line, &dret) < 0) return -1;
	b->n_smpl = h->n_smpl;
	for (p = kstrtok(v->line.s, "\t", &aux), k = 0; p; p = kstrtok(0, 0, &aux), ++k) {
		*(char*)aux.p = 0;
		if (k == 0) { // ref
			int tid = bcf_str2id(v->refhash, p);
			if (tid < 0) {
				tid = bcf_str2id_add(v->refhash, strdup(p));
				kputs(p, &rn); kputc('\0', &rn);
				sync = 1;
			}
			b->tid = tid;
		} else if (k == 1) { // pos
			b->pos = atoi(p) - 1;
		} else if (k == 5) { // qual
			b->qual = (p[0] >= '0' && p[0] <= '9')? atof(p) : 0;
		} else if (k <= 8) { // variable length strings
			kputs(p, &str); kputc('\0', &str);
			b->l_str = str.l; b->m_str = str.m; b->str = str.s;
			if (k == 8)
			{
				bcf_sync(b);
				sync_done=true; 
			}
		} else { // k > 9
			if (strncmp(p, "./.", 3) == 0) {
				for (i = 0; i < b->n_gi; ++i) {
					if (b->gi[i].fmt == bcf_str2int("GT", 2)) {
						((uint8_t*)b->gi[i].data)[k-9] = 1<<7;
					} else if (b->gi[i].fmt == bcf_str2int("GQ", 2)) {
						((uint8_t*)b->gi[i].data)[k-9] = 0;
					} else if (b->gi[i].fmt == bcf_str2int("AD", 2)) 
					{
						// new field // alloc mem here instead of in bcf_sync
						b->gi[i].len = 2;
						b->gi[i].data = realloc(b->gi[i].data, b->n_smpl*b->gi[i].len);
						((uint8_t*)b->gi[i].data)[(k-9)*2] = 0;// k-9 is the current sample with k==9 for the first sample
						((uint8_t*)b->gi[i].data)[(k-9)*2+1] = 0;
					}
					else if (b->gi[i].fmt == bcf_str2int("SP", 2)) {
						((int32_t*)b->gi[i].data)[k-9] = 0;
					} else if (b->gi[i].fmt == bcf_str2int("DP", 2) || b->gi[i].fmt == bcf_str2int("DV", 2)) {
						((uint16_t*)b->gi[i].data)[k-9] = 0;
					} else if (b->gi[i].fmt == bcf_str2int("PL", 2)) {
						int y = b->n_alleles * (b->n_alleles + 1) / 2;
						memset((uint8_t*)b->gi[i].data + (k - 9) * y, 0, y);
					} else if (b->gi[i].fmt == bcf_str2int("GL", 2)) {
						int y = b->n_alleles * (b->n_alleles + 1) / 2;
						memset((float*)b->gi[i].data + (k - 9) * y, 0, y * 4);
					}
				}
				goto endblock;
			}
			for (q = kstrtok(p, ":", &a2), i = 0; q && i < b->n_gi; q = kstrtok(0, 0, &a2), ++i) {
				if (b->gi[i].fmt == bcf_str2int("GT", 2)) {
					((uint8_t*)b->gi[i].data)[k-9] = (q[0] - '0')<<3 | (q[2] - '0') | (q[1] == '/'? 0 : 1) << 6;
				} else if (b->gi[i].fmt == bcf_str2int("GQ", 2)) {
					double _x = strtod(q, &q);
					int x = (int)(_x + .499);
					if (x > 255) x = 255;
					((uint8_t*)b->gi[i].data)[k-9] = x;
				} else if (b->gi[i].fmt == bcf_str2int("AD", 2)) 
				{
					b->gi[i].len = 2;
					b->gi[i].data = realloc(b->gi[i].data, b->n_smpl*b->gi[i].len);

					int x = strtol(q, &q, 10);
					if (x > 255) x = 255;
					assert(*q==','); 
					q++; 
					int x2 = strtol(q, &q, 10);
					if (x2 > 255) x2 = 255;

					((uint8_t*)b->gi[i].data)[(k-9)*2] = x;
					((uint8_t*)b->gi[i].data)[(k-9)*2+1] = x2;
					//printf("p:%s q:%s x:%i, x2:%i\n", p, q, x, x2); 
				}
				else if (b->gi[i].fmt == bcf_str2int("SSC", 3)) 
				{
					b->gi[i].len = 1;
					b->gi[i].data = realloc(b->gi[i].data, b->n_smpl*b->gi[i].len);

					int x = strtol(q, &q, 10);
					if (x > 255) x = 255;
					//assert(*q==','); 
					//q++; 
					//int x2 = strtol(q, &q, 10);
					//if (x2 > 255) x2 = 255;

					((uint8_t*)b->gi[i].data)[(k-9)] = x;
					//printf("p:%s q:%s x:%i, x2:%i\n", p, q, x, x2); 
				}
				else if (b->gi[i].fmt == bcf_str2int("DP4", 3)) 
				{
					b->gi[i].len = 4;
					b->gi[i].data = realloc(b->gi[i].data, b->n_smpl*b->gi[i].len);

					for (int xx=0; xx<4; xx++)
					{
						int x = strtol(q, &q, 10);
						if (x > 255) x = 255;
						if (xx<3)
							assert(*q==','); 
						q++; 
						((uint8_t*)b->gi[i].data)[(k-9)*b->gi[i].len+xx] = x;
					}
					//printf("p:%s q:%s x:%i, x2:%i\n", p, q, x, x2); 
				}
				else if (b->gi[i].fmt == bcf_str2int("SP", 2)) {
					int x = strtol(q, &q, 10);
					if (x > 0xffff) x = 0xffff;
					((uint32_t*)b->gi[i].data)[k-9] = x;
				} else if (b->gi[i].fmt == bcf_str2int("DP", 2) || b->gi[i].fmt == bcf_str2int("DV", 2)) {
					int x = strtol(q, &q, 10);
					if (x > 0xffff) x = 0xffff;
					((uint16_t*)b->gi[i].data)[k-9] = x;
				} else if (b->gi[i].fmt == bcf_str2int("PL", 2)) {
					int x, y, j;
					uint8_t *data = (uint8_t*)b->gi[i].data;
					y = b->n_alleles * (b->n_alleles + 1) / 2;
					for (j = 0; j < y; ++j) {
						x = strtol(q, &q, 10);
						if (x > 255) x = 255;
						data[(k-9) * y + j] = x;
						++q;
					}
				} else if (b->gi[i].fmt == bcf_str2int("GL", 2)) {
					int j, y;
					float x, *data = (float*)b->gi[i].data;
					y = b->n_alleles * (b->n_alleles + 1) / 2;
					for (j = 0; j < y; ++j) {
						x = strtod(q, &q);
						data[(k-9) * y + j] = x > 0? -x/10. : x;
						++q;
					}
				}
			}
		endblock: i = i;
		}
	}
	if (!sync_done)
		my_bcf_sync(b); 
	h->l_nm = rn.l; h->name = rn.s;
	if (sync) bcf_hdr_sync(h);
	return v->line.l + 1;
}



int read_bcf(char* fn_bcf, vector<variant>* variants, float cutoff, float freq, int conf_cnt , const char* normal_tag )
{
	bcf_t *bp;
	bcf_hdr_t *h;

	if (strstr(fn_bcf, ".vcf"))
	{
		printf("[%s] open as vcf file: %s\n", __func__, fn_bcf); 
		bp = vcf_open(fn_bcf, "r"); 
	}
	else
	{
		printf("[%s] open as bcf file: %s\n", __func__, fn_bcf); 
		bp = bcf_open(fn_bcf, "r"); 
	}

	if (bp == 0) 
	{
		if (strstr(fn_bcf, ".vcf"))
			fprintf(stderr, "[%s] fail to open the VCF file.\n", __func__);
		else
			fprintf(stderr, "[%s] fail to open the BCF file.\n", __func__);
		return -1;
	}

	bcf1_t* b = new bcf1_t();
	if (bp->is_vcf)
		h = vcf_hdr_read(bp);
	else
		h = bcf_hdr_read(bp);
	int ret; 
	int cnt = 0;
	int discard=0; 
	int max_mq0 = 4; 
	int discard_mq0=0; 
	int discard_dp=0; 
	int discard_qual=0; 
	int discard_low_qd=0; 
	int discard_fischer=0; 
	int discard_conf_cnt = 0; 
	while (true)
	//while (cnt<100)
	{
		if (!bp->is_vcf)
		{
			if ((ret = bcf_read(bp, h, b)) <= 0)
				break; 
		}
		else
		{
			if ((ret = my_vcf_read(bp, h, b)) <= 0)
				break; 
		}


		char* chr = h->ns[b->tid];
		cnt++;



		int mq0; 
		if (get_from_info("MQ0=", b->info, &mq0))
		{
			/*
			// from http://seqanswers.com/wiki/How-to/exome_analysis
			//
			//java -Xmx4g -jar GenomeAnalysisTK.jar \
			//-R hg19.fa \
			//-T VariantFiltration \
			//-B:variant,VCF snp.vcf.recalibrated \
			//-o snp.recalibrated.filtered.vcf \
			//--clusterWindowSize 10 \
			//--filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
			//--filterName "HARD_TO_VALIDATE" \
			//--filterExpression "DP < 5 " \
			//--filterName "LowCoverage" \
			//--filterExpression "QUAL < 30.0 " \
			//--filterName "VeryLowQual" \
			//--filterExpression "QUAL > 30.0 && QUAL < 50.0 " \
			//--filterName "LowQual" \
			//--filterExpression "QD < 1.5 " \
			//--filterName "LowQD" \
			//--filterExpression "SB > -10.0 " \
			//--filterName "StrandBias"
			*/

			//gatk stuff
			int dp; 
			assert(get_from_info("DP=", b->info, &dp)); 
			if (false)
			{
				if (mq0>max_mq0 && ((float)mq0)/dp>0.1 )
				{
					discard_mq0++; 
					continue; 
				}
				if (dp<5)
				{
					discard_dp++; 
					continue; 
				}
				if ( b->qual < 40.0)
				{
					discard_qual++; 
					continue; 
				}
				float qd; 
				assert(get_from_info("QD=", b->info, &qd));
				if (qd<1.5)
				{
					discard_low_qd++; 
					continue;
				}
				float fs; 
				assert(get_from_info("FS=", b->info, &fs));
				if (fs>150.0)
				{
					discard_fischer++; 
					continue; 
				}
			}
		}

		//if (is_indel(b))
		//	continue; 


		variant v = init_var();
		v.chr = (char*) (new string(chr))->c_str();

		// bcftools subtracts 1 from the pos when reading in vcf files. 
		v.pos = b->pos+1; 
		v.conf_cnt = 0;
		v.non_conf_cnt = 0;
		anno16_t a;
		int ret = get_DP4(b, &a); 


		if (ret>0) 
		{
			// samvar
			v.conf_cnt = a.d[2]+a.d[3];
			v.non_conf_cnt = a.d[0]+a.d[1];
			//printf("conf_cnt: %i non_conf_cnt:%i\n", v.conf_cnt, v.non_conf_cnt); 

			if (b->n_smpl==2 && normal_tag)
			{
				// samtools version with two samples
				uint8_t* ad = (uint8_t*) locate_field(b, "DP4", 3);

				int g_no=0; 
				int g_tu=1; 
				if (!strstr(h->sns[g_no], "muscle"))
				{
					assert(strstr(h->sns[g_no], "NO")); 
					assert(strstr(h->sns[g_tu], "TU")); 
				}

				if (ad)
				{
					//printf("samtools1.2: found field DP4: %i %i %i %i\n", ad[0], ad[1], ad[2], ad[3]); 
					//printf("samtools1.2: found field DP4: %i %i %i %i\n", ad[4*g_tu+0], ad[4*g_tu+1], ad[4*g_tu+2], ad[4*g_tu+3]); 

					//##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
					v.conf_cnt_normal = ad[g_no*4+2] + ad[g_no*4+3]; 
					v.non_conf_cnt_normal = ad[g_no*4+0] + ad[g_no*4+1]; 

					v.conf_cnt = ad[g_tu*4+2] + ad[g_tu*4+3]; 
					v.non_conf_cnt = ad[g_tu*4+0] + ad[g_tu*4+1]; 
					
					//int random = rand(); 
					//if (random%100>99)
					//{
					//	exit(0); 
					//}
				}
				else
				{
					printf("Samtools1.2: did not find field DP4\n"); 
					exit(0); 
				}

			}
		}	
		else if (get_from_info("SRF=", b->info, &v.non_conf_cnt))
		{
			//freebayes
			
			//##INFO=<ID=SRF,Number=1,Type=Integer,Description="Number of reference observations on the forward strand">
			float non_conf_forward; 
			get_from_info("SRF=", b->info, &non_conf_forward); 

			//##INFO=<ID=SRR,Number=1,Type=Integer,Description="Number of reference observations on the reverse strand">
			float non_conf_reverse; 
			get_from_info("SRR=", b->info, &non_conf_reverse); 

			//##INFO=<ID=SAF,Number=A,Type=Integer,Description="Number of alternate observations on the forward strand">
			float conf_forward; 
			get_from_info("SAF=", b->info, &conf_forward); 

			//##INFO=<ID=SAR,Number=A,Type=Integer,Description="Number of alternate observations on the reverse strand">
			float conf_reverse; 
			get_from_info("SAR=", b->info, &conf_reverse); 

			v.conf_cnt = conf_forward + conf_reverse; 
			v.non_conf_cnt = non_conf_forward + non_conf_reverse; 
		}
		else if(get_from_info("TR=", b->info, &v.non_conf_cnt))
		{
			// varscan2
			assert(get_from_info("TV=", b->info, &v.conf_cnt));
			assert(get_from_info("NV=", b->info, &v.conf_cnt_normal)); 
			assert(get_from_info("NR=", b->info, &v.non_conf_cnt_normal)); 

			float somatic_p_val; 
			if (get_from_info("SPV=", b->info, &somatic_p_val))
			{
				//printf("somatic_p_val: %f\n", somatic_p_val); 
				v.qual = -10*log(somatic_p_val+1e-255)/log(10); 
				if (v.qual>250)
					v.qual = 250.0; 
			}
			else 
				v.qual = b->qual; 

		}
		else if (get_from_info("NV=", b->info, &v.conf_cnt_normal))
		{
			//deepSNV: INFO: NV=5;CN=126;TV=1;CT=01;VF=-0.039265579005139;SF=0.000385218292299463;RP=0.0430496558732555
			assert(get_from_info("TV=", b->info, &v.conf_cnt)); 

			int cov_normal; 
			assert(get_from_info("CN=", b->info, &cov_normal)); 
			int cov_tumor; 
			assert(get_from_info("CT=", b->info, &cov_tumor)); 

			v.non_conf_cnt = cov_tumor - v.conf_cnt; 
			v.non_conf_cnt_normal = cov_normal - v.conf_cnt_normal; 
			float pval; 
			assert(get_from_info("RP=", b->info, &pval)); 

			v.qual = -10*log(pval+1e-255)/log(10); 
			
			if (v.qual>250) 
				v.qual=250.0; 
			
			//printf("bcf_fmt: %s\n", bcf_fmt(h, b)); 
			//printf("%i %i %i %i %f \n", v.conf_cnt, v.non_conf_cnt, v.conf_cnt_normal, v.non_conf_cnt_normal, v.qual); 
			//exit(0);  
		}
		else
		{
			// GATK and Mutect: Allel specific depth

			//printf("bcf_fmt: %s\n", bcf_fmt(h, b)); 
			uint8_t* ad = (uint8_t*) locate_field(b, "AD", 2);
			for (int i=0; i<b->n_smpl; i++)
			{
				//printf("retrieve count for sample %i: %i, %i\n", i, ad[i*2], ad[i*2+1]);  
			}
			if (!ad && h->sns)
			{
				// SomaticSniper
				ad = (uint8_t*) locate_field(b, "DP4", 3);

				int g_no=0; 
				int g_tu=1; 
				assert(strstr(h->sns[g_no], "NORMAL")); 
				assert(strstr(h->sns[g_tu], "TUMOR")); 
				if (ad)
				{
					//printf("SomaticSniper: found field DP4: %i %i %i %i\n", ad[0], ad[1], ad[2], ad[3]); 
					//printf("SomaticSniper: found field DP4: %i %i %i %i\n", ad[4*g+0], ad[4*g+1], ad[4*g+2], ad[4*g+3]); 

					//##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
					v.conf_cnt_normal = ad[g_no*4+2] + ad[g_no*4+3]; 
					v.non_conf_cnt_normal = ad[g_no*4+0] + ad[g_no*4+1]; 

					v.conf_cnt = ad[g_tu*4+2] + ad[g_tu*4+3]; 
					v.non_conf_cnt = ad[g_tu*4+0] + ad[g_tu*4+1]; 
				}
				else
				{
					printf("SomaticSniper: did not find field DP4\n"); 
					exit(0); 
				}
				//printf("SomaticSniper: fmt: %s\n", b->fmt); 
				//printf("SomaticSniper: flt: %s\n", b->flt); 
				//printf("SomaticSniper: info: %s\n", b->info); 
				//printf("SomaticSniper: gi[0]: %s\n", (char*) b->gi[0].data); 

				uint8_t* score = (uint8_t*) locate_field(b, "SSC", 3);


				if (score)
				{
					//printf("SomaticSniper: SCC: %i %i\n", score[0], score[1]); 
					v.qual = score[g_tu]; 
				}
			}
			else if (b->n_smpl==2 && normal_tag)
			{
				// check which of the two is the normal	
				if (strstr(h->sns[0], normal_tag) || strstr(h->sns[0], "NORMAL") || strstr(h->sns[0], "normal") || strstr(h->sns[0], "muscle"))
				{
					v.conf_cnt_normal = ad[1]; 
					v.non_conf_cnt_normal = ad[0]; 

					v.conf_cnt = ad[3]; 
					v.non_conf_cnt = ad[2]; 
				}
				else if (strstr(h->sns[1], normal_tag) || strstr(h->sns[1], "NORMAL") || strstr(h->sns[1], "normal") || strstr(h->sns[1], "muscle"))
				{
					v.conf_cnt_normal = ad[3]; 
					v.non_conf_cnt_normal = ad[2]; 

					v.conf_cnt = ad[1]; 
					v.non_conf_cnt = ad[0]; 
				}
				else
				{
					printf("none of the sample names match the normal tag\n"); 
					printf("%s %s, tag:%s\n", h->sns[0], h->sns[1], normal_tag); 
					exit(-1); 
				}
				if (v.pos==80727829)
					printf("variant Map4k3: %i %i %i %i %s %s %p\n", v.non_conf_cnt, v.conf_cnt, v.non_conf_cnt_normal, v.conf_cnt_normal, h->sns[0], h->sns[1], strstr(h->sns[0], normal_tag)); 
					
				v.qual = b->qual; 
				if (v.qual>2000) 
					v.qual=2000;  
			}
		}

		v.sample = -1;
		if (v.qual<0)
		{
			v.qual = b->qual;
			if (v.qual>255) 
			 	v.qual = 255.0; 
		}

		//if (v.qual<=0)
		{
			// this is for Mutect
			// read the filter field of the vcf file
			if (strcmp(b->flt, "REJECT")==0)
			{
				v.qual = - ((float) (rand()%1000))/100; 		
			}
			else if (strcmp(b->flt, "PASS")==0)
			{
				v.qual = 10.0+((float) (rand()%1000))/100; 
			}
		}

		v.source = BCF;
		v.ref = (char*) (new string(b->ref))->c_str();
		v.alt = (char*) (new string(b->alt))->c_str();
		v.freq = -1;
		if (v.conf_cnt>0)
			v.freq = ((float) v.conf_cnt)/(v.conf_cnt + v.non_conf_cnt); 

		v.gt_normal = (char*) (new string("x/x"))->c_str(); 
		v.gt_cancer = (char*) (new string("x/x"))->c_str(); 

		uint8_t* gt = (uint8_t*) locate_field(b, "GT", 2);
		if (gt && b->n_smpl==2)
		{
			v.gt_normal[0] = gt_int2char(gt[0], 0);  
			v.gt_normal[2] = gt_int2char(gt[0], 1);  
			v.gt_cancer[0] = gt_int2char(gt[1], 0);  
			v.gt_cancer[2] = gt_int2char(gt[1], 1);  

		}
		bool take_freq = true || v.conf_cnt<0 || v.non_conf_cnt<0 ; 
		take_freq = (take_freq || (v.non_conf_cnt == 0 && v.conf_cnt>1)); 
		take_freq = (take_freq || (v.non_conf_cnt>0 && ((float) v.conf_cnt)/v.non_conf_cnt>freq)); 
	
		take_freq = v.conf_cnt>conf_cnt; 
		//if (v.qual>cutoff && ! take_freq)
		//	printf("discard: conf %i,  non conf%i\n",  v.conf_cnt, v.non_conf_cnt); 
		discard_conf_cnt += !take_freq; 
		if (cutoff==0.0 && freq==0.0 && conf_cnt==0)
			variants->push_back(v); 
		else if (v.qual>cutoff && take_freq)
			variants->push_back(v); 
		else
		{
			//printf("%.4f <= %.4f, || %i <= %i\n", v.qual, cutoff, v.conf_cnt, conf_cnt); 
			discard++; 
		}
	}
	printf("return %lu variants; discard: %i (%i) \n", variants->size(),  discard, discard_conf_cnt); 
	if (discard_mq0+discard_dp+discard_qual+discard_low_qd+discard_fischer>0)
	{
		printf("GATK filter: discard_mq0>%i:%i discard_dp:%i discard_qual:%i discard_low_qd:%i discard_fischer:%i\n", max_mq0, discard_mq0, discard_dp, discard_qual, discard_low_qd, discard_fischer);
	}
	return 1; 
}

int read_bcf(char* fn_bcf, vector<variant>* variants)
{
	return read_bcf(fn_bcf, variants, 0.0, 0.0, 0, NULL); 
}


int read_bcf_case_and_control(vector<char*>* control_bcf, vector<char*>* cases_bcf,  vector<vector<variant>* >* all_var, float bcf_qual_cutoff, float conf_frac)
{
	vector<variant> normal;
	for (unsigned int i=0; i<control_bcf->size(); i++)
	{
		read_bcf(control_bcf->at(i), &normal, 0.0, 0.0, 0, NULL);	
	}
	int orig_len = normal.size(); 
	for (int i=0; i<orig_len; i++)
	{
		int len1 = strlen(normal[i].ref); 
		int len2 = strlen(normal[i].alt); 
		if (len1>1 || len2>1)
		{
			//for (int j=normal[i].pos-1; j<normal[i].pos+len1 || j<normal[i].pos+len2; j++)
			//{
			//	variant v = init_var();
			//	v.chr = normal[i].chr;
			//	v.pos = j; 
			//	normal.push_back(v); 	
			//}
		}
	}
	printf("sorting %i variants in normal (+ %lu dummys)\n", orig_len, normal.size()-orig_len); 
	sort(normal.begin(), normal.end(), variant_compare); 	

	for (unsigned int i=0; i<cases_bcf->size(); i++)
	{
		printf("read variants from bcf file %s\n", cases_bcf->at(i));
		vector<variant> var;
		read_bcf(cases_bcf->at(i), &var, bcf_qual_cutoff, conf_frac, 5, NULL); 

		sort(var.begin(), var.end(), variant_compare); 	

		vector<variant>* aquired = new vector<variant>(var.size()+normal.size());
		vector<variant>::iterator it;
		it = set_difference(var.begin(), var.end(), normal.begin(), normal.end(), aquired->begin(), variant_compare); 
		aquired->resize(it-aquired->begin());
		all_var->push_back(aquired); 
	}
	return 0; 
}


#endif
