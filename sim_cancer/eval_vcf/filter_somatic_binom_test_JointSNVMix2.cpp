#include <stdio.h>
#include <string.h>
#include <string>
#include <cstdlib>
#include <iostream>
#include <iomanip>
using std::cout;
using std::cerr;
#include "stats_tests_binom.h"
#include <fstream>


std::vector<char*> separate(char* str, char sep)
{
	std::vector<char*> ret;
	int i=0;
	ret.push_back(str);
	while (str[i]!=0)
	{
		if (str[i]==sep)
		{
			str[i]=0;
			ret.push_back(str+i+1);
		}
		i++;
	}
	return ret;
}

std::vector<char*> get_fields(char* line)
{
	std::vector<char*> ret;
	int i=0;
	ret.push_back(line);
	while (line[i]!=0)
	{
		if (line[i]=='\t')
		{
			line[i]=0;
			ret.push_back(line+i+1);
		}
		i++;
	}
	return ret;
}


int main(int argc, char *args[]){

        if (argc < 2) {
                printf("Usage: %s <JointSNVMix2_Raw> <JointSNVMix2_Filtered>\n", args[0]);
                return 0;
        }
	
	std::ifstream JSM2_input (args[1]);
	std::ofstream JSM2_output;
	JSM2_output.open(args[2]);
	std::string line;

	float var_rate = 0.005; // this is the sequencing error rate we assume

	int discardCnt=0;
	int totalCnt=0;
	int includeCnt=0;
	if(JSM2_input.is_open()){
		while ( getline (JSM2_input,line)){ // loop over all lines

			char lineAsCharArray[10000];
			strncpy(lineAsCharArray, line.c_str(), sizeof(lineAsCharArray));
			if(strncmp("#",lineAsCharArray,1) == 0 ){ // if first character is '#' then it is a header line and is just printed out
				JSM2_output << lineAsCharArray << "\n";
				continue;
			}

			totalCnt++;
			int NR=-10;
			int NV=-10;
			std::vector<char*> lineAsCharVector = get_fields(lineAsCharArray); // split up the line according to tabs
			for(int i=0;i<lineAsCharVector.size();i++){ // loop over columns in that line
				
				if(i==7){
					std::vector<char*>  eighthColAsVector = separate(lineAsCharVector[i],';');
					for(int j=0; j<eighthColAsVector.size();j++){
						//cout << eighthColAsVector[j] << "\n";
						if(j==0){
							std::vector<char*> NRfield = separate(eighthColAsVector[j],'=');
							NR=atoi(NRfield[1]);
						} else if (j==1){
							std::vector<char*> NVfield = separate(eighthColAsVector[j],'=');
							NV=atoi(NVfield[1]);
						}
					}
				}
			}
			assert(NR!=-10 && NV!=-10);
			int Ncov=NV+NR; // coverage in normal sample is number of reference reads + number of variant reads
			float pval = -1.0;
			if(Ncov == 0 || NV == 0){
				pval=1;
			} else if(Ncov == NV){
				pval=0;
			} else {
				pval = binom_cdf(NV, Ncov, var_rate,"greater"); // perform binomial test; null-hypothesis: only seq error; alt-hypothesis: germline mut
			}
			if(pval >= 0.05){ // if pval is less than 0.05, reject the null-hypoth. and classify as germline -> filter out
				JSM2_output << line << "\n";
				includeCnt++;
			} else {
				discardCnt++; // discard it, because probably germline
			}
		}
	} else {
		cerr << "Unable to find open file: " << args[1] << "\n";
	}

	assert(discardCnt+includeCnt==totalCnt);
	cout << "We found in total " << totalCnt << " variants, of which "<< discardCnt << " where discarded because they are likely germline variants. And "<< includeCnt << " where kept becuase nothing spoke against seeing them as somatic.\n" ;
			
	/*
	estimated illumina error rate ~0.1% -> rate = 0.001
	here we use 0.005 to be more on the conservative side, i.e. the varCount needs to be higher in order to have a rejection = germline classification = filter out variant
	NullHypothesis: varCount in normal is just sequencing error
	AltHypothesis: varCount in normal is higher and cannot be due to seq error -> must be germline
	in R: p.val=binom.test(varCount, coverage, p = 0.005, alternative = "greater",conf.level = 0.95)$p.val
	*/
	
	JSM2_input.close();
	JSM2_output.close();
	return 0;
}


