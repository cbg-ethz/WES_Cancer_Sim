
#ifndef __UTILS_H__
#define __UTILS_H__
#include <string.h>
#include <assert.h>
#include <string>
	using std::string;


bool check_triplet(char* DNA, const char* triplet);
char triplet_to_AA(char* DNA);


bool check_triplet(char* DNA, const char* triplet)
{
	assert(strlen(DNA)>=3);
	assert(strlen(triplet)>=3);

	bool match=true;
	for (int i=0; i<3; i++)
		match = match && toupper(DNA[i])==toupper(triplet[i]);
	
	return match;
}

int IUB_code(char a, char*& code)
{
	char x = toupper(a);

	if (x == 'A' || x == 'C' || x == 'T' || x == 'G')
	{
		sprintf(code, "%c", a); 
	}
	else if ( x == 'R' )
	{
		sprintf(code, "AG"); 
	}
	else if ( x == 'Y' )
	{
		sprintf(code, "CT"); 
	}
	else if ( x == 'K' )
	{
		sprintf(code, "GT"); 
	}
	else if ( x == 'M' )
	{
		sprintf(code, "AC"); 
	}
	else if ( x == 'S' )
	{
		sprintf(code, "GC"); 
	}
	else if ( x == 'W' )
	{
		sprintf(code, "AT"); 
	}
	else if ( x == 'B' )
	{
		sprintf(code, "CGT"); 
	}
	else if ( x == 'D' )
	{
		sprintf(code, "AGT"); 
	}
	else if ( x == 'H' )
	{
		sprintf(code, "ACT"); 
	}
	else if ( x == 'V' )
	{
		sprintf(code, "ACG"); 
	}
	else if ( x == 'N' )
	{
		sprintf(code, "ACGT"); 
	}
	else
	{
		(code)[0] = '\0'; 
	}
	return strlen(code); 
}

bool IUB_compare(char c1, char c2)
{
	char* code1 = new char(5); 
	char* code2 = new char(5); 

	int len1 = IUB_code(c1, code1); 
	int len2 = IUB_code(c2, code2); 

	for (int i=0; i<len1; i++)
	{
		for (int j=0; j<len2; j++)
		{
			if (code1[i]==code2[j])
			{
				delete code1; 
				delete code2; 
				return true; 
			}
		}
	}
	delete code1; 
	delete code2; 
	return false; 
}

char triplet_to_AA(char* DNA)
{
	if (check_triplet(DNA, "TCA")) return 'S';
	if (check_triplet(DNA, "TCC")) return 'S';
	if (check_triplet(DNA, "TCG")) return 'S';
	if (check_triplet(DNA, "TCT")) return 'S';
	if (check_triplet(DNA, "TTC")) return 'F';
	if (check_triplet(DNA, "TTT")) return 'F';
	if (check_triplet(DNA, "TTA")) return 'L';
	if (check_triplet(DNA, "TTG")) return 'L';
	if (check_triplet(DNA, "TAC")) return 'Y';
	if (check_triplet(DNA, "TAT")) return 'Y';
	if (check_triplet(DNA, "TAA")) return '_';
	if (check_triplet(DNA, "TAG")) return '_';
	if (check_triplet(DNA, "TGC")) return 'C';
	if (check_triplet(DNA, "TGT")) return 'C';
	if (check_triplet(DNA, "TGA")) return '_';
	if (check_triplet(DNA, "TGG")) return 'W';
	if (check_triplet(DNA, "CTA")) return 'L';
	if (check_triplet(DNA, "CTC")) return 'L';
	if (check_triplet(DNA, "CTG")) return 'L';
	if (check_triplet(DNA, "CTT")) return 'L';
	if (check_triplet(DNA, "CCA")) return 'P';
	if (check_triplet(DNA, "CCC")) return 'P';
	if (check_triplet(DNA, "CCG")) return 'P';
	if (check_triplet(DNA, "CCT")) return 'P';
	if (check_triplet(DNA, "CAC")) return 'H';
	if (check_triplet(DNA, "CAT")) return 'H';
	if (check_triplet(DNA, "CAA")) return 'Q';
	if (check_triplet(DNA, "CAG")) return 'Q';
	if (check_triplet(DNA, "CGA")) return 'R';
	if (check_triplet(DNA, "CGC")) return 'R';
	if (check_triplet(DNA, "CGG")) return 'R';
	if (check_triplet(DNA, "CGT")) return 'R';
	if (check_triplet(DNA, "ATA")) return 'I';
	if (check_triplet(DNA, "ATC")) return 'I';
	if (check_triplet(DNA, "ATT")) return 'I';
	if (check_triplet(DNA, "ATG")) return 'M';
	if (check_triplet(DNA, "ACA")) return 'T';
	if (check_triplet(DNA, "ACC")) return 'T';
	if (check_triplet(DNA, "ACG")) return 'T';
	if (check_triplet(DNA, "ACT")) return 'T';
	if (check_triplet(DNA, "AAC")) return 'N';
	if (check_triplet(DNA, "AAT")) return 'N';
	if (check_triplet(DNA, "AAA")) return 'K';
	if (check_triplet(DNA, "AAG")) return 'K';
	if (check_triplet(DNA, "AGC")) return 'S';
	if (check_triplet(DNA, "AGT")) return 'S';
	if (check_triplet(DNA, "AGA")) return 'R';
	if (check_triplet(DNA, "AGG")) return 'R';
	if (check_triplet(DNA, "GTA")) return 'V';
	if (check_triplet(DNA, "GTC")) return 'V';
	if (check_triplet(DNA, "GTG")) return 'V';
	if (check_triplet(DNA, "GTT")) return 'V';
	if (check_triplet(DNA, "GCA")) return 'A';
	if (check_triplet(DNA, "GCC")) return 'A';
	if (check_triplet(DNA, "GCG")) return 'A';
	if (check_triplet(DNA, "GCT")) return 'A';
	if (check_triplet(DNA, "GAC")) return 'D';
	if (check_triplet(DNA, "GAT")) return 'D';
	if (check_triplet(DNA, "GAA")) return 'E';
	if (check_triplet(DNA, "GAG")) return 'E';
	if (check_triplet(DNA, "GGA")) return 'G';
	if (check_triplet(DNA, "GGC")) return 'G';
	if (check_triplet(DNA, "GGG")) return 'G';
	if (check_triplet(DNA, "GGT")) return 'G';

	return '.';
}
#endif
