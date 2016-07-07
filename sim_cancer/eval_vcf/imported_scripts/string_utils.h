#ifndef __STRING_UTILS_H__
#define __STRING_UTILS_H__
#include <cstring>
#include <stdio.h>
#include <stdarg.h>

char* paste( const char* Format, ... )
{
	char* str = new char[10000]; 
	char* pstr = str; 
	va_list Arguments;
	va_start(Arguments, Format);
	double FArg;
	int IArg;
	int off = 0; 
	for (int i = 0; Format[i] != '\0'; ++i )
	{
		if (Format[i] == 'f')
		{
			FArg=va_arg(Arguments, double);
			off = sprintf(pstr, "%.3f",FArg);
		}
		else if (Format[i] == 'i')
		{
			IArg=va_arg(Arguments, int);
			off = sprintf(pstr, "%d",IArg);
		}
		else if (Format[i] == 's')
		{
			char* SArg=va_arg(Arguments, char*);
			off = sprintf(pstr, "%s", SArg); 
		}
		else
		{
			printf("[%s] Error: unknown format: %c, (%s, %i)\n", __func__, Format[i], Format, i); 
			exit(0);
		}
		pstr += off;
	}
	va_end(Arguments);


	return str; 
}

void substitute(char* str1, char* str2, const char* p1, const char* p2)
{
	// find pattern and replace in str2
	sprintf(str2, "%s", str1); 
	char* p = strstr(str2, p1); 
	assert(p); 
	int len = sprintf(p, "%s", p2); 

	// append the rest of str1
	char* pp = strstr(str1, p1); 
	assert(pp); 
	sprintf(p+len, "%s", pp+strlen(p1)); 
}

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

std::vector<char*> my_str_tok(char* line, const char* sep)
{
	std::vector<char*> ret;
	char* tok = strtok(line, sep); 
	while (tok)
	{
		ret.push_back(tok); 
		tok = strtok(NULL, sep);
	}
	return ret;
}
#endif
