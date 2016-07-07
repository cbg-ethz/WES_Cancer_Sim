
#ifndef _BASIC_TOOLS_H__
#define _BASIC_TOOLS_H__

#include <vector>

std::vector<char*> separate(char* str, char sep);
std::vector<char*> get_fields(char* line);
std::vector<char*> my_str_tok(char* line, const char* sep);

#endif
