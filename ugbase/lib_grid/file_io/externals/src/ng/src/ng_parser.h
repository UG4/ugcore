#ifndef NG_PARSER_H_
#define NG_PARSER_H_

#include "../include/ng.h"
#include "../include/ng_info.h"

#define NG_SEP " \t\r\n"
#define NG_DELIM "#;"


int ng_parse(const char* file, struct ng* n, struct ng_info* fileinfo);


#endif /*NG_PARSER_H_*/
