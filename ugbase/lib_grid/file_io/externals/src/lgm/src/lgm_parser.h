#ifndef LGM_PARSER_H_
#define LGM_PARSER_H_

#include "../include/lgm.h"
#include "../include/lgm_info.h"

/* separators which are ignored as whitespace */
#define LGM_SEP " \t\r\n"

/* delimiters are read as tokens */
#define LGM_DELIM "#=:;"

/* name of global domain */
#define LGM_DOMAIN_NAME "DOMAIN"


int lgm_parse(const char* file, struct lgm* l, struct lgm_info* fileinfo);
int lgm_parse_2d(const char* file, struct lgm* l, struct lgm_info* fileinfo);

#endif