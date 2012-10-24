#ifndef NG_ERROR_H_
#define NG_ERROR_H_

#include "../include/ng_info.h"

#include "../../../include/tokstream/tokstream.h"


/* length of temporary error message buffer */
#define NG_ERRLEN 512


int ng_error(struct ng_info* info, const char* msg);
int ng_error_string(struct ng_info* info, const char* msg, const char* str);
int ng_error_parse(struct ng_info* info, const char* msg, tokstream* ts);


#endif /*NG_ERROR_H_*/
