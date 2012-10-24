#ifndef LGM_ERROR_H_
#define LGM_ERROR_H_

#include "../include/lgm_info.h"

#include "../../../include/tokstream/tokstream.h"


/* length of temporary error message buffer */
#define LGM_ERRLEN 512


int lgm_error(struct lgm_info* info, const char* msg);
int lgm_error_string(struct lgm_info* info, const char* msg, const char* str);
int lgm_error_parse(struct lgm_info* info, const char* msg, tokstream* ts);


#endif /*LGM_ERROR_H_*/
