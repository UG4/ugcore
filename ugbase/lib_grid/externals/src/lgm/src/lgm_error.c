#include "lgm_error.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int lgm_error(struct lgm_info* info, const char* msg)
{
    char err_buf[LGM_ERRLEN];

    /* report error */
    if(info)
    {
        info->error = 1;
        sprintf(err_buf, "%s", msg);
        if(info->err_msg)
            free((char*)info->err_msg);
        info->err_msg = malloc(strlen(err_buf)+1);
        strcpy((char*)info->err_msg, err_buf);
    }

    /* return error */
    return 1;
}

int lgm_error_string(struct lgm_info* info, const char* msg, const char* str)
{
    char err_buf[LGM_ERRLEN];

    /* report error */
    if(info)
    {
        info->error = 1;
        sprintf(err_buf, msg, str);
        if(info->err_msg)
            free((char*)info->err_msg);
        info->err_msg = malloc(strlen(err_buf)+1);
        strcpy((char*)info->err_msg, err_buf);
    }

    /* return error */
    return 1;
}

int lgm_error_parse(struct lgm_info* info, const char* msg, tokstream* ts)
{
    char err_buf[LGM_ERRLEN];

    /* report error */
    if(info)
    {
        info->error = 1;
        sprintf(err_buf, msg, ts_line(ts), ts_char(ts));
        if(info->err_msg)
            free((char*)info->err_msg);
        info->err_msg = malloc(strlen(err_buf)+1);
        strcpy((char*)info->err_msg, err_buf);
    }

    /* return error */
    return 1;
}
