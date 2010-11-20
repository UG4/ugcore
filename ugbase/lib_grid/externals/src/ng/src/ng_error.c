#include "ng_error.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int ng_error(struct ng_info* info, const char* msg)
{
    char err_buf[NG_ERRLEN];

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

int ng_error_string(struct ng_info* info, const char* msg, const char* str)
{
    char err_buf[NG_ERRLEN];

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

int ng_error_parse(struct ng_info* info, const char* msg, tokstream* ts)
{
    char err_buf[NG_ERRLEN];

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
