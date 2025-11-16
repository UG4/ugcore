#include "../include/lgm_info.h"

#include <stdlib.h>

struct lgm_info* lgm_info_new(void)
{
    struct lgm_info* linfo = malloc(sizeof(struct lgm_info));
    
    linfo->error = 0;
    linfo->err_msg = nullptr;
    
    return linfo;
}

void lgm_info_delete(struct lgm_info* linfo)
{
    if(linfo->err_msg)
        free((char*)linfo->err_msg);
    
    free(linfo);
}
