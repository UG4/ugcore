#include "../include/ng_info.h"

#include <stdlib.h>

struct ng_info* ng_info_new(void)
{
    struct ng_info* ninfo = malloc(sizeof(struct ng_info));
    
    ninfo->error = 0;
    ninfo->err_msg = nullptr;
    
    return ninfo;
}

void ng_info_delete(struct ng_info* ninfo)
{
    if(ninfo->err_msg)
        free((char*)ninfo->err_msg);
    
    free(ninfo);
}
