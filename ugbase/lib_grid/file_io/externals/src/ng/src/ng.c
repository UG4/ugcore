#include "../include/ng.h"

#include <stdlib.h>
#include <stdio.h>

#include "ng_parser.h"


struct ng* ng_new(void)
{
    /* allocate */
    struct ng* n = malloc(sizeof(struct ng));
    
    /* initialize */
    n->num_bnodes = 0;
    n->bnodes = NULL;
    
    n->num_inodes = 0;
    n->inodes = NULL;
    
    n->num_elements = 0;
    n->elements = NULL;
    
	n->dim = 3;
	
    /* return */
    return n;
}

void ng_delete(struct ng* n)
{
    int i, j;
    
    /* clean up bnodes */
    for(i = 0; i < n->num_bnodes; ++i)
    {
        free(n->bnodes[i].spos);
        free(n->bnodes[i].lpos);
    }
    
    /* clean up elements */
    for(i = 0; i < n->num_elements; ++i)
    {
        for(j = 0; j < n->elements[i].num_faces; ++j)
        {
            free(n->elements[i].faces[j].nodes);
        }
        
        free(n->elements[i].nodes);
        free(n->elements[i].faces);
    }
    
    /* clean up children */
    free(n->bnodes);
    free(n->inodes);
    free(n->elements);
    
    /* clean up ng */
    free(n);
}


int ng_read(const char* filename, struct ng* n, struct ng_info* fileinfo)
{
    /* initialize info if given */
    if(fileinfo)
    {
        fileinfo->error = 0;
        fileinfo->err_msg = NULL;
    }
    
    /* parse ng file */
    return ng_parse(filename, n, fileinfo);
}
