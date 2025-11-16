#include "../include/lgm.h"

#include <stdlib.h>
#include <stdio.h>

#include "lgm_parser.h"


struct lgm* lgm_new(void)
{
    /* allocate */
    struct lgm* l = malloc(sizeof(struct lgm));

    /* initialize */
    l->name = nullptr;
    l->problemname = nullptr;
    l->convex = 0;

    l->dim = 3;

    l->num_subdomains = 0;
    l->subdomains = nullptr;

    l->num_lines = 0;
    l->lines = nullptr;

    l->num_surfaces = 0;
    l->surfaces = nullptr;

    l->num_points = 0;
    l->points = nullptr;

    /* return */
    return l;
}

void lgm_delete(struct lgm* l)
{
    int i;

    /* clean up domain info */
    free((char*)l->name);
    free((char*)l->problemname);

    /* clean up subdomains */
    for(i = 0; i < l->num_subdomains; ++i)
    {
        free((char*)l->subdomains[i]);
    }

    /* clean up lines */
    for(i = 0; i < l->num_lines; ++i)
    {
        free(l->lines[i].points);
    }

    /* clean up surfaces */
    for(i = 0; i < l->num_surfaces; ++i)
    {
        free(l->surfaces[i].points);
        free(l->surfaces[i].lines);
        free(l->surfaces[i].triangles);
    }

    /* clean up children */
    free(l->subdomains);
    free(l->lines);
    free(l->surfaces);

    /* clean up points */
    if(l->points)
    {
    	free(l->points[0]);
    	free(l->points);
    }

    /* clean up lgm */
    free(l);
}


int lgm_read(const char* filename, struct lgm* l, struct lgm_info* fileinfo)
{
    if(fileinfo)
    {
        fileinfo->error = 0;
        fileinfo->err_msg = nullptr;
    }

    /* parse lgm file */
    return lgm_parse(filename, l, fileinfo);
}
