#include "../include/lgm.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../include/lgm_line.h"
#include "../include/lgm_surface.h"
#include "../include/lgm_info.h"
#include "lgm_error.h"


/****
 * section headers
 */

static const char* header_dom = "Domain-Info";
static const char* header_sub = "Unit-Info";
static const char* header_lin = "Line-Info";
static const char* header_sfc = "Surface-Info";
static const char* header_pnt = "Point-Info";


/****
 * inlining macros
 */

/* try-like construct, if cmd returns 1 (error), return 1 too */
#define $(cmd) do { if(cmd) return 1; } while (0)


/****
 * error messages
 */

static const char* err_msg_fil = "Error opening file %s.";


/****
 * lgm file writing
 */

int lgm_write_domain(FILE* fp, const struct lgm* l /*, struct lgm_info* fileinfo*/);
int lgm_write_subdomains(FILE* fp, const struct lgm* l /*, struct lgm_info* fileinfo*/);
int lgm_write_lines(FILE* fp, const struct lgm* l /*, struct lgm_info* fileinfo*/);
int lgm_write_surfaces(FILE* fp, const struct lgm* l /*, struct lgm_info* fileinfo*/);
int lgm_write_points(FILE* fp, const struct lgm* l /*, struct lgm_info* fileinfo*/);


int lgm_write(const char* file, const struct lgm* l, struct lgm_info* fileinfo)
{
    /* open file */
    FILE* fp = fopen(file, "w");
    if(!fp)
        return lgm_error_string(fileinfo, err_msg_fil, file);

    /* write lgm */
    $(lgm_write_domain(fp, l));
    $(lgm_write_subdomains(fp, l));
    $(lgm_write_lines(fp, l));
    $(lgm_write_surfaces(fp, l));
    $(lgm_write_points(fp, l));

    /* close file */
    fclose(fp);

    /* success */
    return 0;
}

int lgm_write_domain(FILE* fp, const struct lgm* l /*, struct lgm_info* fileinfo*/)
{
    char c = l->convex ? '1' : '0';

    /* write header */
    fputs("# ", fp);
    fputs(header_dom, fp);
    fputc('\n', fp);

    /* write name */
    fputs("name = ", fp);
    fputs(l->name, fp);
    fputc('\n', fp);

    /* write problemname */
    fputs("problemname = ", fp);
    fputs(l->problemname, fp);
    fputc('\n', fp);

    /* write convex */
    fputs("convex = ", fp);
    fputc(c, fp);
    fputc('\n', fp);

    /* write newline */
    fputc('\n', fp);

    /* success */
    return 0;
}

int lgm_write_subdomains(FILE* fp, const struct lgm* l /*, struct lgm_info* fileinfo*/)
{
    int i;

    /* write header */
    fputs("# ", fp);
    fputs(header_sub, fp);
    fputc('\n', fp);

    /* write subdomains, skip global domain (i = 1) */
    for(i = 1; i < l->num_subdomains; ++i)
        fprintf(fp, "unit %d %s\n", i, l->subdomains[i]);

    /* write newline */
    fputc('\n', fp);

    /* success */
    return 0;
}

int lgm_write_lines(FILE* fp, const struct lgm* l /*, struct lgm_info* fileinfo*/)
{
    int i, j;

    /* write header */
    fputs("# ", fp);
    fputs(header_lin, fp);
    fputc('\n', fp);

    /* write lines */
    for(i = 0; i < l->num_lines; ++i)
    {
        fprintf(fp, "line %d: ", i);

        /* write points */
        fputs("points:", fp);
        for(j = 0; j < l->lines[i].num_points; ++j)
            fprintf(fp, " %d", l->lines[i].points[j]);
        fputc(';', fp);

        fputc('\n', fp);
    }

    /* write newline */
    fputc('\n', fp);

    /* success */
    return 0;
}

int lgm_write_surfaces(FILE* fp, const struct lgm* l /*, struct lgm_info* fileinfo*/)
{
    int i, j;

    /* write header */
    fputs("# ", fp);
    fputs(header_sfc, fp);
    fputc('\n', fp);

    /* write surfaces */
    for(i = 0; i < l->num_surfaces; ++i)
    {
        fprintf(fp, "surface %d: ", i);
        fprintf(fp, "left=%d; ", l->surfaces[i].left);
        fprintf(fp, "right=%d; ", l->surfaces[i].right);

        /* write points */
        if(l->surfaces[i].num_points)
        {
            fputs("points:", fp);
            for(j = 0; j < l->surfaces[i].num_points; ++j)
                fprintf(fp, " %d", l->surfaces[i].points[j]);
            fputs("; ", fp);
        }

        /* write lines */
        if(l->surfaces[i].num_lines)
        {
            fputs("lines:", fp);
            for(j = 0; j < l->surfaces[i].num_lines; ++j)
                fprintf(fp, " %d", l->surfaces[i].lines[j]);
            fputs("; ", fp);
        }

        /* write surfaces */
        if(l->surfaces[i].num_triangles)
        {
            fputs("triangles:", fp);
            for(j = 0; j < l->surfaces[i].num_triangles; ++j)
                fprintf(fp, " %d %d %d;", l->surfaces[i].triangles[j][0], l->surfaces[i].triangles[j][1], l->surfaces[i].triangles[j][2]);
        }

        fputc('\n', fp);
    }

    /* write newline */
    fputc('\n', fp);

    /* success */
    return 0;
}

int lgm_write_points(FILE* fp, const struct lgm* l /*, struct lgm_info* fileinfo*/)
{
    int i, j;

    /* write header */
    fputs("# ", fp);
    fputs(header_pnt, fp);
    fputc('\n', fp);

    /* write points */
    for(i = 0; i < l->num_points; ++i)
    {
        for(j = 0; j < l->dim; ++j)
        {
            if(j > 0)
                fputc(' ', fp);
            fprintf(fp, "%.15G", l->points[i][j]);
        }
        fputc('\n', fp);
    }

    /* success */
    return 0;
}
