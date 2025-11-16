/* TODO: Use carray for dynamic arrays. */

#include "lgm_parser.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../../../include/tokstream/tokstream.h"

#include "../include/lgm.h"
#include "../include/lgm_line.h"
#include "../include/lgm_surface.h"
#include "../include/lgm_info.h"
#include "lgm_error.h"


/*
 * default sizes for array allocations
 * always counts per entry
 */
static constexpr int init_num_subdomains = 64;
static constexpr int init_num_lines = 128;
static constexpr int init_num_surfaces = 256;
static constexpr int init_num_coords = 512;

static constexpr int init_num_line_points = 32;

static constexpr int init_num_surface_points = 32;
static constexpr int init_num_surface_lines = 32;
static constexpr int init_num_surface_triangles = 32;


/*
 * error messages
 */

static const char* err_msg_tok = "Could not read token at line %d, char %d.";
static const char* err_msg_end = "Could not find end of entry from line %d, char %d.";
static const char* err_msg_mem = "Failed to allocate memory for LGM data.";
static const char* err_msg_fil = "Error opening file \"%s\".";
static const char* err_msg_eof = "Could not read all of file, stopped at line %d, char %d.";

static const char* err_msg_neq = "Expected equal sign at line %d, char %d.";
static const char* err_msg_ndc = "Expected double colon at line %d, char %d.";
static const char* err_msg_nsc = "Expected semicolon at line %d, char %d.";
static const char* err_msg_idm = "Index and id do not match at line %d, char %d.";

static const char* err_msg_dom = "Could not read domain info at line %d, char %d.";
static const char* err_msg_sub = "Could not read subdomain info at line %d, char %d.";
static const char* err_msg_lin = "Could not read line info at line %d, char %d.";
static const char* err_msg_sur = "Could not read surface info at line %d, char %d.";
static const char* err_msg_pnt = "Could not read point info at line %d, char %d.";

static const char* err_msg_cvx = "Convex setting must be '0' or '1' at line %d, char %d.";

static const char* err_msg_nou = "Not a subdomain [internal error], line %d, char %d.";

static const char* err_msg_nol = "Not a line [internal error], line %d, char %d.";
static const char* err_msg_lid = "Could not read line id at line %d, char %d.";
static const char* err_msg_nlp = "Not line points [internal error], line %d, char %d.";

static const char* err_msg_nos = "Not a surface [internal error], line %d, char %d.";
static const char* err_msg_sid = "Could not read surface id at line %d, char %d.";
static const char* err_msg_sle = "Could not read surface left at line %d, char %d.";
static const char* err_msg_sri = "Could not read surface right at line %d, char %d.";
static const char* err_msg_nsp = "Not surface points [internal error], line %d, char %d.";
static const char* err_msg_nsl = "Not surface lines [internal error], line %d, char %d.";
static const char* err_msg_nst = "Not surface triangles [internal error], line %d, char %d.";
static const char* err_msg_tri = "Could not read triangle index at line %d, char %d.";


/*
 * inlining macros
 */

/* try-like construct, if cmd returns 1 (error), return 1 too */
#define $(cmd) do { if(cmd) return 1; } while (0)

/* short for str != nullptr && strcmp == 0 */
#define is_token(str, tok) (str != nullptr && strcmp(str, tok) == 0)


/*
 * internal functions
 */

int lgm_parser_strtoi(const char* str, int* i);
int lgm_parser_strtod(const char* str, double* d);


/****
 * lgm file parsing
 */

/* parsing functions */
int lgm_parse_domain_info(tokstream* ts, struct lgm* l, struct lgm_info* fileinfo);

int lgm_parse_subdomain_info(tokstream* ts, struct lgm* l, struct lgm_info* fileinfo);
int lgm_parse_subdomain(tokstream* ts, struct lgm* l, int* alloc_subdomains, struct lgm_info* fileinfo);

int lgm_parse_line_info(tokstream* ts, struct lgm* l, struct lgm_info* fileinfo);
int lgm_parse_line(tokstream* ts, int id, struct lgm_line* line, struct lgm_info* fileinfo);
int lgm_parse_line_points(tokstream* ts, struct lgm_line* line, struct lgm_info* fileinfo);

int lgm_parse_surface_info(tokstream* ts, struct lgm* l, struct lgm_info* fileinfo);
int lgm_parse_surface(tokstream* ts, int id, struct lgm_surface* surface, struct lgm_info* fileinfo);
int lgm_parse_surface_points(tokstream* ts, struct lgm_surface* surface, struct lgm_info* fileinfo);
int lgm_parse_surface_lines(tokstream* ts, struct lgm_surface* surface, struct lgm_info* fileinfo);
int lgm_parse_surface_triangles(tokstream* ts, struct lgm_surface* surface, struct lgm_info* fileinfo);

int lgm_parse_point_info(tokstream* ts, struct lgm* l, struct lgm_info* fileinfo);


int lgm_parse(const char* file, struct lgm* l, struct lgm_info* fileinfo)
{
    /* open token stream */
    tokstream* ts = ts_open(file);
    if(!ts)
        return lgm_error_string(fileinfo, err_msg_fil, file);

    /* set separators */
    ts_sep(ts, LGM_SEP);

    /* set delimiters */
    ts_delim(ts, LGM_DELIM);

    /* set buffer size */
    ts_bufsiz(ts, 8192);

    /* read file */
    $(lgm_parse_domain_info(ts, l, fileinfo));
    $(lgm_parse_subdomain_info(ts, l, fileinfo));
    $(lgm_parse_line_info(ts, l, fileinfo));
	if(l->dim == 3)
    	$(lgm_parse_surface_info(ts, l, fileinfo));
    $(lgm_parse_point_info(ts, l, fileinfo));

    /* make sure whole file was read */
    if(!ts_eof(ts))
        return lgm_error_parse(fileinfo, err_msg_eof, ts);

    /* close token stream */
    ts_close(ts);

    /* success */
    return 0;
}

int lgm_parse_domain_info(tokstream* ts, struct lgm* l, struct lgm_info* fileinfo)
{
    const char* tok;

    /* domain header */
    tok = ts_get(ts);
    if(!is_token(tok, "#"))
        return lgm_error_parse(fileinfo, err_msg_dom, ts);
    tok = ts_get(ts);
    if(!is_token(tok, "Domain-Info"))
        return lgm_error_parse(fileinfo, err_msg_dom, ts);

    /* initialize domain properties */
    l->name = nullptr;
    l->problemname = nullptr;
    l->convex = 0;

    /* read domain properties */
    while((tok = ts_get(ts)))
    {
        /* check for name token */
        if(is_token(tok, "name"))
        {
            /* read equal sign */
            tok = ts_get(ts);
            if(!is_token(tok, "="))
                return lgm_error_parse(fileinfo, err_msg_neq, ts);

            /* read name property */
            tok = ts_getline(ts);
            l->name = strcpy(malloc(strlen(tok)+1), tok);

            continue;
        }

        /* check for problemname token */
        if(is_token(tok, "problemname"))
        {
            /* read equal sign */
            tok = ts_get(ts);
            if(!is_token(tok, "="))
                return lgm_error_parse(fileinfo, err_msg_neq, ts);

            /* read problemname property */
            tok = ts_getline(ts);
            l->problemname = strcpy(malloc(strlen(tok)+1), tok);

            continue;
        }

        /* check for convex token */
        if(is_token(tok, "convex"))
        {
            /* read equal sign */
            tok = ts_get(ts);
            if(!is_token(tok, "="))
                return lgm_error_parse(fileinfo, err_msg_neq, ts);

            /* read convex property */
            /* tok = ts_getline(ts); */ /*replaced with the next line. Caused problems with windows line-end and with spaces behind the number. */
            tok = ts_get(ts);
            if(lgm_parser_strtoi(tok, &l->convex))
                return lgm_error_parse(fileinfo, err_msg_cvx, ts);

            /* read the rest of the line */
            ts_getline(ts);
            continue;
        }

        /* end reading properties */
        $(ts_unget(ts));
        break;
    }

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return lgm_error_parse(fileinfo, err_msg_tok, ts);

    /* success */;
    return 0;
}

int lgm_parse_subdomain_info(tokstream* ts, struct lgm* l, struct lgm_info* fileinfo)
{
    const char* tok;

    int alloc_subdomains = init_num_subdomains;

    /* initialize subdomains */
    l->num_subdomains = 0;
    l->subdomains = malloc(sizeof(const char*) * alloc_subdomains);
    if(!l->subdomains)
        return lgm_error(fileinfo, err_msg_mem);

    /* subdomain header */
    tok = ts_get(ts);
    if(!is_token(tok, "#"))
        return lgm_error_parse(fileinfo, err_msg_sub, ts);
    tok = ts_get(ts);
    if(!is_token(tok, "Unit-Info"))
        return lgm_error_parse(fileinfo, err_msg_sub, ts);

    /* insert global domain */
    l->subdomains[0] = strcpy(malloc(strlen(LGM_DOMAIN_NAME)+1), LGM_DOMAIN_NAME);
    ++l->num_subdomains;

    /* scan subdomains */
    while((tok = ts_get(ts)))
    {
        /* check for subdomain */
        if(is_token(tok, "unit"))
        {
            /* read subdomain, also does reallocation */
            $(lgm_parse_subdomain(ts, l, &alloc_subdomains, fileinfo));

            continue;
        }

        /* end subdomains */
        $(ts_unget(ts));
        break;
    }

    /* finalize subdomains */
    l->subdomains = realloc(l->subdomains, sizeof(const char*) * l->num_subdomains);

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return lgm_error_parse(fileinfo, err_msg_tok, ts);

    /* success */
    return 0;
}

int lgm_parse_subdomain(tokstream* ts, struct lgm* l, int* alloc_subdomains, struct lgm_info* fileinfo)
{
    const char* tok;
    int subdomain_start;
    int id;

    /* check that this is a subdomain */
    tok = ts_tok(ts);
    if(!is_token(tok, "unit"))
        return lgm_error_parse(fileinfo, err_msg_nou, ts);

    /* store start of subdomains which will all get the same name */
    subdomain_start = l->num_subdomains;

    /* scan ids */
    while((tok = ts_get(ts)))
    {
        /* check for id */
        if(lgm_parser_strtoi(tok, &id) == 0)
        {
            /* check for correct id */
            if(id != l->num_subdomains)
                return lgm_error_parse(fileinfo, err_msg_idm, ts);

            /* check for reallocation */
            if(l->num_subdomains == (*alloc_subdomains))
            {
                (*alloc_subdomains) *= 2;
                l->subdomains = realloc(l->subdomains, sizeof(const char*) * (*alloc_subdomains));
                if(!l->subdomains)
                    return lgm_error(fileinfo, err_msg_mem);
            }

            /* store no name, so far */
            l->subdomains[l->num_subdomains] = nullptr;
            ++l->num_subdomains;

            continue;
        }

        /* end ids */
        $(ts_unget(ts));
        break;
    }

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return lgm_error_parse(fileinfo, err_msg_tok, ts);

    /* read name */
    tok = ts_getline(ts);

    /* assign name to all subdomains read */
    while(subdomain_start < l->num_subdomains)
        l->subdomains[subdomain_start++] = strcpy(malloc(strlen(tok)+1), tok);

    /* success */
    return 0;
}

int lgm_parse_line_info(tokstream* ts, struct lgm* l, struct lgm_info* fileinfo)
{
    const char* tok;

    int alloc_lines = init_num_lines;

    /* initialize lines */
    l->num_lines = 0;
    l->lines = malloc(sizeof(struct lgm_line) * alloc_lines);
    if(!l->lines)
        return lgm_error(fileinfo, err_msg_mem);

    /* lines header */
    tok = ts_get(ts);
    if(!is_token(tok, "#"))
        return lgm_error_parse(fileinfo, err_msg_lin, ts);
    tok = ts_get(ts);
    if(!is_token(tok, "Line-Info"))
        return lgm_error_parse(fileinfo, err_msg_lin, ts);

    /* scan lines */
    while((tok = ts_get(ts)))
    {
        /* check for line */
        if(is_token(tok, "line"))
        {
            /* check for reallocation */
            if(l->num_lines == alloc_lines)
            {
                alloc_lines *= 2;
                l->lines = realloc(l->lines, sizeof(struct lgm_line) * alloc_lines);
                if(!l->lines)
                    return lgm_error(fileinfo, err_msg_mem);
            }

            /* read line */
            $(lgm_parse_line(ts, l->num_lines, &l->lines[l->num_lines], fileinfo));
            ++l->num_lines;

            continue;
        }

        /* end lines */
        $(ts_unget(ts));
        break;
    }

    /* finalize lines */
    l->lines = realloc(l->lines, sizeof(struct lgm_line) * l->num_lines);

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return lgm_error_parse(fileinfo, err_msg_tok, ts);

    /* success */
    return 0;
}

int lgm_parse_line(tokstream* ts, int id, struct lgm_line* line, struct lgm_info* fileinfo)
{
    const char* tok;
    int line_id;
    int left, right;

    /* check that this is a line */
    tok = ts_tok(ts);
    if(!is_token(tok, "line"))
        return lgm_error_parse(fileinfo, err_msg_nol, ts);

    /* read line id */
    tok = ts_get(ts);
    if(lgm_parser_strtoi(tok, &line_id))
        return lgm_error_parse(fileinfo, err_msg_lid, ts);
    if(line_id != id)
        return lgm_error_parse(fileinfo, err_msg_idm, ts);

    /* read double colon */
    tok = ts_get(ts);
    if(!is_token(tok, ":"))
        return lgm_error_parse(fileinfo, err_msg_ndc, ts);

    /* read entries */
    while((tok = ts_get(ts)))
    {
    	/* check for left / right*/
    	if(is_token(tok, "left"))
    	{
    		tok = ts_get(ts);
            if(!is_token(tok, "="))
                return lgm_error_parse(fileinfo, err_msg_neq, ts);

            /* read left property */
            tok = ts_get(ts);
            if(lgm_parser_strtoi(tok, &left))
                return lgm_error_parse(fileinfo, err_msg_sle, ts);
            //todo: add line->left
            //line->left = left;

            /* read delimiter */
            tok = ts_get(ts);
            if(!is_token(tok, ";"))
                return lgm_error_parse(fileinfo, err_msg_nsc, ts);

            continue;
    	}
    	else if(is_token(tok, "right"))
    	{
    		/* read equal sign */
            tok = ts_get(ts);
            if(!is_token(tok, "="))
                return lgm_error_parse(fileinfo, err_msg_neq, ts);

            /* read left property */
            tok = ts_get(ts);
            if(lgm_parser_strtoi(tok, &right))
                return lgm_error_parse(fileinfo, err_msg_sri, ts);
            //todo: add line->right
            //line->right = right;

            /* read delimiter */
            tok = ts_get(ts);
            if(!is_token(tok, ";"))
                return lgm_error_parse(fileinfo, err_msg_nsc, ts);

            continue;
    	}
        /* check for points */
    	else if(is_token(tok, "points"))
        {
            /* read points */
            $(lgm_parse_line_points(ts, line, fileinfo));

            continue;
        }

        /* end line */
        $(ts_unget(ts));
        break;
    }

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return lgm_error_parse(fileinfo, err_msg_tok, ts);

    /* success */
    return 0;
}

int lgm_parse_line_points(tokstream* ts, struct lgm_line* line, struct lgm_info* fileinfo)
{
    const char* tok;
    int point;

    int alloc_line_points = init_num_line_points;

    /* check that this is line points */
    tok = ts_tok(ts);
    if(!is_token(tok, "points"))
        return lgm_error_parse(fileinfo, err_msg_nlp, ts);

    /* read double colon */
    tok = ts_get(ts);
    if(!is_token(tok, ":"))
        return lgm_error_parse(fileinfo, err_msg_ndc, ts);

    /* initialize points */
    line->num_points = 0;
    line->points = malloc(sizeof(int) * alloc_line_points);
    if(!line->points)
        return lgm_error(fileinfo, err_msg_mem);

    /* scan points */
    while((tok = ts_get(ts)))
    {
        /* check for point */
        if(!lgm_parser_strtoi(tok, &point))
        {
            /* check for reallocation */
            if(line->num_points == alloc_line_points)
            {
                alloc_line_points *= 2;
                line->points = realloc(line->points, sizeof(int) * alloc_line_points);
                if(!line->points)
                    return lgm_error(fileinfo, err_msg_mem);
            }

            /* store point */
            line->points[line->num_points] = point;
            ++line->num_points;

            continue;
        }

        /* end points */
        $(ts_unget(ts));
        break;
    }

    /* finalize points */
    line->points = realloc(line->points, sizeof(int) * line->num_points);

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return lgm_error_parse(fileinfo, err_msg_tok, ts);

    /* read delimiter */
    tok = ts_get(ts);
    if(!is_token(tok, ";")){
        //return lgm_error_parse(fileinfo, err_msg_end, ts);
        $(ts_unget(ts));
    }

    /* success */
    return 0;
}

int lgm_parse_surface_info(tokstream* ts, struct lgm* l, struct lgm_info* fileinfo)
{
    const char* tok;

    int alloc_surfaces = init_num_surfaces;

    /* initialize surfaces */
    l->num_surfaces = 0;
    l->surfaces = malloc(sizeof(struct lgm_surface) * alloc_surfaces);
    if(!l->surfaces)
        return lgm_error(fileinfo, err_msg_mem);

    /* surfaces header */
    tok = ts_get(ts);
    if(!is_token(tok, "#"))
        return lgm_error_parse(fileinfo, err_msg_sur, ts);
    tok = ts_get(ts);
    if(!is_token(tok, "Surface-Info"))
        return lgm_error_parse(fileinfo, err_msg_sur, ts);

    /* scan surfaces */
    while((tok = ts_get(ts)))
    {
        /* check for surface */
        if(is_token(tok, "surface"))
        {
            /* check for reallocation */
            if(l->num_surfaces == alloc_surfaces)
            {
                alloc_surfaces *= 2;
                l->surfaces = realloc(l->surfaces, sizeof(struct lgm_surface) * alloc_surfaces);
                if(!l->surfaces)
                    return lgm_error(fileinfo, err_msg_mem);
            }

            /* read surface */
            $(lgm_parse_surface(ts, l->num_surfaces, &l->surfaces[l->num_surfaces], fileinfo));
            ++l->num_surfaces;

            continue;
        }

        /* end surfaces */
        $(ts_unget(ts));
        break;
    }

    /* finalize surfaces */
    l->surfaces = realloc(l->surfaces, sizeof(struct lgm_surface) * l->num_surfaces);

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return lgm_error_parse(fileinfo, err_msg_tok, ts);

    /* success */
    return 0;
}

int lgm_parse_surface(tokstream* ts, int id, struct lgm_surface* surface, struct lgm_info* fileinfo)
{
    const char* tok;
    int surface_id;
    int left, right;

    /* check that this is a surface */
    tok = ts_tok(ts);
    if(!is_token(tok, "surface"))
        return lgm_error_parse(fileinfo, err_msg_nos, ts);

    /* read surface id */
    tok = ts_get(ts);
    if(lgm_parser_strtoi(tok, &surface_id))
        return lgm_error_parse(fileinfo, err_msg_sid, ts);
    if(surface_id != id)
        return lgm_error_parse(fileinfo, err_msg_idm, ts);

    /* read double colon */
    tok = ts_get(ts);
    if(!is_token(tok, ":"))
        return lgm_error_parse(fileinfo, err_msg_ndc, ts);

    /* read entries */
    while((tok = ts_get(ts)))
    {
        /* check for left */
        if(is_token(tok, "left"))
        {
            /* read equal sign */
            tok = ts_get(ts);
            if(!is_token(tok, "="))
                return lgm_error_parse(fileinfo, err_msg_neq, ts);

            /* read left property */
            tok = ts_get(ts);
            if(lgm_parser_strtoi(tok, &left))
                return lgm_error_parse(fileinfo, err_msg_sle, ts);
            surface->left = left;

            /* read delimiter */
            tok = ts_get(ts);
            if(!is_token(tok, ";"))
                return lgm_error_parse(fileinfo, err_msg_nsc, ts);

            continue;

        }

        /* check for right */
        if(is_token(tok, "right"))
        {
            /* read equal sign */
            tok = ts_get(ts);
            if(!is_token(tok, "="))
                return lgm_error_parse(fileinfo, err_msg_neq, ts);

            /* read left property */
            tok = ts_get(ts);
            if(lgm_parser_strtoi(tok, &right))
                return lgm_error_parse(fileinfo, err_msg_sri, ts);
            surface->right = right;

            /* read delimiter */
            tok = ts_get(ts);
            if(!is_token(tok, ";"))
                return lgm_error_parse(fileinfo, err_msg_nsc, ts);

            continue;

        }

        /* check for points */
        if(is_token(tok, "points"))
        {
            /* read points */
            $(lgm_parse_surface_points(ts, surface, fileinfo));

            continue;
        }

        /* check for lines */
        if(is_token(tok, "lines"))
        {
            /* read lines */
            $(lgm_parse_surface_lines(ts, surface, fileinfo));

            continue;
        }

        /* check for triangles */
        if(is_token(tok, "triangles"))
        {
            /* read triangles */
            $(lgm_parse_surface_triangles(ts, surface, fileinfo));

            continue;
        }

        /* end line */
        $(ts_unget(ts));
        break;
    }

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return lgm_error_parse(fileinfo, err_msg_tok, ts);

    /* success */
    return 0;
}

int lgm_parse_surface_points(tokstream* ts, struct lgm_surface* surface, struct lgm_info* fileinfo)
{
    const char* tok;
    int point;

    int alloc_surface_points = init_num_surface_points;

    /* check that this is surface points */
    tok = ts_tok(ts);
    if(!is_token(tok, "points"))
        return lgm_error_parse(fileinfo, err_msg_nsp, ts);

    /* read double colon */
    tok = ts_get(ts);
    if(!is_token(tok, ":"))
        return lgm_error_parse(fileinfo, err_msg_ndc, ts);

    /* initialize points */
    surface->num_points = 0;
    surface->points = malloc(sizeof(int) * alloc_surface_points);
    if(!surface->points)
        return lgm_error(fileinfo, err_msg_mem);

    /* scan points */
    while((tok = ts_get(ts)))
    {
        /* check for point */
        if(!lgm_parser_strtoi(tok, &point))
        {
            /* check for reallocation */
            if(surface->num_points == alloc_surface_points)
            {
                alloc_surface_points *= 2;
                surface->points = realloc(surface->points, sizeof(int) * alloc_surface_points);
                if(!surface->points)
                    return lgm_error(fileinfo, err_msg_mem);
            }

            /* store point */
            surface->points[surface->num_points] = point;
            ++surface->num_points;

            continue;
        }

        /* end points */
        $(ts_unget(ts));
        break;
    }

    /* finalize points */
    surface->points = realloc(surface->points, sizeof(int) * surface->num_points);

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return lgm_error_parse(fileinfo, err_msg_tok, ts);

    /* read delimiter */
    tok = ts_get(ts);
    if(!is_token(tok, ";"))
        return lgm_error_parse(fileinfo, err_msg_end, ts);

    /* success */
    return 0;
}

int lgm_parse_surface_lines(tokstream* ts, struct lgm_surface* surface, struct lgm_info* fileinfo)
{
    const char* tok;
    int line;

    int alloc_surface_lines = init_num_surface_lines;

    /* check that this is surface lines */
    tok = ts_tok(ts);
    if(!is_token(tok, "lines"))
        return lgm_error_parse(fileinfo, err_msg_nsl, ts);

    /* read double colon */
    tok = ts_get(ts);
    if(!is_token(tok, ":"))
        return lgm_error_parse(fileinfo, err_msg_ndc, ts);

    /* initialize lines */
    surface->num_lines = 0;
    surface->lines = malloc(sizeof(int) * alloc_surface_lines);
    if(!surface->lines)
        return lgm_error(fileinfo, err_msg_mem);

    /* scan lines */
    while((tok = ts_get(ts)))
    {
        /* check for line */
        if(!lgm_parser_strtoi(tok, &line))
        {
            /* check for reallocation */
            if(surface->num_lines == alloc_surface_lines)
            {
                alloc_surface_lines *= 2;
                surface->lines = realloc(surface->lines, sizeof(int) * alloc_surface_lines);
                if(!surface->lines)
                    return lgm_error(fileinfo, err_msg_mem);
            }

            /* store line */
            surface->lines[surface->num_lines] = line;
            ++surface->num_lines;

            continue;
        }

        /* end lines */
        $(ts_unget(ts));
        break;
    }

    /* finalize lines */
    surface->lines = realloc(surface->lines, sizeof(int) * surface->num_lines);

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return lgm_error_parse(fileinfo, err_msg_tok, ts);

    /* read delimiter */
    tok = ts_get(ts);
    if(!is_token(tok, ";"))
        return lgm_error_parse(fileinfo, err_msg_end, ts);

    /* success */
    return 0;
}

int lgm_parse_surface_triangles(tokstream* ts, struct lgm_surface* surface, struct lgm_info* fileinfo)
{
    const char* tok;
    int triangle;

    int alloc_surface_triangles = init_num_surface_triangles;

    /* check that this is surface triangles */
    tok = ts_tok(ts);
    if(!is_token(tok, "triangles"))
        return lgm_error_parse(fileinfo, err_msg_nst, ts);

    /* read double colon */
    tok = ts_get(ts);
    if(!is_token(tok, ":"))
        return lgm_error_parse(fileinfo, err_msg_ndc, ts);

    /* initialize triangles */
    surface->num_triangles = 0;
    surface->triangles = malloc(sizeof(int[3]) * alloc_surface_triangles);
    if(!surface->triangles)
        return lgm_error(fileinfo, err_msg_mem);

    /* scan triangles */
    while((tok = ts_get(ts)))
    {
        /* check for triangle */
        if(!lgm_parser_strtoi(tok, &triangle))
        {
            /* check for reallocation */
            if(surface->num_triangles == alloc_surface_triangles)
            {
                alloc_surface_triangles *= 2;
                surface->triangles = realloc(surface->triangles, sizeof(int[3]) * alloc_surface_triangles);
                if(!surface->triangles)
                    return lgm_error(fileinfo, err_msg_mem);
            }

            /* store triangle's first index */
            surface->triangles[surface->num_triangles][0] = triangle;

            /* store triangle's second index */
            tok = ts_get(ts);
            if(lgm_parser_strtoi(tok, &triangle))
                return lgm_error_parse(fileinfo, err_msg_tri, ts);
            surface->triangles[surface->num_triangles][1] = triangle;

            /* store triangle's third index */
            tok = ts_get(ts);
            if(lgm_parser_strtoi(tok, &triangle))
                return lgm_error_parse(fileinfo, err_msg_tri, ts);
            surface->triangles[surface->num_triangles][2] = triangle;

            /* read delimiter */
            tok = ts_get(ts);
            if(!is_token(tok, ";"))
                return lgm_error_parse(fileinfo, err_msg_end, ts);

            ++surface->num_triangles;

            continue;
        }

        /* end triangles */
        $(ts_unget(ts));
        break;
    }

    /* finalize triangles */
    surface->triangles = realloc(surface->triangles, sizeof(int[3]) * surface->num_triangles);

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return lgm_error_parse(fileinfo, err_msg_tok, ts);

    /* success */
    return 0;
}

int lgm_parse_point_info(tokstream* ts, struct lgm* l, struct lgm_info* fileinfo)
{
    const char* tok;
    int num_coords;
    double* coords;
    double coord;
    int i;

    int alloc_coords = init_num_coords;

    /* initialize points */
    l->num_points = 0;
    l->points = nullptr;

    /* initialize coordinate array */
    coords = malloc(sizeof(double) * alloc_coords);
    if(!coords)
        return lgm_error(fileinfo, err_msg_mem);

    /* points header */
    tok = ts_get(ts);
    if(!is_token(tok, "#"))
        return lgm_error_parse(fileinfo, err_msg_pnt, ts);
    tok = ts_get(ts);
    if(!is_token(tok, "Point-Info"))
        return lgm_error_parse(fileinfo, err_msg_pnt, ts);

    /* delimiting ';' is optional, ignore it */
    ts_sep_on(ts, ';');

    /* set newline delimiter to get dimension from first point */
    ts_delim_on(ts, '\n');

    /* scan coordinates */
    num_coords = 0;
    while((tok = ts_get(ts)))
    {
    	/* TODO: I have no clue what this is! */
        if(num_coords == 3*46000)
            coord = 2.0f;

        /* check for coordinate */
        if(lgm_parser_strtod(tok, &coord) == 0)
        {
            /* check for reallocation */
            if(num_coords == alloc_coords)
            {
                alloc_coords *= 2;
                coords = realloc(coords, sizeof(double) * alloc_coords);
                if(!coords)
                    return lgm_error(fileinfo, err_msg_mem);
            }

            /* store coordinate */
            coords[num_coords] = coord;
            ++num_coords;

            continue;
        }

        /* check for next point */
        if(is_token(tok, "\n"))
        {
            /* make sure there is a point yet */
            if(num_coords > 0)
            {
                /* first point specifies lgm dimension */
                l->dim = num_coords;

                /* ignore further newlines again */
                ts_delim_off(ts, '\n');
            }

            continue;
        }

        /* end points */
        $(ts_unget(ts));
        break;
    }

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return lgm_error_parse(fileinfo, err_msg_tok, ts);

    /* finalize coordinates */
    coords = realloc(coords, sizeof(double) * num_coords);

    /* set number of points */
    l->num_points = num_coords / l->dim;

    /* set points */
    l->points = malloc(sizeof(double*) * l->num_points);
    for(i = 0; i < l->num_points; ++i)
    	l->points[i] = coords + (l->dim * i);

    /* success */
    return 0;
}


/*
 * internal implementation
 */

int lgm_parser_strtoi(const char* str, int* i)
{
    long int ii;
    char* endptr;

    /* convert */
    ii = strtol(str, &endptr, 10);

    /* check that nothing remains in string */
    if(*endptr)
        return 1;

    /* copy value */
    *i = (int)ii;

    /* success */
    return 0;
}

int lgm_parser_strtod(const char* str, double* d)
{
    double dd;
    char* endptr;

    /* convert */
    dd = strtod(str, &endptr);

    /* check that nothing remains in string */
    if(*endptr)
		return 1;

    /* copy value */
    *d = dd;

    /* success */
    return 0;
}
