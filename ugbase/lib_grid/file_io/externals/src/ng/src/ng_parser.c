/* TODO: Use carray for dynamic arrays. */

#include "ng_parser.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "../../../include/tokstream/tokstream.h"

#include "../include/ng_node.h"
#include "../include/ng_element.h"
#include "ng_error.h"


/*
 * default sizes for array allocations
 * always counts per entry
 */
static constexpr int init_num_bnodes = 128;
static constexpr int init_num_inodes = 128;
static constexpr int init_num_elements = 64;

static constexpr int init_num_spos = 4;
static constexpr int init_num_lpos = 4;

static constexpr int init_num_nodes = 16;
static constexpr int init_num_faces = 8;

/*
 * error messages
 */

static const char* err_msg_tok = "Could not read token at line %d, char %d.";
static const char* err_msg_BIE = "Expected 'B', 'I' or 'E' token at line %d, char %d.";
static const char* err_msg_mem = "Failed to allocate memory for NG data.";
static const char* err_msg_fil = "Error opening file \"%s\".";

static const char* err_msg_nob = "Not a bnode [internal error], line %d, char %d.";
static const char* err_msg_noi = "Not an inode [internal error], line %d, char %d.";
static const char* err_msg_coo = "Error reading coordinates at line %d, char %d.";
static const char* err_msg_SLS = "Expected 'S', 'L' or ';' token at line %d, char %d.";
static const char* err_msg_sem = "Expected ';' token at line %d, char %d.";

static const char* err_msg_nos = "Not a surface_pos [internal error], line %d, char %d.";
static const char* err_msg_sid = "Error reading surface id at line %d, char %d.";
static const char* err_msg_spo = "Error reading surface position at line %d, char %d.";

static const char* err_msg_nol = "Not a line_pos [internal error], line %d, char %d.";
static const char* err_msg_lid = "Error reading line id at line %d, char %d.";
static const char* err_msg_lpo = "Error reading line position at line %d, char %d.";

static const char* err_msg_noe = "Not an element [internal error], line %d, char %d.";
static const char* err_msg_uid = "Error reading subdomain id at line %d, char %d.";
static const char* err_msg_SNF = "Expected node id or 'F' token at line %d, char %d.";

static const char* err_msg_nof = "Not a face nor a side [internal error], line %d, char %d.";


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

int ng_parser_strtoi(const char* str, int* i);
int ng_parser_strtod(const char* str, double* d);


/****
 * ng file parsing
 */

/* parsing functions */
int ng_parse_file(tokstream* ts, struct ng* n, struct ng_info* fileinfo);
int ng_parse_bnode(tokstream* ts, struct ng_bnode* bnode, struct ng* n, struct ng_info* fileinfo);
int ng_parse_inode(tokstream* ts, struct ng_inode* inode, struct ng* n, struct ng_info* fileinfo);
int ng_parse_element(tokstream* ts, struct ng_element* element, struct ng_info* fileinfo, int dim);
int ng_parse_surface_pos(tokstream* ts, struct ng_surface_pos* spos, struct ng_info* fileinfo);
int ng_parse_line_pos(tokstream* ts, struct ng_line_pos* lpos, struct ng_info* fileinfo);
int ng_parse_face(tokstream* ts, struct ng_face* face, struct ng_info* fileinfo);


int ng_parse(const char* file, struct ng* n, struct ng_info* fileinfo)
{
    /* open token stream */
    tokstream* ts = ts_open(file);
    if(!ts)
        return ng_error_string(fileinfo, err_msg_fil, file);

    /* set whitespace separators */
    ts_sep(ts, NG_SEP);

    /* set semicolon delimiter */
    ts_delim(ts, NG_DELIM);

    /* set buffer size */
    ts_bufsiz(ts, 8192);

    /* read file */
    $(ng_parse_file(ts, n, fileinfo));

    /* close token stream */
    ts_close(ts);

    /* success */
    return 0;
}

int ng_parse_file(tokstream* ts, struct ng* n, struct ng_info* fileinfo)
{
    const char* tok;

    int alloc_bnodes = init_num_bnodes;
    int alloc_inodes = init_num_inodes;
    int alloc_elements = init_num_elements;

    /* initialize bnodes */
    n->num_bnodes = 0;
    n->bnodes = malloc(sizeof(struct ng_bnode) * alloc_bnodes);
    if(!n->bnodes)
        return ng_error(fileinfo, err_msg_mem);

    /* scan bnodes */
    while((tok = ts_get(ts)))
    {
        /* skip over comments */
        if(tok[0] == '#')
        {
            ts_skipline(ts);
            continue;
        }

        /* check for bnode */
        if(is_token(tok, "B"))
        {
            /* check for reallocation */
            if(n->num_bnodes == alloc_bnodes)
            {
                alloc_bnodes *= 2;
                n->bnodes = realloc(n->bnodes, sizeof(struct ng_bnode) * alloc_bnodes);
                if(!n->bnodes)
                    return ng_error(fileinfo, err_msg_mem);
            }

            /* read bnode */
            $(ng_parse_bnode(ts, &n->bnodes[n->num_bnodes], n, fileinfo));
            ++n->num_bnodes;

            continue;
        }

        /* end bnodes */
        $(ts_unget(ts));
        break;
    }

    /* finalize bnodes */
    n->bnodes = realloc(n->bnodes, sizeof(struct ng_bnode) * n->num_bnodes);

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return ng_error_parse(fileinfo, err_msg_tok, ts);

    /* initialize inodes */
    n->num_inodes = 0;
    n->inodes = malloc(sizeof(struct ng_inode) * alloc_inodes);
    if(!n->inodes)
        return ng_error(fileinfo, err_msg_mem);

    /* scan inodes */
    while((tok = ts_get(ts)))
    {
        /* skip over comments */
        if(tok[0] == '#')
        {
            ts_skipline(ts);
            continue;
        }

        /* check for inode */
        if(is_token(tok, "I"))
        {
            /* check for reallocation */
            if(n->num_inodes == alloc_inodes)
            {
                alloc_inodes *= 2;
                n->inodes = realloc(n->inodes, sizeof(struct ng_inode) * alloc_inodes);
                if(!n->inodes)
                    return ng_error(fileinfo, err_msg_mem);
            }

            /* read inode */
            $(ng_parse_inode(ts, &n->inodes[n->num_inodes], n, fileinfo));
            ++n->num_inodes;

            continue;
        }

        /* end inodes */
        $(ts_unget(ts));
        break;
    }

    /* finalize inodes */
    n->inodes = realloc(n->inodes, sizeof(struct ng_inode) * n->num_inodes);

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return ng_error_parse(fileinfo, err_msg_tok, ts);

    /* initialize elements */
    n->num_elements = 0;
    n->elements = malloc(sizeof(struct ng_element) * alloc_elements);
    if(!n->elements)
        return ng_error(fileinfo, err_msg_mem);

    /* scan elements */
    while((tok = ts_get(ts)))
    {
        /* skip over comments */
        if(tok[0] == '#')
        {
            ts_skipline(ts);
            continue;
        }

        /* check for element */
        if(is_token(tok, "E"))
        {
            /* check for reallocation */
            if(n->num_elements == alloc_elements)
            {
                alloc_elements *= 2;
                n->elements = realloc(n->elements, sizeof(struct ng_element) * alloc_elements);
                if(!n->elements)
                    return ng_error(fileinfo, err_msg_mem);
            }

            /* read element */
            $(ng_parse_element(ts, &n->elements[n->num_elements], fileinfo, n->dim));
            ++n->num_elements;

            continue;
        }

        /* end elements */
        $(ts_unget(ts));
        break;
    }

    /* finalize elements */
    n->elements = realloc(n->elements, sizeof(struct ng_element) * n->num_elements);

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return ng_error_parse(fileinfo, err_msg_tok, ts);

    /* make sure whole file was read */
    if(!ts_eof(ts))
        return ng_error_parse(fileinfo, err_msg_BIE, ts);

    /* success */
    return 0;
}

int ng_parse_bnode(tokstream* ts, struct ng_bnode* bnode,
				   struct ng* n, struct ng_info* fileinfo)
{
    const char* tok;
    int i;
    double coord;

    int alloc_spos = init_num_spos;
    int alloc_lpos = init_num_lpos;

    /* check that this is a bnode */
    tok = ts_tok(ts);
    if(!is_token(tok, "B"))
        return ng_error_parse(fileinfo, err_msg_nob, ts);

    /* read coordinates */
    for(i = 0; i < n->dim; ++i)
    {
        tok = ts_get(ts);
        if(ng_parser_strtod(tok, &coord))
            return ng_error_parse(fileinfo, err_msg_coo, ts);
        bnode->coords[i] = coord;
    }

    /* initialize surface and line positions */
    bnode->num_spos = 0;
    bnode->spos = malloc(sizeof(struct ng_surface_pos) * alloc_spos);
    if(!bnode->spos)
        return ng_error(fileinfo, err_msg_mem);

    /* initialize line positions */
    bnode->num_lpos = 0;
    bnode->lpos = malloc(sizeof(struct ng_line_pos) * alloc_lpos);
    if(!bnode->lpos)
        return ng_error(fileinfo, err_msg_mem);

    /* scan surface and line positions */
    while((tok = ts_get(ts)))
    {
        /* skip over comments */
        if(tok[0] == '#')
        {
            ts_skipline(ts);
            continue;
        }

        /* check for surface position */
        if(is_token(tok, "S"))
        {
            /* check for reallocation */
            if(bnode->num_spos == alloc_spos)
            {
                alloc_spos *= 2;
                bnode->spos = realloc(bnode->spos, sizeof(struct ng_surface_pos) * alloc_spos);
                if(!bnode->spos)
                    return ng_error(fileinfo, err_msg_mem);
            }

            /* read surface position */
            $(ng_parse_surface_pos(ts, &bnode->spos[bnode->num_spos], fileinfo));
            ++bnode->num_spos;

            continue;
        }

        /* check for line position */
        if(is_token(tok, "L"))
        {
            /* check for reallocation */
            if(bnode->num_lpos == alloc_lpos)
            {
                alloc_lpos *= 2;
                bnode->lpos = realloc(bnode->lpos, sizeof(struct ng_line_pos) * alloc_lpos);
                if(!bnode->lpos)
                    return ng_error(fileinfo, err_msg_mem);
            }

            /* read line position */
            $(ng_parse_line_pos(ts, &bnode->lpos[bnode->num_lpos], fileinfo));
            ++bnode->num_lpos;

            continue;
        }

        /* end surface and line positions */
        break;
    }

    /* finalize surface positions */
    bnode->spos = realloc(bnode->spos, sizeof(struct ng_surface_pos) * bnode->num_spos);

    /* finalize line positions */
    bnode->lpos = realloc(bnode->lpos, sizeof(struct ng_line_pos) * bnode->num_lpos);

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return ng_error_parse(fileinfo, err_msg_tok, ts);

    /* read delimiter */
    tok = ts_tok(ts);
    if(!is_token(tok, ";"))
        return ng_error_parse(fileinfo, err_msg_SLS, ts);

    /* success */
    return 0;
}

int ng_parse_surface_pos(tokstream* ts, struct ng_surface_pos* spos, struct ng_info* fileinfo)
{
    const char* tok;
    int i;

    /* check that this is a surface_pos */
    tok = ts_tok(ts);
    if(!is_token(tok, "S"))
        return ng_error_parse(fileinfo, err_msg_nos, ts);

    /* read surface id */
    tok = ts_get(ts);
    if(ng_parser_strtoi(tok, &spos->surface))
        return ng_error_parse(fileinfo, err_msg_sid, ts);

    /* read surface pos */
    for(i = 0; i < 3; ++i)
    {
        tok = ts_get(ts);
        if(ng_parser_strtod(tok, &spos->pos[i]))
            return ng_error_parse(fileinfo, err_msg_spo, ts);
    }

    /* success */
    return 0;
}

int ng_parse_line_pos(tokstream* ts, struct ng_line_pos* lpos, struct ng_info* fileinfo)
{
    const char* tok;

    /* check that this is a line_pos */
    tok = ts_tok(ts);
    if(!is_token(tok, "L"))
        return ng_error_parse(fileinfo, err_msg_nol, ts);

    /* read line id */
    tok = ts_get(ts);
    if(ng_parser_strtoi(tok, &lpos->line))
        return ng_error_parse(fileinfo, err_msg_lid, ts);

    /* read line pos */
    tok = ts_get(ts);
    if(ng_parser_strtod(tok, &lpos->pos))
        return ng_error_parse(fileinfo, err_msg_lpo, ts);

    /* success */
    return 0;
}

int ng_parse_inode(tokstream* ts, struct ng_inode* bnode,
				   struct ng* n, struct ng_info* fileinfo)
{
    const char* tok;
    int i;

    /* check that this is an inode */
    tok = ts_tok(ts);
    if(!is_token(tok, "I"))
        return ng_error_parse(fileinfo, err_msg_noi, ts);

    /* read coordinates */
    for(i = 0; i < n->dim; ++i)
    {
        tok = ts_get(ts);
        if(ng_parser_strtod(tok, &bnode->coords[i]))
            return ng_error_parse(fileinfo, err_msg_coo, ts);
    }

    /* read delimiter */
    tok = ts_get(ts);
    if(!is_token(tok, ";"))
        return ng_error_parse(fileinfo, err_msg_sem, ts);

    /* success */
    return 0;
}

int ng_parse_element(tokstream* ts, struct ng_element* element, struct ng_info* fileinfo, int dim)
{
    const char* tok;
    int subdomain_id;
    int node_id;

    int alloc_nodes = init_num_nodes;
    int alloc_faces = init_num_faces;

    /* check that this is an element */
    tok = ts_tok(ts);
    if(!is_token(tok, "E"))
        return ng_error_parse(fileinfo, err_msg_noe, ts);

    /* read subdomain id */
    tok = ts_get(ts);
    if(ng_parser_strtoi(tok, &subdomain_id))
        return ng_error_parse(fileinfo, err_msg_uid, ts);
    element->subdomain = subdomain_id;

    /* initialize nodes */
    element->num_nodes = 0;
    element->nodes = malloc(sizeof(int) * alloc_nodes);
    if(!element->nodes)
        return ng_error(fileinfo, err_msg_mem);

    /* scan nodes */
    while((tok = ts_get(ts)))
    {
        /* skip over comments */
        if(tok[0] == '#')
        {
            ts_skipline(ts);
            continue;
        }

        /* check for node id */
        if(!ng_parser_strtoi(tok, &node_id))
        {
            /* check for reallocation */
            if(element->num_nodes == alloc_nodes)
            {
                alloc_nodes *= 2;
                element->nodes = realloc(element->nodes, sizeof(int) * alloc_nodes);
                if(!element->nodes)
                    return ng_error(fileinfo, err_msg_mem);
            }

            /* store node id */
            element->nodes[element->num_nodes] = node_id;
            ++element->num_nodes;

            continue;
        }

        /* or end nodes */
        ts_unget(ts);
        break;
    }

    /* finalize nodes */
    element->nodes = realloc(element->nodes, sizeof(int) * element->num_nodes);

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return ng_error_parse(fileinfo, err_msg_tok, ts);

    /* initialize faces */
    element->num_faces = 0;
    element->faces = malloc(sizeof(struct ng_face) * alloc_faces);
    if(!element->faces)
        return ng_error(fileinfo, err_msg_mem);

    /* scan faces */
    while((tok = ts_get(ts)))
    {
        /* skip over comments */
        if(tok[0] == '#')
        {
            ts_skipline(ts);
            continue;
        }

        /* check for face */
        if(is_token(tok, "F") || is_token(tok, "S"))
        {
            /* check for reallocation */
            if(element->num_faces == alloc_faces)
            {
                alloc_faces *= 2;
                element->faces = realloc(element->faces, sizeof(struct ng_face) * alloc_faces);
                if(!element->faces)
                    return ng_error(fileinfo, err_msg_mem);
            }

            /* read face */
            $(ng_parse_face(ts, &element->faces[element->num_faces], fileinfo));
            ++element->num_faces;

            continue;
        }

        /* or end faces */
        ts_unget(ts);
        break;
    }

    /* finalize faces */
    element->faces = realloc(element->faces, sizeof(struct ng_face) * element->num_faces);

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return ng_error_parse(fileinfo, err_msg_tok, ts);

    /* read delimiter */
    tok = ts_get(ts);
    if(!is_token(tok, ";"))
		return ng_error_parse(fileinfo, err_msg_SNF, ts);

    /* success */
    return 0;
}

int ng_parse_face(tokstream* ts, struct ng_face* face, struct ng_info* fileinfo)
{
    const char* tok;
    int node_id;

    int alloc_nodes = init_num_nodes;

    /* check that this is a face or a side*/
    tok = ts_tok(ts);
    if((!is_token(tok, "F")) && (!is_token(tok, "S")))
        return ng_error_parse(fileinfo, err_msg_nof, ts);

    /* initialize nodes */
    face->num_nodes = 0;
    face->nodes = malloc(sizeof(int) * alloc_nodes);
    if(!face->nodes)
        return ng_error(fileinfo, err_msg_mem);

    /* scan nodes */
    while((tok = ts_get(ts)))
    {
        /* skip over comments */
        if(tok[0] == '#')
        {
            ts_skipline(ts);
            continue;
        }

        /* check for node id */
        if(!ng_parser_strtoi(tok, &node_id))
        {
            /* check for reallocation */
            if(face->num_nodes == alloc_nodes)
            {
                alloc_nodes *= 2;
                face->nodes = realloc(face->nodes, sizeof(int) * alloc_nodes);
                if(!face->nodes)
                    return ng_error(fileinfo, err_msg_mem);
            }

            /* store node id */
            face->nodes[face->num_nodes] = node_id;
            ++face->num_nodes;

            continue;
        }

        /* or end face */
        ts_unget(ts);
        break;
    }

    /* finalize nodes */
    face->nodes = realloc(face->nodes, sizeof(int) * face->num_nodes);

    /* check for errors */
    if(!tok && !ts_eof(ts))
        return ng_error_parse(fileinfo, err_msg_tok, ts);

    /* success */
    return 0;
}


/*
 * internal implementation
 */

int ng_parser_strtoi(const char* str, int* i)
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

int ng_parser_strtod(const char* str, double* d)
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
