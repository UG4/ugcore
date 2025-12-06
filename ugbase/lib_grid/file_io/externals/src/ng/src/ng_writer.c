#include "../include/ng.h"

//ø #include <stdlib.h>
#include <stdio.h>
//ø #include <string.h>

#include "../include/ng_node.h"
#include "../include/ng_element.h"
#include "../include/ng_info.h"
#include "ng_error.h"


/*
 * entry tags
 */

static const char* tag_bnode = "B";
static const char* tag_inode = "I";
static const char* tag_element = "E";
static const char* tag_spos = "S";
static const char* tag_lpos = "L";
static const char* tag_face = "F";


/*
 * inlining macros
 */

/* try-like construct, if cmd returns 1 (error), return 1 too */
#define $(cmd) do { if(cmd) return 1; } while (0)


/*
 * error messages
 */

static const char* err_msg_fil = "Error opening file %s.";


/****
 * ng file writing
 */

int ng_write_bnodes(FILE* fp, const struct ng* n /*, struct ng_info* fileinfo*/);
int ng_write_inodes (FILE* fp, const struct ng* n /*, struct ng_info* fileinfo*/);
int ng_write_elements(FILE* fp, const struct ng* n /*, struct ng_info* fileinfo*/);


int ng_write(const char* file, const struct ng* n, struct ng_info* fileinfo)
{
    /* open file */
    FILE* fp = fopen(file, "w");
    if(!fp)
        return ng_error_string(fileinfo, err_msg_fil, file);
    
    /* write ng */
    $(ng_write_bnodes(fp, n));
    $(ng_write_inodes(fp, n));
    $(ng_write_elements(fp, n));
    
    /* close file */
    fclose(fp);
    
    /* success */
    return 0;
}

int ng_write_bnodes(FILE* fp, const struct ng* n /*, struct ng_info* fileinfo*/)
{
    int i, j;
    
    /* write header */
    fputs("# boundary nodes\n\n", fp);
        
    /* write bnodes */
    for(i = 0; i < n->num_bnodes; ++i)
    {
        fprintf(fp, "%s %G %G %G # node %d\n", tag_bnode, n->bnodes[i].coords[0], n->bnodes[i].coords[1], n->bnodes[i].coords[2], i);
        for(j = 0; j < n->bnodes[i].num_spos; ++j)
            fprintf(fp, "%s %d %G %G %G ", tag_spos, n->bnodes[i].spos[j].surface, n->bnodes[i].spos[j].pos[0], n->bnodes[i].spos[j].pos[1], n->bnodes[i].spos[j].pos[2]);
        for(j = 0; j < n->bnodes[i].num_lpos; ++j)
            fprintf(fp, "%s %d %G ", tag_lpos, n->bnodes[i].lpos[j].line, n->bnodes[i].lpos[j].pos);
        fputs(";\n\n", fp);
    }
    
    /* success */
    return 0;
}

int ng_write_inodes(FILE* fp, const struct ng* n /*, struct ng_info* fileinfo*/)
{
    int i;
    
    /* write header */
    fputs("# internal nodes\n\n", fp);
        
    /* write bnodes */
    for(i = 0; i < n->num_inodes; ++i)
        fprintf(fp, "%s %G %G %G ; # node %d\n\n", tag_inode, n->inodes[i].coords[0], n->inodes[i].coords[1], n->inodes[i].coords[2], n->num_bnodes + i);
    
    /* success */
    return 0;
}

int ng_write_elements(FILE* fp, const struct ng* n /*, struct ng_info* fileinfo*/)
{
    int i, j, k;
    
    /* write header */
    fputs("# elements\n\n", fp);
    
    /* write lines */
    for(i = 0; i < n->num_elements; ++i)
    {
        fprintf(fp, "%s %d ", tag_element, n->elements[i].subdomain);
        
        /* write nodes */
        for(j = 0; j < n->elements[i].num_nodes; ++j)
            fprintf(fp, "%d ", n->elements[i].nodes[j]);
        
        /* write faces */
        for(j = 0; j < n->elements[i].num_faces; ++j)
        {
            fprintf(fp, "%s ", tag_face);
            for(k = 0; k < n->elements[i].faces[j].num_nodes; ++k)
                fprintf(fp, "%d ", n->elements[i].faces[j].nodes[k]);
        }
        
        fprintf(fp, "; # element %d\n\n", i);
    }
    
    /* success */
    return 0;
}
