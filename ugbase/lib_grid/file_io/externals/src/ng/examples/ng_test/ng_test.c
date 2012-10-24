#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <ng.h>

/*
 * get means for measuring milliseconds
 */
#ifdef WIN32
/* windows has GetTickCount */
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
int mtime(void)
{
    return GetTickCount();
}
#else
/* unix has gettimeofday */
#include <sys/time.h>
int mtime(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (int)(tv.tv_sec*1000+tv.tv_usec/1000);
}
#endif

const char usage[] =
    "ng_test - ng read timing\n"
    "------------------------\n"
    "usage: ng_test [ng_file]\n"
;

int main(int argc, char* argv[])
{
    char filename[128];

    struct ng* n;
    struct ng_info* ninfo;

    int t1, t2, dtRead, dtWrite;

    if(argc < 2)
    {
        int len;
        fprintf(stdout, "File> ");
        fgets(filename, 128, stdin);
        len = strlen(filename);
        if(filename[len-1] == '\n')
            filename[len-1] = '\0';
    }
    else
        strncpy(filename, argv[1], 128);

    n = ng_new();
    ninfo = ng_info_new();

    fprintf(stdout, "Reading...\n");

    t1 = mtime();
    if(ng_read(filename, n, ninfo))
    {
        fprintf(stdout, "error!\n");
        fprintf(stderr, "Could not read file '%s'.\n", filename);
        fprintf(stderr, "Error: %s\n", ninfo->err_msg);

        ng_delete(n);
        ng_info_delete(ninfo);

        return EXIT_FAILURE;
    }
    t2 = mtime();

    fprintf(stdout, "done!\n\n");

    dtRead = t2 - t1;

    fprintf(stdout, "Reading took %d ms.\n\n", dtRead);

    fprintf(stdout, "Writing...\n");

    t1 = mtime();
    if(ng_write("ng_test.ng", n, ninfo))
    {
        fprintf(stdout, "error!\n");
        fprintf(stderr, "Could not write file '%s'.\n", "ng_test.ng");
        fprintf(stderr, "Error: %s\n", ninfo->err_msg);

        ng_delete(n);
        ng_info_delete(ninfo);

        return EXIT_FAILURE;
    }
    t2 = mtime();

    fprintf(stdout, "done!\n\n");

    dtWrite = t2 - t1;

    fprintf(stdout, "Writing took %d ms.\n\n", dtWrite);

    fprintf(stdout, "Total time: %d ms\n\n", dtRead + dtWrite);

    fprintf(stdout, "ng bnodes:    %u\n", n->num_bnodes);
    fprintf(stdout, "ng inodes:    %u\n", n->num_inodes);
    fprintf(stdout, "ng elements:  %u\n", n->num_elements);

    fprintf(stdout, "\n");

    if(n->num_bnodes > 0)
    {
        struct ng_bnode* bn = &n->bnodes[0];

        fprintf(stdout, "bnode 0:\n");
        fprintf(stdout, "  Coords:   %lf %lf %lf\n", bn->coords[0], bn->coords[1], bn->coords[2]);
    }

    if(n->num_inodes > 0)
    {
        struct ng_inode* in = &n->inodes[0];

        fprintf(stdout, "inode 0:\n");
        fprintf(stdout, "  Coords:   %lf %lf %lf\n", in->coords[0], in->coords[1], in->coords[2]);
    }

    if(n->num_elements > 0)
    {
        struct ng_element* el = &n->elements[0];

        fprintf(stdout, "element 0:\n");
        fprintf(stdout, "  Node #:   %d\n", el->num_nodes);
        fprintf(stdout, "  Face #:   %d\n", el->num_faces);

        if(el->num_faces > 0)
        {
            struct ng_face* f = &el->faces[0];

            fprintf(stdout, "  Face 0:\n");
            fprintf(stdout, "    Node #:   %d\n", f->num_nodes);
        }
    }

    ng_delete(n);
    ng_info_delete(ninfo);

    return EXIT_SUCCESS;
}
