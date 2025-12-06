
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <lgm.h>
#include <lgm_info.h>

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
    "lgm_test - lgm read timing\n"
    "--------------------------\n"
    "usage: lgm_test [lgm_file]\n"
;

int main(int argc, char* argv[])
{
    char filename[128];

    struct lgm* l;
    struct lgm_info* linfo;

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

    l = lgm_new();
    linfo = lgm_info_new();

    fprintf(stdout, "Reading...\n");

    t1 = mtime();
    if(lgm_read(filename, l, linfo))
    {
        fprintf(stdout, "error!\n");
        fprintf(stderr, "Could not read file '%s'.\n", filename);
        fprintf(stderr, "Error: %s\n", linfo->err_msg);

        lgm_delete(l);
        lgm_info_delete(linfo);

        return EXIT_FAILURE;
    }
    t2 = mtime();

    fprintf(stdout, "done!\n\n");

    dtRead = t2 - t1;

    fprintf(stdout, "Read time:  %d ms.\n\n", dtRead);

    fprintf(stdout, "Writing...\n");

    t1 = mtime();
    if(lgm_write("lgm_test.lgm", l, linfo))
    {
        fprintf(stdout, "error!\n");
        fprintf(stderr, "Could not write file '%s'.\n", "lgm_test.lgm");

        lgm_delete(l);
        lgm_info_delete(linfo);

        return EXIT_FAILURE;
    }
    t2 = mtime();

    fprintf(stdout, "done!\n\n");

    dtWrite = t2 - t1;

    fprintf(stdout, "Write time: %d ms.\n\n", dtWrite);

    fprintf(stdout, "Total time: %d ms\n\n", dtRead + dtWrite);

    fprintf(stdout, "lgm units:    %u\n", l->num_subdomains);
    fprintf(stdout, "lgm lines:    %u\n", l->num_lines);
    fprintf(stdout, "lgm surfaces: %u\n", l->num_surfaces);
    fprintf(stdout, "lgm points:   %u\n", l->num_points);

    lgm_delete(l);
    lgm_info_delete(linfo);

    return EXIT_SUCCESS;
}
