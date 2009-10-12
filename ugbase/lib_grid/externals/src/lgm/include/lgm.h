/**
 * \defgroup liblgm liblgm
 *
 * \brief Module for reading UG's Lgm files.
 *
 * \example lgm_test.c
 * Example of standard lgm operations.
 */

/**
 * \file lgm.h
 *
 * \brief The main liblgm header file.
 * \ingroup liblgm
 */

#ifndef LGM_H_
#define LGM_H_

#include "lgm_info.h"

#include "lgm_line.h"
#include "lgm_surface.h"

/**
 * \brief main lgm data structure
 * \ingroup liblgm
 *
 * This structure contains the top-level lgm data.
 *
 * An instance should only be allocated dynamically by means of lgm_new(),
 * due to internal initialization. Subsequently, the disposal of an instance
 * is to be done by calling lgm_delete(), which triggers internal cleanup.
 *
 * Reading from or writing to a file can be done by using the lgm_read() and
 * lgm_write() pair of functions.
 *
 * \sa lgm_new(), lgm_delete()
 * \sa lgm_read(), lgm_write()
 */
struct lgm
{
    /** \brief name of the lgm */
    const char* name;

    /** \brief name of the problem */
    const char* problemname;

    /** \brief convexity flag, should be 0 or 1 */
    int convex;

    /** \brief dimension of lgm, should be 2 or 3 */
    int dim;

    /** \brief number of subdomains */
    int num_subdomains;

    /** \brief array of subdomain names */
    const char** subdomains;

    /** \brief number of lines */
    int num_lines;

    /** \brief line array */
    struct lgm_line* lines;

    /** \brief number of surfaces */
    int num_surfaces;

    /** \brief surface array */
    struct lgm_surface* surfaces;

    /** \brief number of points */
    int num_points;

    /**
     * \brief points
     *
     * A point consists of the number of coordinates defined by lgm::dim.
     *
     * The points are stored as a two-dimensional array. To access the j'th
     * coordinate of the i'th point, one would use code like this:
     * \code
     * mylgm->points[i][j]
     * \endcode
     *
     * Here, i must be in the range of [0, num_points) and j must be in
     * [0, dim).
     *
     * The allocated memory for the coordinates is contiguous; that means that
     * a one-dimensional array of all coordinates can be accessed like this:
     * \code
     * double* coords = points[0];
     * \endcode
     * This can be useful for passing the coordinates to external algorithms.
     */
    double** points;
};

/**
 * \brief allocate and initialize a new lgm instance
 * \relates lgm
 * \ingroup liblgm
 */
struct lgm* lgm_new(void);

/**
 * \brief dispose and clean a no longer needed lgm instance
 * \relates lgm
 * \ingroup liblgm
 */
void lgm_delete(struct lgm* l);

/**
 * \brief read lgm from file
 * \relates lgm
 * \ingroup liblgm
 *
 * \param filename The file to be read.
 * \param l The lgm object to contain the file data.
 * \param fileinfo (optional) A lgm_info object for status reporting, or NULL.
 *
 * \returns A non-zero value indicates an error.
 */
int lgm_read(const char* filename, struct lgm* l, struct lgm_info* fileinfo);

/**
 * \brief write lgm to file
 * \relates lgm
 * \ingroup liblgm
 *
 * \param filename The file to be written to.
 * \param l The lgm to be written to file.
 * \param fileinfo (optional) A lgm_info object for status reporting, or NULL.
 *
 * \returns A non-zero value indicates an error.
 */
int lgm_write(const char* filename, const struct lgm* l, struct lgm_info* fileinfo);

#endif /*LGM_H_*/
