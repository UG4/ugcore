/**
 * \file lgm_surface.h
 *
 * \brief Surface-related declarations.
 * \ingroup liblgm
 */

#ifndef LGM_SURFACE_H_
#define LGM_SURFACE_H_

/**
 * \brief lgm surface data structure
 * \ingroup liblgm
 */
struct lgm_surface
{
    /** \brief the id of the left bounding subdomain */
    int left;
    /** \brief the id of the right bounding subdomain */
    int right;

    /** \brief the number of points for this surface */
    int num_points;
    /** \brief the points for this surface */
    int* points;

    /** \brief the number of lines for this surface */
    int num_lines;
    /** \brief the lines for this surface */
    int* lines;

    /** \brief the number of triangles for this surface */
    int num_triangles;
    /**
     * \brief the triangles for this surface
     *
     * \note This is an array of int[3], access to it must use
     * twofold indication.
     */
    int (*triangles)[3];
};

#endif
