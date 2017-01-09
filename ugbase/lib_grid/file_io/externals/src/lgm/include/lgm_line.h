/**
 * \file lgm_line.h
 *
 * \brief Line-related declarations.
 * \ingroup liblgm
 */

#ifndef LGM_LINE_H_
#define LGM_LINE_H_

/**
 * \brief lgm line data structure
 * \ingroup liblgm
 */
struct lgm_line
{
    /** \brief the number of line points */
    int num_points;
    /** \brief the line point array */
    int* points;
};

#endif
