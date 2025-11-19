/**
 * \file ng_node.h
 *
 * \brief Node-related declarations.
 * \ingroup libng
 */

#ifndef NG_NODE_H_
#define NG_NODE_H_

/**
 * \brief surface position data structure
 * \ingroup libng
 *
 * This structure is used to describe the position of a boundary node
 * on a lgm surface.
 */
struct ng_surface_pos
{
    /** \brief the id of the lgm surface */
    int surface;

    /** \brief the coordinates of the node on the surface */
    double pos[3];
};

/**
 * \brief line position data structure
 * \ingroup libng
 *
 * This structure is used to describe the position of a boundary node
 * on a lgm line.
 */
struct ng_line_pos
{
    /** \brief the id of the lgm line */
    int line;

    /**
     * \brief the position of the node on the line
     *
     * A integer value corresponds to a line point, fractional values
     * can be used for linear interpolation between line points.
     */
    double pos;
};

/**
 * \brief boundary node data structure
 * \ingroup libng
 *
 * This structure contains information about nodes that lie on boundaries
 * as described by the lgm file.
 *
 * The mapping between ng nodes and lgm data is done with the surface position
 * and line position data.
 *
 * \sa ng_bnode::spos, ng_bnode::lpos
 */
struct ng_bnode
{
    /** \brief the coordinates of the node */
    double coords[3];

    /** \brief the number of surface positions */
    int num_spos;
    /** \brief the node's surface position data */
    struct ng_surface_pos* spos;

    /** \brief the number of line positions */
    int num_lpos;
    /** \brief the node's line position data */
    struct ng_line_pos* lpos;
};

/**
 * \brief internal node data structure
 * \ingroup libng
 */
struct ng_inode
{
    /** \brief the coordinates of the node */
    double coords[3];
};

#endif