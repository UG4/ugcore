/**
 * \file ng_element.h
 *
 * \brief Element-related declarations.
 * \ingroup libng
 */

#ifndef NG_ELEMENT_H_
#define NG_ELEMENT_H_

/**
 * \brief face data structure
 * \ingroup libng
 */
struct ng_face
{
    /** \brief the number of nodes in this face */
    int num_nodes;
    /** \brief the ids of the nodes that make up the face */
    int* nodes;
};

/**
 * \brief element data structure
 * \ingroup libng
 *
 * An element represents a volume entity for this ng file.
 * The actual type of the volume is determined by the number
 * of nodes and faces.
 */
struct ng_element
{
    /** \brief the subdomain this element belongs to */
    int subdomain;

    /** \brief the number of nodes the element contains */
    int num_nodes;
    /** \brief the array of node ids for this element */
    int* nodes;

    /** \brief the number of faces the element contains */
    int num_faces;
    /** \brief the faces for this element */
    struct ng_face* faces;
};

#endif