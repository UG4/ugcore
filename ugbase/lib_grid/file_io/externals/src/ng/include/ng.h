/**
 * \defgroup libng libng
 *
 * \brief Module for reading Netgen Ng files.
 *
 * \example ng_test.c
 * Example of standard ng operations.
 */

/**
 * \file ng.h
 *
 * \brief The main libng header file.
 * \ingroup libng
 */

#ifndef NG_H_
#define NG_H_

#include "ng_info.h"

#include "ng_node.h"
#include "ng_element.h"

/**
 * \brief main ng data structure
 * \ingroup libng
 *
 * This structure contains the top-level ng data.
 *
 * An instance should only be allocated dynamically by means of ng_new(),
 * due to internal initialization. Subsequently, the disposal of an instance
 * is to be done by calling ng_delete(), which triggers internal cleanup.
 *
 * Reading from or writing to a file can be done by using the ng_read() and
 * ng_write() pair of functions.
 *
 * \sa ng_new(), ng_delete()
 * \sa ng_read(), ng_write()
 */
struct ng
{
    /** \brief dimension of lgm, should be 2 or 3 */
    int dim;
	
    /** \brief the number of boundary nodes */
    int num_bnodes;

    /** \brief the boundary nodes */
    struct ng_bnode* bnodes;

    /** \brief the number of internal nodes */
    int num_inodes;

    /** \brief the internal nodes */
    struct ng_inode* inodes;

    /** \brief the number of elements */
    int num_elements;

    /** \brief the elements */
    struct ng_element* elements;
};

/**
 * \brief allocate and initialize a new ng instance
 * \relates ng
 * \ingroup libng
 */
struct ng* ng_new(void);

/**
 * \brief dispose and clean a no longer needed ng instance
 * \relates ng
 * \ingroup libng
 */
void ng_delete(struct ng* n);

/**
 * \brief read ng from file
 * \relates ng
 * \ingroup libng
 *
 * \param filename The file to be read.
 * \param n The ng object to contain the file data.
 * \param fileinfo (optional) A ng_info object for status reporting, or nullptr.
 *
 * \returns A non-zero value indicates an error.
 */
int ng_read(const char* filename, struct ng* n, struct ng_info* fileinfo);

/**
 * \brief write ng to file
 * \relates ng
 * \ingroup libng
 *
 * \param filename The file to be written to.
 * \param n The ng to be written to file.
 * \param fileinfo (optional) A ng_info object for status reporting, or nullptr.
 *
 * \returns A non-zero value indicates an error.
 */
int ng_write(const char* filename, const struct ng* n, struct ng_info* fileinfo);

#endif
