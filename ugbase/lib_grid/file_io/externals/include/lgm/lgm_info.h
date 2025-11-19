/**
 * \file lgm_info.h
 *
 * \brief Header file for Lgm Status Information.
 * \ingroup liblgm
 */

#ifndef LGM_INFO_H_
#define LGM_INFO_H_

/**
 * \brief status information about lgm objects
 * \ingroup liblgm
 *
 * The purpose of this structure is to contain information about a lgm object
 * as it occurs while ie. reading or writing files.
 *
 * Instances of this struct should be created and deleted with
 * the lgm_info_new() and lgm_info_delete() set of functions.
 *
 * \sa lgm_info_new(), lgm_info_delete()
 */
struct lgm_info
{
    /**
     * \brief error flag
     *
     * A non-zero value indicates that an error occured.
     *
     * \sa err_msg
     */
    int error;

    /**
     * \brief error message
     *
     * In case of an error, this is a human-readable error description.
     */
    const char* err_msg;
};

/**
 * \brief allocate and initialize a new lgm_info instance
 * \relates lgm_info
 * \ingroup liblgm
 */
struct lgm_info* lgm_info_new(void);

/**
 * \brief dispose and clean a no longer needed lgm_info instance
 * \relates lgm_info
 * \ingroup liblgm
 */
void lgm_info_delete(struct lgm_info* linfo);

#endif