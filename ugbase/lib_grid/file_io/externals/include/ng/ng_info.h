/**
 * \file ng_info.h
 *
 * \brief Header file for Ng Status Information.
 * \ingroup libng
 */

#ifndef NG_INFO_H_
#define NG_INFO_H_

/**
 * \brief status information about ng objects
 * \ingroup libng
 *
 * The purpose of this structure is to contain information about a ng object
 * as it occurs while ie. reading or writing files.
 *
 * Instances of this struct should be created and deleted with
 * the ng_info_new() and ng_info_delete() set of functions.
 *
 * \sa ng_info_new(), ng_info_delete()
 */
struct ng_info
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
 * \brief allocate and initialize a new ng_info instance
 * \relates ng_info
 * \ingroup libng
 */
struct ng_info* ng_info_new();

/**
 * \brief dispose and clean a no longer needed ng_info instance
 * \relates ng_info
 * \ingroup libng
 */
void ng_info_delete(ng_info* ninfo);

#endif
