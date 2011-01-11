/*
 * spatial_discretization.h
 *
 *  Created on: 29.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__SPACIAL_DISCRETIZATION__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__SPACIAL_DISCRETIZATION__

///////////////////////////////////////////////////////////////
// This file is intended to include all spatial discretizations
///////////////////////////////////////////////////////////////

/**
 * \brief Domain Discretization.
 *
 *	Provides the domain discretization interface and special implementations.
 *
 * \defgroup lib_disc_domain_assemble Domain Assembling
 * \ingroup lib_disc_assemble
 */

// user data
#include "./ip_data/user_data.h"

// import / export
#include "./ip_data/data_export.h"
#include "./ip_data/data_import.h"

// interface
#include "./domain_discretization_interface.h"

// domain discretization
#include "./domain_discretization.h"

// post process
#include "./post_process/post_process.h"

// element discs
#include "./elem_disc/elem_disc.h"

// coupled element discs
//#include "./coupled_elem_disc/coupled_elem_disc.h"

// disc helper
#include "./disc_helper/disc_helper.h"

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__SPACIAL_DISCRETIZATION__ */
