/*
 * err_est_data.h
 *
 *  Created on: 25.02.2014
 *      Author: Dmitriy Logashenko
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ERR_EST_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ERR_EST_DATA__

// extern headers
#include <vector>
#include <string>

// intern headers
#include "lib_grid/tools/surface_view.h"

namespace ug{

/// Base class for error estimator data
/**
 * This virtual class should be the base of particular error estimators
 * implemented in the elem_disc's. Every elem_disc class (not object!) should
 * declare its own derived class for keeping intermediate information that
 * should be accumulated in the computation of the local error estimators.
 * Several objects of the elem_disc class may share the same object of the
 * derived class for a consistent computation of the error estimator.
 */
class IErrEstData
{
	public:

	///	virtual function to allocate data structures for the error estimator
		virtual void alloc_err_est_data (ConstSmartPtr<SurfaceView> spSV, const GridLevel& gl) = 0;
		
	///	virtual function called after the computation of the error estimator data in all the elements
		virtual void summarize_err_est_data () = 0;
		
	///	virtual function to release data structures for the error estimator
		virtual void release_err_est_data () = 0;
};

} // end of namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ERR_EST_DATA__ */

/* End of File */
