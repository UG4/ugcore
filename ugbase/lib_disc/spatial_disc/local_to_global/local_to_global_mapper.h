/*
 * local_to_global_mapper.h
 *
 *  Created on: 19.02.2013
 *      Author: susannehoellbacher
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__LOCAL_TO_GLOBAL_MAPPER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__LOCAL_TO_GLOBAL_MAPPER__

// extern headers
#include <vector>

// intern headers
#include "lib_disc/assemble_interface.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/dof_manager/dof_distribution.h"

namespace ug{


/// interface for definition of special LocalToGlobal mappings
/**
 * \tparam	TAlgebra			type of Algebra
 */

template <typename TAlgebra>
class ILocalToGlobalMapper
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	///	default Constructor
		ILocalToGlobalMapper() {}

	///	send local entries to global matrix
		virtual void AddLocalMatrixToGlobal(ConstSmartPtr<DoFDistribution> dd, matrix_type& mat, const LocalMatrix& lmat) = 0;

	///	send local entries to global rhs
		virtual void AddLocalVector(ConstSmartPtr<DoFDistribution> dd, vector_type& vec, const LocalVector& lvec) = 0;

	///	virtual destructor
		virtual ~ILocalToGlobalMapper() {};
};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__LOCAL_TO_GLOBAL_MAPPER__*/
