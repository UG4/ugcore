/*
 * function_base.h
 *
 *  Created on: 30.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__FUNCTION_BASE__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__FUNCTION_BASE__

#include "common/common.h"
#include <vector>

namespace ug{

/// base class for functions used in operator interface
/**
 * This function is used as a base class for functions/vectors in the IOperator
 * interfaces. It currently provides a single norm-method. Thus, this base class
 * can be used, whenever only norm dependent algorithm are to be implemented
 * in a non-templated way (e.g. IConvergenceCheck)
 */
class IFunctionBase
{
	public:
	///	returns the two norm of the function
		virtual number two_norm() = 0;

	/// returns the two norm of the function on the specified indices
		virtual number two_norm(std::vector<size_t>& ind) = 0;

	///	virtual destructor
		virtual ~IFunctionBase() {}
};

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__INTERFACE__FUNCTION_BASE__ */
