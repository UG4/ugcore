/*
 * operator_base_interface.h
 *
 *  Created on: 30.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__OPERATOR_BASE_INTERFACE__
#define __H__LIB_ALGEBRA__OPERATOR__OPERATOR_BASE_INTERFACE__

#include "common/common.h"

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

	///	virtual destructor
		virtual ~IFunctionBase() {}
};

} // end namespace ug

#endif /* __H__LIB_ALGEBRA__OPERATOR__OPERATOR_BASE_INTERFACE__ */
