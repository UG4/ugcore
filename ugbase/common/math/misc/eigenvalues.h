/**
 * \file eigenvalues.h
 *
 * Calculation of Eigenvalues.
 */

//	created by Daniel Jungblut, converted to ugmath by Sebastian Reiter

#ifndef __H__UG_MATH__EIGENVALUES__
#define __H__UG_MATH__EIGENVALUES__

#include "../ugmath_types.h"

namespace ug
{

bool CalculateEigenvalues(const ug::matrix33& mat, number& lambdaMinOut,
						number& lambdaMedOut, number& lambdaMaxOut,
						ug::vector3& evMinOut, ug::vector3& evMedOut,
						ug::vector3& evMaxOut);

}//	end of namespace

#endif
