//	created by Daniel Jungblut, converted to ugmath by Sebastian Reiter

#ifndef __H__UG_MATH__EIGENVALUES__
#define __H__UG_MATH__EIGENVALUES__

#include "../ugmath_types.h"

namespace ug
{

bool CalculateEigenvalues(const matrix33& mat, number& lambdaMinOut,
						number& lambdaMedOut, number& lambdaMaxOut,
						vector3& evMinOut, vector3& evMedOut,
						vector3& evMaxOut);

}//	end of namespace

#endif
