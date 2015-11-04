#ifndef __H__UG_MATH__EIGENVALUES__
#define __H__UG_MATH__EIGENVALUES__

#include "../ugmath_types.h"

namespace ug
{

/// \addtogroup ugbase_math
/// \{

bool CalculateEigenvalues(const ug::matrix33& mat, number& lambdaMinOut,
						number& lambdaMedOut, number& lambdaMaxOut,
						ug::vector3& evMinOut, ug::vector3& evMedOut,
						ug::vector3& evMaxOut);

// end group ugbase_math
/// \}

}//	end of namespace

#endif
