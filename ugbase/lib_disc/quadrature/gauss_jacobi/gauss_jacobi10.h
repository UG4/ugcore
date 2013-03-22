/*
 * gauss_jacobi.h
 *
 *  Created on: 15.03.2013
 *      Author: lisagrau
 */

#include "../quadrature.h"

#ifndef __H__UG__LIB_DISC__QUADRATURE__GAUSS_JACOBI10__
#define __H__UG__LIB_DISC__QUADRATURE__GAUSS_JACOBI10__

namespace ug{

	/**This class provides GaussJacobi integrals up to order 70
	 * with alpha = 1 and beta = 0. Further information about these
	 * quadrature rules can be found on wikipedia.
	 */

	class GaussJacobi10 : public QuadratureRule<1>
	{
		public:
		//constructor
		GaussJacobi10(int order);

		//destructor
		~GaussJacobi10();
	};

} // namespace ug

#endif /* __H__UG__LIB_DISC__QUADRATURE__GAUSS_JACOBI10__ */
