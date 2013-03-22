/*
 * gauss_legendre.h
 *
 *  Created on: 12.03.2013
 *      Author: lisagrau
 */

#include "../quadrature.h"

#ifndef __H__UG__LIB_DISC__QUADRATURE__GAUSS_LEGENDRE__
#define __H__UG__LIB_DISC__QUADRATURE__GAUSS_LEGENDRE__

namespace ug{

	/** This class provides GaussLegendre integrals up to order 70
	 *  furher information considering these rules an be found on wikipedia.
	 */

	class GaussLegendre : public QuadratureRule<1>
	{
		public:
		//constructor
		GaussLegendre(int order);

		//destructor
		~GaussLegendre();
	};

} // namespace ug

#endif /* __H__UG__LIB_DISC__QUADRATURE__GAUSS_LEGENDRE__ */
