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

//This class provides GaussLegendre integrals up to order 6
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
