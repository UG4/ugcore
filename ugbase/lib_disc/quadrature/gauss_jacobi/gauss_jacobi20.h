/*
 * gauss_jacobi.h
 *
 *  Created on: 15.03.2013
 *      Author: lisagrau
 */

#include "../quadrature.h"

#ifndef __H__UG__LIB_DISC__QUADRATURE__GAUSS_JACOBI20__
#define __H__UG__LIB_DISC__QUADRATURE__GAUSS_JACOBI20__

namespace ug{

//This class provides GaussLegendre integrals up to order 6
class GaussJacobi20 : public QuadratureRule<1>
{
	public:
	//constructor
	GaussJacobi20(int order);

	//destructor
	~GaussJacobi20();
};

} // namespace ug

#endif /* __H__UG__LIB_DISC__QUADRATURE__GAUSS_JACOBI20__ */
