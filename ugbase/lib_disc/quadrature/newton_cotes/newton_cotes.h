/*
 * newton_cotes.h
 *
 *  Created on: 12.03.2013
 *      Author: lisagrau
 */

#include "../quadrature.h"

#ifndef __H__UG__LIB_DISC__QUADRATURE__NEWTON_COTES__
#define __H__UG__LIB_DISC__QUADRATURE__NEWTON_COTES__

namespace ug{
	//This class provides Newton-Cotes integrals up to order 6
class NewtonCotes : public QuadratureRule<1>
{
	public:
	//constructor
		NewtonCotes(int order);
	//destructor
		~NewtonCotes();
};

} // namespace ug

#endif /* __H__UG__LIB_DISC__QUADRATURE__NEWTON_COTES__ */
