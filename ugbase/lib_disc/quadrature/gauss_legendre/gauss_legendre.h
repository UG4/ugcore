
#include "../quadrature.h"

#ifndef __H__UG__LIB_DISC__QUADRATURE__GAUSS_LEGENDRE__
#define __H__UG__LIB_DISC__QUADRATURE__GAUSS_LEGENDRE__

namespace ug{

/**
 * This class provides GaussLegendre integrals up to order 70.
 * For further information, see e.g.,
 * http://en.wikipedia.org/wiki/Gaussian_quadrature
 */
class GaussLegendre : public QuadratureRule<1>
{
	public:
	///	constructor
		GaussLegendre(size_t order);

	///	destructor
		~GaussLegendre();
};

} // namespace ug

#endif /* __H__UG__LIB_DISC__QUADRATURE__GAUSS_LEGENDRE__ */
