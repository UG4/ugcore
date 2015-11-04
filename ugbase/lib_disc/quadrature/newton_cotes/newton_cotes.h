
#include "../quadrature.h"

#ifndef __H__UG__LIB_DISC__QUADRATURE__NEWTON_COTES__
#define __H__UG__LIB_DISC__QUADRATURE__NEWTON_COTES__

namespace ug{


/**
 * This class provides Newton-Cotes integrals for the 1d line [0,1]. See e.g.
 * http://en.wikipedia.org/wiki/Newton–Cotes_formulas for details.
 *
 * The implemented rules are auto-generate using mathematica. If higher orders
 * are needed, rerun the corresponding file.
 */
class NewtonCotes : public QuadratureRule<1>
{
	public:
	/// constructor
		NewtonCotes(size_t order);

	/// destructor
		~NewtonCotes();
};

} // namespace ug

#endif /* __H__UG__LIB_DISC__QUADRATURE__NEWTON_COTES__ */
