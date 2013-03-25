/*
 * quadrature_provider_impl.h
 *
 *  Created on: 15.02.2010
 *      Author: andreasvogel
 */

#include "quadrature_provider.h"

namespace ug{

template <int TDim>
template <typename TRefElem>
const QuadratureRule<TDim>&
QuadratureRuleProvider<TDim>::get(size_t order, QuadType type)
{
//	check that dimension is correct
	if(TRefElem::dim != dim)
		UG_THROW("QuadratureRuleProvider: requesting by reference"
				" element, but at provider of different dimension.");

//	get reference object id
	ReferenceObjectID roid = TRefElem::REFERENCE_OBJECT_ID;

//	forward request
	return instance().get_quad_rule(roid, order, type);
}

} // namespace ug
