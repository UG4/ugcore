/*
 * quadrature_provider.cpp
 *
 *  Created on: 15.02.2010
 *      Author: andreasvogel
 */

#include "quadrature_provider.h"
#include "gauss_quad/gauss_quad.h"

namespace ug{

template <int TDim>
QuadratureRuleProvider<TDim>::QuadratureRuleProvider()
{
//	register standard rules
	RegisterGaussQuadRuleDim<dim>(*this);
}

template <int TDim>
const QuadratureRule<TDim>&
QuadratureRuleProvider<TDim>::get_quad_rule(ReferenceObjectID roid,
                                            size_t order,
                                            QuadratureType type)
{
	//	check if order or higerh order registered
	if(order >= m_vRule[type][roid].size())
		UG_THROW("QuadratureRuleProvider: Quadrature Rule not found for "
				<<roid<<" (dim="<<dim<<") and order "<<order);

	//	look for rule of order or next higher one
	if(m_vRule[type][roid][order] == NULL)
	{
		for(size_t i = order + 1; i < m_vRule[type][roid].size(); ++i)
		{
			//	return higher order than requested
			if(m_vRule[type][roid][i] != NULL) return *m_vRule[type][roid][i];
		}
		UG_THROW("QuadratureRuleProvider: Quadrature Rule not found for "
				<<roid<<" (dim="<<dim<<") and order "<<order);
	}

	//	return correct order
	return *m_vRule[type][roid][order];
}

template <int TDim>
void
QuadratureRuleProvider<TDim>::register_rule(ReferenceObjectID roid,
                                            const QuadratureRule<dim>& rule,
                                            QuadratureType type)
{
//	get order of rule to register
	size_t order = rule.order();

//	resize vector if needed
	if(m_vRule[type][roid].size() <= order) m_vRule[type][roid].resize(order+1, NULL);

//	set or override rule
	m_vRule[type][roid][order] = &rule;
}

template <int TDim>
const QuadratureRule<TDim>&
QuadratureRuleProvider<TDim>::get_rule(ReferenceObjectID roid, size_t order,
                                       QuadratureType type)
{
//	forward request
	return instance().get_quad_rule(roid, order, type);
}

/// writes the Identifier to the output stream
template <int TDim>
std::ostream& operator<<(std::ostream& out,	const QuadratureType& v)
{
	std::stringstream ss;

	switch(v)
	{
		case GAUSS: out << "Gauss"; break;
		case NEWTON_COTES: out << "Newton-Cotes"; break;
		default: out << "invalid";
	}
	return out;
}

////////////////////////////////////////////////////////////////////////////////
// explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class QuadratureRuleProvider<0>;
template class QuadratureRuleProvider<1>;
template class QuadratureRuleProvider<2>;
template class QuadratureRuleProvider<3>;

} // namespace ug
