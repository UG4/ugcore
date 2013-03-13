/*
 * quadrature_provider.cpp
 *
 *  Created on: 15.02.2010
 *      Author: andreasvogel
 */

#include "quadrature_provider.h"
#include "gauss/gauss_quad.h"
#include "gauss/gauss_quad_vertex.h"
#include "lib_disc/reference_element/reference_element.h"

namespace ug{

template <int TDim>
QuadratureRuleProvider<TDim>::QuadratureRuleProvider()
{
	for(int type = 0; type < NUM_QUADRATURE_TYPES; ++type)
		for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid)
			m_vRule[type][roid].clear();
}

template <int TDim>
QuadratureRuleProvider<TDim>::~QuadratureRuleProvider()
{
	for(int type = 0; type < NUM_QUADRATURE_TYPES; ++type)
		for(int roid = 0; roid < NUM_REFERENCE_OBJECTS; ++roid)
			for(size_t order = 0; order < m_vRule[type][roid].size(); ++order)
				if(m_vRule[type][roid][order] != NULL)
					delete m_vRule[type][roid][order];
}

template <int TDim>
const QuadratureRule<TDim>&
QuadratureRuleProvider<TDim>::get_quad_rule(ReferenceObjectID roid,
                                            size_t order,
                                            QuadratureType type)
{
	//	check if order present, else resize and create
	if(order >= m_vRule[type][roid].size() ||
			m_vRule[type][roid][order] == NULL)
		create_rule(roid, order, type);

	//	return correct order
	return *m_vRule[type][roid][order];
}

template <>
const QuadratureRule<0>*
QuadratureRuleProvider<0>::create_gauss_rule(ReferenceObjectID roid,
                                                size_t order)
{
	QuadratureRule<0>* q = NULL;
	try{
	switch(roid){
		case ROID_VERTEX: q = new GaussQuadratureVertex(); break;
		default: UG_THROW("QuadratureRuleProvider<"<<dim<<">: "<<roid<<" not supported.");
	}
	}catch(...){return NULL;}
	return q;
}

template <>
const QuadratureRule<1>*
QuadratureRuleProvider<1>::create_gauss_rule(ReferenceObjectID roid,
                                                size_t order)
{
	QuadratureRule<1>* q = NULL;
	try{
	switch(roid){
		case ROID_EDGE: q = new FlexGaussQuadrature<ReferenceEdge>(order); break;
		default: UG_THROW("QuadratureRuleProvider<"<<dim<<">: "<<roid<<" not supported.");
	}
	}catch(...){return NULL;}
	return q;
}

template <>
const QuadratureRule<2>*
QuadratureRuleProvider<2>::create_gauss_rule(ReferenceObjectID roid,
                                                size_t order)
{
	QuadratureRule<2>* q = NULL;
	try{
	switch(roid){
		case ROID_TRIANGLE: q = new FlexGaussQuadrature<ReferenceTriangle>(order); break;
		case ROID_QUADRILATERAL: q = new FlexGaussQuadrature<ReferenceQuadrilateral>(order); break;
		default: UG_THROW("QuadratureRuleProvider<"<<dim<<">: "<<roid<<" not supported.");
	}
	}catch(...){return NULL;}
	return q;
}

template <>
const QuadratureRule<3>*
QuadratureRuleProvider<3>::create_gauss_rule(ReferenceObjectID roid,
                                                size_t order)
{
	QuadratureRule<3>* q = NULL;
	try{
	switch(roid){
		case ROID_TETRAHEDRON: q = new FlexGaussQuadrature<ReferenceTetrahedron>(order); break;
		case ROID_PYRAMID: q = new FlexGaussQuadrature<ReferencePyramid>(order); break;
		case ROID_PRISM: q = new FlexGaussQuadrature<ReferencePrism>(order); break;
		case ROID_HEXAHEDRON: q = new FlexGaussQuadrature<ReferenceHexahedron>(order); break;
		default: UG_THROW("QuadratureRuleProvider<"<<dim<<">: "<<roid<<" not supported.");
	}
	}catch(...){return NULL;}
	return q;
}

template <int TDim>
void
QuadratureRuleProvider<TDim>::create_rule(ReferenceObjectID roid,
                                          size_t order,
                                          QuadratureType type)
{
//	resize vector if needed
	if(m_vRule[type][roid].size() <= order) m_vRule[type][roid].resize(order+1, NULL);
	if(m_vRule[type][roid][order] != NULL)
		delete m_vRule[type][roid][order];
	m_vRule[type][roid][order] = NULL;

	switch(type){
		case BEST: {
			// 1. Try GaussQuad
			m_vRule[type][roid][order] = create_gauss_rule(roid, order);
			if(m_vRule[type][roid][order] != NULL) break;
			// 2. Try others ...


		}break;
		case GAUSS: {
			m_vRule[type][roid][order] = create_gauss_rule(roid, order);
		}break;
		case NEWTON_COTES: {

		}break;
		default: UG_THROW("QuadratureRuleProvider<"<<dim<<">: Cannot create rule for "
		                  <<roid<<", order "<<order<<"and type "<<type);
	}

	if(m_vRule[type][roid][order] == NULL)
		UG_THROW("QuadratureRuleProvider<"<<dim<<">: Cannot create rule for "
				                  <<roid<<", order "<<order<<"and type "<<type);
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
std::ostream& operator<<(std::ostream& out,	const typename QuadratureRuleProvider<TDim>::QuadratureType& v)
{
	std::stringstream ss;

	switch(v)
	{
		case QuadratureRuleProvider<TDim>::GAUSS: out << "Gauss"; break;
		case QuadratureRuleProvider<TDim>::NEWTON_COTES: out << "Newton-Cotes"; break;
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
