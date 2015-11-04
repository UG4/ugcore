
#include "common/util/string_util.h"
#include "quadrature_provider.h"
#include "gauss/gauss_quad.h"
#include "gauss/gauss_quad_vertex.h"
#include "newton_cotes/newton_cotes.h"
#include "gauss_legendre/gauss_legendre.h"
#include "gauss_jacobi/gauss_jacobi10.h"
#include "gauss_jacobi/gauss_jacobi20.h"
#include "gauss_tensor_prod/gauss_tensor_prod.h"
#include "lib_disc/reference_element/reference_element.h"
#include <algorithm>
#include <locale>

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// gauss
////////////////////////////////////////////////////////////////////////////////

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
		case ROID_OCTAHEDRON: q = new FlexGaussQuadrature<ReferenceOctahedron>(order); break;
		default: UG_THROW("QuadratureRuleProvider<"<<dim<<">: "<<roid<<" not supported.");
	}
	}catch(...){return NULL;}
	return q;
}

////////////////////////////////////////////////////////////////////////////////
// gauss-legendre
////////////////////////////////////////////////////////////////////////////////

template <>
const QuadratureRule<0>*
QuadratureRuleProvider<0>::create_gauss_legendre_rule(ReferenceObjectID roid, size_t order)
{
	return NULL;
}

template <>
const QuadratureRule<1>*
QuadratureRuleProvider<1>::create_gauss_legendre_rule(ReferenceObjectID roid, size_t order)
{
	QuadratureRule<1>* q = NULL;
	try{
		q = new GaussLegendre(order);
	}catch(...){return NULL;}
	return q;
}

template <>
const QuadratureRule<2>*
QuadratureRuleProvider<2>::create_gauss_legendre_rule(ReferenceObjectID roid,
                                                size_t order)
{
	QuadratureRule<2>* q = NULL;
	try{
	switch(roid){
		case ROID_QUADRILATERAL: q = new GaussQuadratureQuadrilateral(order); break;
		case ROID_TRIANGLE: q = new GaussQuadratureTriangle(order); break;
		default: UG_THROW("QuadratureRuleProvider<"<<dim<<">: "<<roid<<" not supported.");
	}
	}catch(...){return NULL;}
	return q;
}

template <>
const QuadratureRule<3>*
QuadratureRuleProvider<3>::create_gauss_legendre_rule(ReferenceObjectID roid,
                                                size_t order)
{
	QuadratureRule<3>* q = NULL;
	try{
	switch(roid){
		case ROID_TETRAHEDRON: q = new GaussQuadratureTetrahedron(order); break;
		case ROID_PRISM: q = new GaussQuadraturePrism(order); break;
		case ROID_PYRAMID: q = new GaussQuadraturePyramid(order); break;
		case ROID_HEXAHEDRON: q = new GaussQuadratureHexahedron(order); break;
		case ROID_OCTAHEDRON: q = new GaussQuadratureOctahedron(order); break;
		default: UG_THROW("QuadratureRuleProvider<"<<dim<<">: "<<roid<<" not supported.");
	}
	}catch(...){return NULL;}
	return q;
}

////////////////////////////////////////////////////////////////////////////////
// newton-cotes
////////////////////////////////////////////////////////////////////////////////

template <>
const QuadratureRule<1>*
QuadratureRuleProvider<1>::create_newton_cotes_rule(ReferenceObjectID roid, size_t order)
{
	QuadratureRule<1>* q = NULL;
	try{
		q = new NewtonCotes(order);
	}catch(...){return NULL;}
	return q;
}

template <int TDim>
const QuadratureRule<TDim>*
QuadratureRuleProvider<TDim>::create_newton_cotes_rule(ReferenceObjectID roid, size_t order)
{
	return NULL;
}

////////////////////////////////////////////////////////////////////////////////
// general
////////////////////////////////////////////////////////////////////////////////

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
                                            QuadType type)
{
	//	check if order present, else resize and create
	if(order >= m_vRule[type][roid].size() ||
			m_vRule[type][roid][order] == NULL)
		create_rule(roid, order, type);

	//	return correct order
	return *m_vRule[type][roid][order];
}

template <int TDim>
void
QuadratureRuleProvider<TDim>::create_rule(ReferenceObjectID roid,
                                          size_t order,
                                          QuadType type)
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

			// 2. Try Newton-Cotes
				m_vRule[type][roid][order] = create_newton_cotes_rule(roid, order);
			if(m_vRule[type][roid][order] != NULL) break;

			// 3. Try Gauss-Legendre
				m_vRule[type][roid][order] = create_gauss_legendre_rule(roid, order);
			if(m_vRule[type][roid][order] != NULL) break;
		}break;
		case GAUSS: {
			m_vRule[type][roid][order] = create_gauss_rule(roid, order);
		}break;
		case GAUSS_LEGENDRE: {
			m_vRule[type][roid][order] = create_gauss_legendre_rule(roid, order);
		}break;
		case NEWTON_COTES: {
			m_vRule[type][roid][order] = create_newton_cotes_rule(roid, order);
		}break;
		default: UG_THROW("QuadratureRuleProvider<"<<dim<<">: Cannot create rule for "
		                  <<roid<<", order "<<order<<" and type "<<type);
	}

	if(m_vRule[type][roid][order] == NULL)
		UG_THROW("QuadratureRuleProvider<"<<dim<<">: Cannot create rule for "
				                  <<roid<<", order "<<order<<" and type "<<type);
}

template <int TDim>
const QuadratureRule<TDim>&
QuadratureRuleProvider<TDim>::get(ReferenceObjectID roid, size_t order,
                                       QuadType type)
{
//	forward request
	return instance().get_quad_rule(roid, order, type);
}

std::ostream& operator<<(std::ostream& out,	const QuadType& v)
{
	switch(v)
	{
		case BEST: out << "Best"; break;
		case GAUSS: out << "Gauss"; break;
		case NEWTON_COTES: out << "Newton-Cotes"; break;
		case GAUSS_LEGENDRE: out << "Gauss-Legendre"; break;
		default: out << "invalid";
	}
	return out;
}

QuadType GetQuadratureType(const std::string& name)
{
	std::string n = TrimString(name);
	std::transform(n.begin(), n.end(), n.begin(), ::tolower);
	if(n == "best") return BEST;
	if(n == "gauss") return GAUSS;
	if(n == "gauss-legendre") return GAUSS_LEGENDRE;
	if(n == "newton-cotes") return NEWTON_COTES;

	UG_THROW("GetQuadratureType: type '"<<name<<"' not recognized. Options "
	         "are: best, gauss, gauss-legendre, newton-cotes.");
}


////////////////////////////////////////////////////////////////////////////////
// explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

template class QuadratureRuleProvider<0>;
template class QuadratureRuleProvider<1>;
template class QuadratureRuleProvider<2>;
template class QuadratureRuleProvider<3>;

} // namespace ug
