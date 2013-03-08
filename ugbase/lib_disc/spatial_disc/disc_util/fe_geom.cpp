/*
 * fe_geom.cpp
 *
 *  Created on: 04.09.2010
 *      Author: andreasvogel
 */

#include "fe_geom.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// DimFEGeometry
////////////////////////////////////////////////////////////////////////////////

template <int TWorldDim, int TRefDim>
DimFEGeometry<TWorldDim,TRefDim>::
DimFEGeometry() :
	m_roid(ROID_UNKNOWN), m_quadOrder(0),
	m_lfeID(LFEID(LFEID::NONE, LFEID::INVALID)),
	m_vIPLocal(NULL), m_vQuadWeight(NULL)
{}

template <int TWorldDim, int TRefDim>
DimFEGeometry<TWorldDim,TRefDim>::
DimFEGeometry(size_t order, LFEID lfeid) :
	m_roid(ROID_UNKNOWN), m_quadOrder(order), m_lfeID(lfeid),
	m_vIPLocal(NULL), m_vQuadWeight(NULL)
{}

template <int TWorldDim, int TRefDim>
DimFEGeometry<TWorldDim,TRefDim>::
DimFEGeometry(ReferenceObjectID roid, size_t order, LFEID lfeid) :
	m_roid(roid), m_quadOrder(order), m_lfeID(lfeid),
	m_vIPLocal(NULL), m_vQuadWeight(NULL)
{}

template <int TWorldDim, int TRefDim>
void
DimFEGeometry<TWorldDim,TRefDim>::
update_local(ReferenceObjectID roid, LFEID lfeID, size_t orderQuad)
{
//	remember current setting
	m_roid = roid;
	m_lfeID = lfeID;
	m_quadOrder = orderQuad;

//	request for quadrature rule
	try{
	const QuadratureRule<dim>& quadRule
			= QuadratureRuleProvider<dim>::get_rule(roid, orderQuad);

//	copy quad informations
	m_nip = quadRule.size();
	m_vIPLocal = quadRule.points();
	m_vQuadWeight = quadRule.weights();

	}UG_CATCH_THROW("FEGeometry::update: Quadrature Rule error.");

//	resize for number of integration points
	m_vIPGlobal.resize(m_nip);
	m_vJTInv.resize(m_nip);
	m_vDetJ.resize(m_nip);

	m_vvGradGlobal.resize(m_nip);
	m_vvGradLocal.resize(m_nip);
	m_vvShape.resize(m_nip);

//	request for trial space
	try{
	const LocalShapeFunctionSet<dim>& lsfs
		 = LocalShapeFunctionSetProvider::get<dim>(roid, m_lfeID);

//	copy shape infos
	m_nsh = lsfs.num_sh();

//	resize for number of shape functions
	for(size_t ip = 0; ip < m_nip; ++ip)
	{
		m_vvGradGlobal[ip].resize(m_nsh);
		m_vvGradLocal[ip].resize(m_nsh);
		m_vvShape[ip].resize(m_nsh);
	}

//	get all shapes by on call
	for(size_t ip = 0; ip < m_nip; ++ip)
	{
		lsfs.shapes(&(m_vvShape[ip][0]), m_vIPLocal[ip]);
		lsfs.grads(&(m_vvGradLocal[ip][0]), m_vIPLocal[ip]);
	}

	}UG_CATCH_THROW("FEGeometry::update: Shape Function error.");
}

template <int TWorldDim, int TRefDim>
void
DimFEGeometry<TWorldDim,TRefDim>::
update(GeometricObject* pElem, const MathVector<worldDim>* vCorner,
			LFEID lfeID, size_t orderQuad)
{
//	check if same element
	if(pElem == m_pElem) return;
	else m_pElem = pElem;

//	get reference element type
	ReferenceObjectID roid = pElem->reference_object_id();

//	if already prepared for this roid, skip update of local values
	if(m_roid != roid || lfeID != m_lfeID || (int)orderQuad != m_quadOrder)
		update_local(roid, lfeID, orderQuad);

//	get reference element mapping
	try{
	DimReferenceMapping<dim, worldDim>& map
		= ReferenceMappingProvider::get<dim, worldDim>(roid, vCorner);

//	compute global integration points
	map.local_to_global(&(m_vIPGlobal[0]), &(m_vIPLocal[0]), m_nip);

// 	compute transformation inverse and determinate at ip
	map.jacobian_transposed_inverse(&(m_vJTInv[0]), &(m_vDetJ[0]),
	                                &(m_vIPLocal[0]), m_nip);

// 	compute global gradients
	for(size_t ip = 0; ip < m_nip; ++ip)
		for(size_t sh = 0; sh < m_nsh; ++sh)
			MatVecMult(m_vvGradGlobal[ip][sh],
			           m_vJTInv[ip], m_vvGradLocal[ip][sh]);

	}UG_CATCH_THROW("FEGeometry::update: Reference Mapping error.");
}

////////////////////////////////////////////////////////////////////////////////
// explicit instantiations
////////////////////////////////////////////////////////////////////////////////

template class DimFEGeometry<1, 1>;
template class DimFEGeometry<1, 2>;
template class DimFEGeometry<1, 3>;

template class DimFEGeometry<2, 2>;
template class DimFEGeometry<2, 3>;

template class DimFEGeometry<3, 3>;

} // end namespace ug
