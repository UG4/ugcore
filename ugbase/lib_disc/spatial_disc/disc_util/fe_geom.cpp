/*
 * fe_geom.cpp
 *
 *  Created on: 04.09.2010
 *      Author: andreasvogel
 */

#include "fe_geom.h"
#include "lib_disc/domain_traits.h"
#include "lib_disc/common/geometry_util.h"
#include "lib_disc/quadrature/quadrature_provider.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// DimFEGeometry
////////////////////////////////////////////////////////////////////////////////

template <int TWorldDim, int TRefDim>
DimFEGeometry<TWorldDim,TRefDim>::
DimFEGeometry() :
	m_roid(ROID_UNKNOWN), m_quadOrder(0),
	m_lfeID(),
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
update_local(ReferenceObjectID roid, const LFEID& lfeID, size_t orderQuad)
{
//	remember current setting
	m_roid = roid;
	m_lfeID = lfeID;
	m_quadOrder = orderQuad;

//	request for quadrature rule
	try{
	const QuadratureRule<dim>& quadRule
			= QuadratureRuleProvider<dim>::get(roid, orderQuad);

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
		 = LocalFiniteElementProvider::get<dim>(roid, m_lfeID);

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
update(GridObject* pElem, const MathVector<worldDim>* vCorner,
		const LFEID& lfeID, size_t orderQuad)
{
//	check if same element
	if(pElem == m_pElem) return;
	else m_pElem = pElem;

//	get reference element type
	ReferenceObjectID roid = pElem->reference_object_id();

//	if already prepared for this roid, skip update of local values
	if(roid != m_roid || lfeID != m_lfeID || (int)orderQuad != m_quadOrder)
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

template <int TWorldDim, int TRefDim>
void
DimFEGeometry<TWorldDim,TRefDim>::
update_boundary_faces(GridObject* pElem,
                      const MathVector<worldDim>* vCornerCoords,
                      size_t quadOrder,
                      const ISubsetHandler* ish)
{
	typedef typename domain_traits<dim>::side_type Side;

	//	get reference element type
	const ReferenceObjectID roid = pElem->reference_object_id();

	//	get grid
	Grid& grid = *(ish->grid());

	//	vector of subset indices of side
	std::vector<int> vSubsetIndex;

	//	get subset indices for sides (i.e. edge in 2d, faces in 3d)
	std::vector<Side*> vSide;
	CollectAssociated(vSide, grid, pElem);
	vSubsetIndex.resize(vSide.size());
	for(size_t i = 0; i < vSide.size(); ++i)
		vSubsetIndex[i] = ish->get_subset_index(vSide[i]);

	//	get reference element mapping
	try{
		DimReferenceMapping<dim, worldDim>& rMapping
		= ReferenceMappingProvider::get<dim, worldDim>(roid);

		//	update reference mapping
		rMapping.update(vCornerCoords);

		const DimReferenceElement<dim>& rRefElem
			= ReferenceElementProvider::get<dim>(roid);

		const LocalShapeFunctionSet<dim>& rTrialSpace =
			LocalFiniteElementProvider::get<dim>(m_roid, m_lfeID);

		//	loop requested subset
		typename std::map<int, std::vector<BF> >::iterator it;
		for (it=m_mapVectorBF.begin() ; it != m_mapVectorBF.end(); ++it)
		{
			//	get subset index
			const int bndIndex = (*it).first;

			//	get vector of BF for element
			std::vector<BF>& vBF = (*it).second;

			//	clear vector
			vBF.clear();

			//	loop sides of element
			for(size_t side = 0; side < vSubsetIndex.size(); ++side)
			{
				//	skip non boundary sides
				if(vSubsetIndex[side] != bndIndex) continue;

				vBF.resize(vBF.size()+1);
				BF& bf = vBF.back();

				const ReferenceObjectID sideRoid = rRefElem.roid(dim-1,side);

				std::vector<MathVector<worldDim> > vSideCorner(rRefElem.num(dim-1, side, 0));
				std::vector<MathVector<dim> > vLocalSideCorner(rRefElem.num(dim-1, side, 0));
				for(size_t co = 0; co < vSideCorner.size(); ++co){
					vSideCorner[co] = vCornerCoords[rRefElem.id(dim-1, side, 0, co)];
					vLocalSideCorner[co] = rRefElem.corner(rRefElem.id(dim-1, side, 0, co));
				}

				const QuadratureRule<dim-1>& rSideQuadRule
						= QuadratureRuleProvider<dim-1>::get(sideRoid, quadOrder);

				// 	normal on scvf
				ElementNormal<worldDim>(sideRoid, bf.Normal, &vSideCorner[0]);

				//	compute volume
				bf.Vol = VecTwoNorm(bf.Normal);

				//	compute local integration points
				bf.vWeight = rSideQuadRule.weights();
				bf.nip = rSideQuadRule.size();
				bf.vLocalIP.resize(bf.nip);
				bf.vGlobalIP.resize(bf.nip);

				DimReferenceMapping<dim-1, dim>& map
					= ReferenceMappingProvider::get<dim-1, dim>(sideRoid, vLocalSideCorner);
				for(size_t ip = 0; ip < rSideQuadRule.size(); ++ip)
					map.local_to_global(bf.vLocalIP[ip], rSideQuadRule.point(ip));

				//	compute global integration points
				for(size_t ip = 0; ip < bf.num_ip(); ++ip)
					rMapping.local_to_global(bf.vGlobalIP[ip], bf.vLocalIP[ip]);

				bf.vJtInv.resize(bf.nip);
				bf.vDetJ.resize(bf.nip);

				rMapping.jacobian_transposed_inverse(&bf.vJtInv[0], &bf.vLocalIP[0], bf.num_ip());
				DimReferenceMapping<dim-1, worldDim>& map2
					= ReferenceMappingProvider::get<dim-1, worldDim>(sideRoid, vSideCorner);
				map2.sqrt_gram_det(&bf.vDetJ[0], rSideQuadRule.points(), bf.num_ip());

				bf.nsh = rTrialSpace.num_sh();
				bf.vvShape.resize(bf.nip);
				bf.vvLocalGrad.resize(bf.nip);
				bf.vvGlobalGrad.resize(bf.nip);
				for(size_t ip = 0; ip < bf.num_ip(); ++ip)
				{
					bf.vvShape[ip].resize(bf.nsh);
					bf.vvLocalGrad[ip].resize(bf.nsh);
					bf.vvGlobalGrad[ip].resize(bf.nsh);
				}

			//	compute shapes and gradients
				for(size_t ip = 0; ip < bf.num_ip(); ++ip)
				{
					rTrialSpace.shapes(&(bf.vvShape[ip][0]), bf.local_ip(ip));
					rTrialSpace.grads(&(bf.vvLocalGrad[ip][0]), bf.local_ip(ip));
				}

			//	compute global gradient
				for(size_t ip = 0; ip < bf.num_ip(); ++ip)
					for(size_t sh = 0 ; sh < bf.num_sh(); ++sh)
						MatVecMult(bf.vvGlobalGrad[ip][sh],
						           bf.vJtInv[ip], bf.vvLocalGrad[ip][sh]);

			}
		}
	}
	UG_CATCH_THROW("DimFEGeometry: update failed.");

}
////////////////////////////////////////////////////////////////////////////////
// explicit instantiations
////////////////////////////////////////////////////////////////////////////////

template class DimFEGeometry<1, 1>;
template class DimFEGeometry<2, 1>;
template class DimFEGeometry<3, 1>;

template class DimFEGeometry<2, 2>;
template class DimFEGeometry<3, 2>;

template class DimFEGeometry<3, 3>;

} // end namespace ug
