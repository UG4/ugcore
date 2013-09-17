/*
 * hfvcr_geom.cpp
 *
 *  Created on: 17.09.2013
 *      Author: Christian Wehner
 *
 * Node centered finite volume geometry for Crouzeix-Raviart-Elements
 */

#include "common/util/provider.h"
#include "hfvcr_geom.h"
#include "lib_disc/reference_element/reference_element.h"
#include "lib_disc/reference_element/reference_mapping.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/quadrature/quadrature.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//// Methods for HCRFVGeometry class
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <	typename TElem, int TWorldDim>
HCRFVGeometry<TElem, TWorldDim>::
HCRFVGeometry()
	: m_pElem(NULL), m_rRefElem(Provider<ref_elem_type>::get()),
	  m_rTrialSpace(Provider<local_shape_fct_set_type>::get())
{
	update_local_data();
}

template <	typename TElem, int TWorldDim>
void HCRFVGeometry<TElem, TWorldDim>::
update_local_data()
{
//  compute barycenter coordinates
	localBary = m_rRefElem.corner(0);
	for (size_t j=1;j<m_rRefElem.num(0);j++){
	   	localBary+=m_rRefElem.corner(j);
	}
	localBary*=1./(number)m_rRefElem.num(0);

// 	set up local informations for SubControlVolumeFaces (scvf)
// 	each scvf is associated to one edge of the element
	for(size_t i = 0; i < numNaturalSCVF; ++i)
	{
	//	this scvf separates the given edges/faces
		m_vSCVF[i].From = m_rRefElem.id(dim-2, i, dim-1, 0);// todo handle dim==1
		m_vSCVF[i].To = m_rRefElem.id(dim-2, i, dim-1, 1);

		for (size_t j=0;j<m_vSCVF[i].numCo-1;j++){
			m_vSCVF[i].vLocPos[j]=m_rRefElem.corner(m_rRefElem.id(dim-2,i,0,j));
		}

		m_vSCVF[i].vLocPos[m_vSCVF[i].numCo-1]=localBary;

		AveragePositions(m_vSCVF[i].localIP, m_vSCVF[i].vLocPos, m_vSCVF[i].numCo);
	}

// 	set up local informations for SubControlVolumes (scv)
// 	each scv is associated to one corner of the element
	for(size_t i = 0; i < numNaturalSCV; ++i)
	{
	//	store associated node
		m_vSCV[i].nodeID = i;

		m_vSCV[i].numCorners = m_rRefElem.num(dim-1,i,0)+1;
		for (int j=0;j<m_vSCV[i].numCorners-1;j++){
			m_vSCV[i].vLocPos[m_vSCV[i].numCorners-2-j]=m_rRefElem.corner(m_rRefElem.id(dim-1,i,0,j));
		}
		AveragePositions(m_vLocUnkCoords[i], m_vSCV[i].vLocPos, m_vSCV[i].numCorners-1);
		m_vSCV[i].vLocIP=m_vLocUnkCoords[i];

		m_vSCV[i].vLocPos[m_vSCV[i].numCorners-1]=localBary;
	}

	/////////////////////////
	// Shapes and Derivatives
	/////////////////////////

// 	compute Shapes and Derivatives
	for(size_t i = 0; i < numNaturalSCVF; ++i)
	{
		m_rTrialSpace.shapes(&(m_vSCVF[i].vShape[0]), m_vSCVF[i].local_ip());
		m_rTrialSpace.grads(&(m_vSCVF[i].vLocalGrad[0]), m_vSCVF[i].local_ip());
	}

	for(size_t i = 0; i < numNaturalSCV; ++i)
	{
		m_rTrialSpace.shapes(&(m_vSCV[i].vShape[0]), m_vSCV[i].local_ip());
		m_rTrialSpace.grads(&(m_vSCV[i].vLocalGrad[0]), m_vSCV[i].local_ip());
	}

// 	copy ip positions in a list for Sub Control Volumes Faces (SCVF)
	for(size_t i = 0; i < numNaturalSCVF; ++i)
		m_vLocSCVF_IP[i] = scvf(i).local_ip();

	localUpdateNecessary=false;
	numConstrDofs = 0;
	numSCV = numNaturalSCV;
	numSCVF = numNaturalSCVF;
}


/// update data for given element
template <	typename TElem, int TWorldDim>
void HCRFVGeometry<TElem, TWorldDim>::
update(GeometricObject* elem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
	UG_ASSERT(dynamic_cast<TElem*>(elem) != NULL, "Wrong element type.");
	TElem* pElem = static_cast<TElem*>(elem);

// 	if already update for this element, do nothing
	if(m_pElem == pElem) return; else m_pElem = pElem;

//  update local data if some of it has been overwritten by constrained object scv/scvf
	if (localUpdateNecessary) update_local_data();

	//  compute barycenter coordinates
	globalBary = vCornerCoords[0];
	m_vCo[0] = vCornerCoords[0];
	for (size_t j=1;j<m_rRefElem.num(0);j++){
	   globalBary+=vCornerCoords[j];
	   m_vCo[j] = vCornerCoords[j];
	}
	globalBary*=1./(number)m_rRefElem.num(0);

// 	compute global informations for scvf
	for(size_t i = 0; i < num_scvf(); ++i)
	{
		for (size_t j=0;j<m_vSCVF[i].numCo-1;j++){
			m_vSCVF[i].vGloPos[j]=vCornerCoords[m_rRefElem.id(dim-2,i,0,j)];
		}
		m_vSCVF[i].vGloPos[m_vSCVF[i].numCo-1]=globalBary;
		AveragePositions(m_vSCVF[i].globalIP, m_vSCVF[i].vGloPos, m_vSCVF[i].numCo);
		ElementNormal<face_type0,worldDim>(m_vSCVF[i].Normal,m_vSCVF[i].vGloPos);// face_type0 identical to scvf type
	}

// 	compute size of scv
	for(size_t i = 0; i < num_scv(); ++i)
	{
		// side nodes in reverse order to fulfill standard element order
		for (int j=0;j<m_vSCV[i].numCorners-1;j++){
			m_vSCV[i].vGloPos[m_vSCV[i].numCorners-2-j]=vCornerCoords[m_rRefElem.id(dim-1,i,0,j)];
		}
		AveragePositions(m_vGlobUnkCoords[i], m_vSCV[i].vGloPos, m_vSCV[i].numCorners-1);
		m_vSCV[i].vGlobIP = m_vGlobUnkCoords[i];

		m_vSCV[i].vGloPos[m_vSCV[i].numCorners-1]=globalBary;
		// 	compute volume of scv and normal to associated element face
		//CRSCVSizeAndNormal<dim>(m_vSCV[i].Vol,m_vSCV[i].Normal,m_vSCV[i].vGloPos,m_vSCV[i].numCorners);
		if (m_vSCV[i].numCorners-1==dim){
		     m_vSCV[i].Vol = ElementSize<scv_type0,worldDim>(m_vSCV[i].vGloPos);
		     ElementNormal<face_type0, worldDim>(m_vSCV[i].Normal, m_vSCV[i].vGloPos);
		} else { // m_vSCV[i].numCorners-2==dim , only possible in 3d (pyramid)
		     m_vSCV[i].Vol = ElementSize<scv_type1,worldDim>(m_vSCV[i].vGloPos);
		     ElementNormal<face_type1, worldDim>(m_vSCV[i].Normal, m_vSCV[i].vGloPos);
		};
		// nodes are in reverse order therefore reverse sign to get outward normal
		m_vSCV[i].Normal*=-1;
	}

//  check for hanging nodes


// 	Shapes and Derivatives
	m_mapping.update(vCornerCoords);

//	compute jacobian for linear mapping
	if(ReferenceMapping<ref_elem_type, worldDim>::isLinear)
	{
		MathMatrix<worldDim,dim> JtInv;
		m_mapping.jacobian_transposed_inverse(JtInv, m_vSCVF[0].local_ip());
		const number detJ = m_mapping.sqrt_gram_det(m_vSCVF[0].local_ip());

		for(size_t i = 0; i < num_scvf(); ++i)
		{
			m_vSCVF[i].JtInv = JtInv;
			m_vSCVF[i].detj = detJ;
		}

		for(size_t i = 0; i < num_scv(); ++i)
		{
			m_vSCV[i].JtInv = JtInv;
			m_vSCV[i].detj = detJ;
		}
	}
//	else compute jacobian for each integration point
	else
	{
		for(size_t i = 0; i < num_scvf(); ++i)
		{
			m_mapping.jacobian_transposed_inverse(m_vSCVF[i].JtInv, m_vSCVF[i].local_ip());
			m_vSCVF[i].detj = m_mapping.sqrt_gram_det(m_vSCVF[i].local_ip());
		}
		for(size_t i = 0; i < num_scv(); ++i)
		{
			m_mapping.jacobian_transposed_inverse(m_vSCV[i].JtInv, m_vSCV[i].local_ip());
			m_vSCV[i].detj = m_mapping.sqrt_gram_det(m_vSCV[i].local_ip());
		}
	}

//	compute global gradients
	for(size_t i = 0; i < num_scvf(); ++i)
		for(size_t sh = 0; sh < scvf(i).num_sh(); ++sh)
			MatVecMult(m_vSCVF[i].vGlobalGrad[sh], m_vSCVF[i].JtInv, m_vSCVF[i].vLocalGrad[sh]);

	for(size_t i = 0; i < num_scv(); ++i)
		for(size_t sh = 0; sh < scv(i).num_sh(); ++sh)
			MatVecMult(m_vSCV[i].vGlobalGrad[sh], m_vSCV[i].JtInv, m_vSCV[i].vLocalGrad[sh]);

// 	copy ip points in list (SCVF)
	for(size_t i = 0; i < num_scvf(); ++i)
		m_vGlobSCVF_IP[i] = scvf(i).global_ip();

}


#ifdef UG_DIM_2
template class HCRFVGeometry<Triangle, 2>;
template class HCRFVGeometry<Quadrilateral, 2>;
#endif

#ifdef UG_DIM_3
template class HCRFVGeometry<Triangle, 3>;
template class HCRFVGeometry<Quadrilateral, 3>;
template class HCRFVGeometry<Tetrahedron, 3>;
template class HCRFVGeometry<Prism, 3>;
template class HCRFVGeometry<Pyramid, 3>;
template class HCRFVGeometry<Hexahedron, 3>;
#endif

} // end namespace ug
