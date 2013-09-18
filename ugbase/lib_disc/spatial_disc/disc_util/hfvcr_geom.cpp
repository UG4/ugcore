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
	numConstrainedDofs = 0;
	numSCV = numNaturalSCV;
	numSCVF = numNaturalSCVF;
	numDofs = numSCV;
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
	
// get grid
	Grid& grid = *(ish->grid());

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
	if (dim==2){
		std::vector<EdgeBase*> vEdges;
		CollectEdgesSorted(vEdges, grid, pElem);
		for(size_t side = 0; side < vEdges.size(); ++side){
			ConstrainingEdge* constrainingObj = dynamic_cast<ConstrainingEdge*>(vEdges[side]);
			if(constrainingObj == NULL) continue;
			
			// found constraining edge
			MathVector<worldDim> globalMidpoint = m_vSCV[side].vGlobIP;
			MathVector<dim> localMidpoint = m_vSCV[side].vLocIP;
			// get edge corners
			size_t edgeCo[2];
			for (size_t j=0;j<2;j++) edgeCo[j] = m_rRefElem.id(1,side,0,j);
			// compute dof positions on constraining edge,
			// replace dof "side" with first and insert second at the end
			for (int d=0;d<worldDim;d++)
				m_vGlobUnkCoords[side][d] = 0.5 * (vCornerCoords[edgeCo[0]][d] + globalMidpoint[d]);
			for (int d=0;d<dim;d++)
				m_vLocUnkCoords[side][d] = 0.5 * (m_rRefElem.corner(edgeCo[0])[d] + localMidpoint[d]);
			for (int d=0;d<worldDim;d++)
				m_vGlobUnkCoords[numSCV][d] = 0.5 * (vCornerCoords[edgeCo[1]][d] + globalMidpoint[d]);
			for (int d=0;d<dim;d++)
				m_vLocUnkCoords[numSCV][d] = 0.5 * (m_rRefElem.corner(edgeCo[1])[d] + localMidpoint[d]);
			// handle corresponding scvfs
			for (size_t j=0;j<2;j++){
				if (m_vSCVF[edgeCo[j]].From==side) m_vSCVF[edgeCo[j]].From=numDofs+j;
				if (m_vSCVF[edgeCo[j]].To==side) m_vSCVF[edgeCo[j]].To=numDofs+j;
			}
			// set up new scvs
			// scv "side" gets replaced by new one
			m_vSCV[side].Normal *= 0.5;
			m_vSCV[side].Vol *= 0.5;
			m_vSCV[side].nodeID = numDofs;
			m_vSCV[side].vGlobIP = m_vGlobUnkCoords[side];
			m_vSCV[side].vLocIP = m_vLocUnkCoords[side];
			m_rTrialSpace.shapes(&(m_vSCV[side].vShape[0]), m_vSCV[side].local_ip());
			m_rTrialSpace.grads(&(m_vSCV[side].vLocalGrad[0]), m_vSCV[side].local_ip());
			// second scv inserted at the end
			m_vSCV[numSCV].Normal = m_vSCV[side].Normal;
			m_vSCV[numSCV].Vol = m_vSCV[side].Vol;
			m_vSCV[numSCV].nodeID = numDofs+1;
			m_vSCV[numSCV].vGlobIP = m_vGlobUnkCoords[numSCV];
			m_vSCV[numSCV].vLocIP = m_vLocUnkCoords[numSCV];
			m_rTrialSpace.shapes(&(m_vSCV[numSCV].vShape[0]), m_vSCV[numSCV].local_ip());
			m_rTrialSpace.grads(&(m_vSCV[numSCV].vLocalGrad[0]), m_vSCV[numSCV].local_ip());
			// insert new scvf
			m_vSCVF[numSCVF].From = numDofs;
			m_vSCVF[numSCVF].To = numDofs+1;
			m_vSCVF[numSCVF].vLocPos[0] = localMidpoint;
			m_vSCVF[numSCVF].vLocPos[1] = localBary;
			m_vSCVF[numSCVF].vGloPos[0] = globalMidpoint;
			m_vSCVF[numSCVF].vGloPos[1] = globalBary;
			for (int d=0;d<dim;d++) m_vSCVF[numSCVF].localIP[d] = 0.5*(localMidpoint[d] + localBary[d]);
			for (int d=0;d<worldDim;d++) m_vSCVF[numSCVF].globalIP[d] = 0.5*(globalMidpoint[d] + globalBary[d]);
			m_rTrialSpace.shapes(&(m_vSCVF[numSCVF].vShape[0]), m_vSCVF[numSCVF].local_ip());
			m_rTrialSpace.grads(&(m_vSCVF[numSCVF].vLocalGrad[0]), m_vSCVF[numSCVF].local_ip());
			ElementNormal<face_type0,worldDim>(m_vSCVF[numSCVF].Normal,m_vSCVF[numSCVF].vGloPos);
			// insert new constrained dof object
			m_vCD[numConstrainedDofs].i = side;
			m_vCD[numConstrainedDofs].numConstrainingDofs = 2;
			m_vCD[numConstrainedDofs].cDofInd[0] = numDofs;
			m_vCD[numConstrainedDofs].cDofInd[1] = numDofs+1;
			m_vCD[numConstrainedDofs].cDofWeights[0] = 0.5;
			m_vCD[numConstrainedDofs].cDofWeights[1] = 0.5;
			numSCV+=1;
			numSCVF+=1;
			numDofs+=2;
			numConstrainedDofs+=1;
			localUpdateNecessary = true;
		}
	}


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

//  debug output
//	if (localUpdateNecessary==true) print();

}


// debug output
template <	typename TElem, int TWorldDim>
void HCRFVGeometry<TElem, TWorldDim>::print()
{
	UG_LOG("\nFVG hanging debug output\n");
	for(size_t i = 0; i < numSCV; ++i)
	{
		UG_LOG(i<<" SCV: ");
		UG_LOG("node_id=" << m_vSCV[i].node_id());
		UG_LOG(", local_pos="<< m_vSCV[i].local_ip());
		UG_LOG(", global_pos="<< m_vSCV[i].global_ip());
		UG_LOG(", vol=" << m_vSCV[i].volume());
	//	UG_LOG("\n    localCorner=" << m_vSCV[i].m_vLocPos[0]);
	//	UG_LOG(", localSide1=" << m_vSCV[i].m_vLocPos[1]);
	//	UG_LOG(", localCenter=" << m_vSCV[i].m_vLocPos[2]);
	//	UG_LOG(", localSide2=" << m_vSCV[i].m_vLocPos[3]);
	//	UG_LOG("\n    globalCorner=" << m_vSCV[i].m_vGloPos[0]);
	//	UG_LOG(", globalSide1=" << m_vSCV[i].m_vGloPos[1]);
	//	UG_LOG(", globalCenter=" << m_vSCV[i].m_vGloPos[2]);
	//	UG_LOG(", globalSide2=" << m_vSCV[i].m_vGloPos[3]);

		UG_LOG("\n");
	}
	UG_LOG("\n");
	for(size_t i = 0; i < numSCVF; ++i)
	{
		UG_LOG(i<<" SCVF: ");
		UG_LOG("from=" << m_vSCVF[i].from()<<", to="<<m_vSCVF[i].to());
		UG_LOG(", local_pos="<< m_vSCVF[i].local_ip());
		UG_LOG(", global_pos="<< m_vSCVF[i].global_ip());
		UG_LOG(", normal=" << m_vSCVF[i].normal());
		UG_LOG("\n    Shapes:\n");
		for(size_t sh=0; sh < m_vSCVF[i].num_sh(); ++sh)
		{
			UG_LOG("         " <<sh << ": shape="<<m_vSCVF[i].shape(sh));
			UG_LOG(", global_grad="<<m_vSCVF[i].global_grad(sh));
			UG_LOG(", local_grad="<<m_vSCVF[i].local_grad(sh));
			UG_LOG("\n");
		}
	}
	UG_LOG("\n");
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
