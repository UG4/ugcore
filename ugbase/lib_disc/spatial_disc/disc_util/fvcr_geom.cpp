/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Christian Wehner
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/*
 * Node centered finite volume geometry for Crouzeix-Raviart-Elements
 */

#include "common/util/provider.h"
#include "fvcr_geom.h"
#include "lib_disc/reference_element/reference_element.h"
#include "lib_disc/reference_element/reference_mapping.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/quadrature/quadrature.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Dim-dependent Crouzeix-Raviart Finite Volume Geometry
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int TDim, int TWorldDim>
void DimCRFVGeometry<TDim, TWorldDim>::
update_local_data()
{
//	get reference element
	try{
	m_rRefElem
		= ReferenceElementProvider::get<dim>(m_roid);

//	set number of scvf / scv of this roid
	m_numSCV = m_rRefElem.num(dim-1); // number of faces
	m_numSCVF = m_rRefElem.num(1); // number of edges

//  compute barycenter coordinates
	localBary = m_rRefElem.corner(0);
	for (size_t j=1;j<m_rRefElem.num(0);j++){
	   	localBary+=m_rRefElem.corner(j);
	}
	localBary*=1./(number)m_rRefElem.num(0);

// 	set up local informations for SubControlVolumeFaces (scvf)
// 	each scvf is associated to one vertex (2d) / edge (3d) of the element
	for(size_t i = 0; i < m_numSCVF; ++i)
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
// 	each scv is associated to one edge(2d) / face(3d) of the element
	for(size_t i = 0; i < m_numSCV; ++i)
	{
	//	store associated node
		m_vSCV[i].nodeID = i;

		m_vSCV[i].numCorners = m_rRefElem.num(dim-1,i,0)+1;
		for (int j=0;j<m_vSCV[i].numCorners-1;j++){
			m_vSCV[i].vLocPos[m_vSCV[i].numCorners-2-j]=m_rRefElem.corner(m_rRefElem.id(dim-1,i,0,j));
		}
		AveragePositions(m_vLocUnkCoords[i], m_vSCV[i].vLocPos, m_vSCV[i].numCorners-1);
		m_vSCV[i].vLocIP = m_vLocUnkCoords[i];
		m_vSCV[i].vLocPos[m_vSCV[i].numCorners-1]=localBary;
	}

	/////////////////////////
	// Shapes and Derivatives
	/////////////////////////

	const LocalShapeFunctionSet<dim>& rTrialSpace =
		LocalFiniteElementProvider::get<dim>(m_roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

	m_nsh = rTrialSpace.num_sh();

	for(size_t i = 0; i < m_numSCVF; ++i)
	{
		m_vSCVF[i].numSH = rTrialSpace.num_sh();
		rTrialSpace.shapes(&(m_vSCVF[i].vShape[0]), m_vSCVF[i].local_ip());
		rTrialSpace.grads(&(m_vSCVF[i].vLocalGrad[0]), m_vSCVF[i].local_ip());
	}

	for(size_t i = 0; i < m_numSCV; ++i)
	{
		m_vSCV[i].numSH = rTrialSpace.num_sh();
		rTrialSpace.shapes(&(m_vSCV[i].vShape[0]), m_vSCV[i].local_ip());
		rTrialSpace.grads(&(m_vSCV[i].vLocalGrad[0]), m_vSCV[i].local_ip());
	}

	}
	UG_CATCH_THROW("DimCRFVGeometry: update failed.");

// 	copy ip positions in a list for Sub Control Volumes Faces (SCVF)
	for(size_t i = 0; i < m_numSCVF; ++i)
		m_vLocSCVF_IP[i] = scvf(i).local_ip();
		
	m_numConstrainedDofs = 0;
	m_numConstrainedSCVF = 0;
}


/// update data for given element
template <int TDim, int TWorldDim>
void DimCRFVGeometry<TDim, TWorldDim>::
update(GridObject* pElem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
// 	If already update for this element, do nothing
	if(m_pElem == pElem) return; else m_pElem = pElem;

//	refresh local data, if different roid given
	if(m_roid != pElem->reference_object_id())
	{
	//	remember new roid
		m_roid = (ReferenceObjectID) pElem->reference_object_id();

	//	update local data
		update_local_data();
	}

//	get reference element
	try{
	const DimReferenceElement<dim>& m_rRefElem
		= ReferenceElementProvider::get<dim>(m_roid);

	//  compute barycenter coordinates
	globalBary = vCornerCoords[0];
	for (size_t j=1;j<m_rRefElem.num(0);j++){
	   globalBary+=vCornerCoords[j];
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
		     ElementNormal<face_type1,worldDim>(m_vSCV[i].Normal, m_vSCV[i].vGloPos);
		};
		// nodes are in reverse order therefore reverse sign to get outward normal
		m_vSCV[i].Normal*=-1;
	}

//	get reference mapping
	DimReferenceMapping<dim, worldDim>& rMapping = ReferenceMappingProvider::get<dim, worldDim>(m_roid);
	rMapping.update(vCornerCoords);

	//\todo compute with on virt. call
//	compute jacobian for linear mapping
	if(rMapping.is_linear())
	{
		MathMatrix<worldDim,dim> JtInv;
		rMapping.jacobian_transposed_inverse(JtInv, m_vSCVF[0].local_ip());
		const number detJ = rMapping.sqrt_gram_det(m_vSCVF[0].local_ip());

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
			rMapping.jacobian_transposed_inverse(m_vSCVF[i].JtInv, m_vSCVF[i].local_ip());
			m_vSCVF[i].detj = rMapping.sqrt_gram_det(m_vSCVF[i].local_ip());
		}
		for(size_t i = 0; i < num_scv(); ++i)
		{
			rMapping.jacobian_transposed_inverse(m_vSCV[i].JtInv, m_vSCV[i].local_ip());
			m_vSCV[i].detj = rMapping.sqrt_gram_det(m_vSCV[i].local_ip());
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
	UG_CATCH_THROW("DimCRFVGeometry: update failed.");

//	if no boundary subsets required, return
	if(num_boundary_subsets() == 0 || ish == NULL) return;
	else update_boundary_faces(pElem, vCornerCoords, ish);
}

/// update data checking for hanging nodes for given element
template <int TDim, int TWorldDim>
void DimCRFVGeometry<TDim, TWorldDim>::
update_hanging(GridObject* pElem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish,bool keepSCV,bool keepSCVF)
{
// 	If already update for this element, do nothing
	if(m_pElem == pElem) return; else m_pElem = pElem;
	
//  get grid
	Grid& grid = *(ish->grid());

//	refresh local data, if different roid given
	if(m_roid != pElem->reference_object_id())
	{
	//	remember new roid
		m_roid = (ReferenceObjectID) pElem->reference_object_id();

	//	update local data
		update_local_data();
		
		m_numDofs	= num_scv();
	}
	
	const LocalShapeFunctionSet<dim>& rTrialSpace =
				LocalFiniteElementProvider::get<dim>(m_roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

	//  compute barycenter coordinates
	globalBary = vCornerCoords[0];
	for (size_t j=1;j<m_rRefElem.num(0);j++){
	   globalBary+=vCornerCoords[j];
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
		     ElementNormal<face_type1,worldDim>(m_vSCV[i].Normal, m_vSCV[i].vGloPos);
		};
		// nodes are in reverse order therefore reverse sign to get outward normal
		m_vSCV[i].Normal*=-1;
	}
	
	//	get reference mapping
	DimReferenceMapping<dim, worldDim>& rMapping = ReferenceMappingProvider::get<dim, worldDim>(m_roid);
	rMapping.update(vCornerCoords);
	
	//  check for hanging nodes
	if (dim==2){
		std::vector<Edge*> vEdges;
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
			// set up new scvs
			// keepSCV parameter specifies if scv "side" gets replaced by new one
			size_t ind=side;
			size_t keepOffset=0;
			MathVector<worldDim> normal = m_vSCV[side].Normal;
			normal*=0.5;
			number vol    = 0.5*m_vSCV[side].Vol;
 			if (keepSCV){ 
				ind=m_numSCV+1;
				keepOffset=1;
			}
			for (int d=0;d<worldDim;d++)
				m_vGlobUnkCoords[ind][d] = 0.5 * (vCornerCoords[edgeCo[0]][d] + globalMidpoint[d]);
			for (int d=0;d<dim;d++)
				m_vLocUnkCoords[ind][d] = 0.5 * (m_rRefElem.corner(edgeCo[0])[d] + localMidpoint[d]);
			for (int d=0;d<worldDim;d++)
				m_vGlobUnkCoords[m_numSCV][d] = 0.5 * (vCornerCoords[edgeCo[1]][d] + globalMidpoint[d]);
			for (int d=0;d<dim;d++)
				m_vLocUnkCoords[m_numSCV][d] = 0.5 * (m_rRefElem.corner(edgeCo[1])[d] + localMidpoint[d]);
			// handle corresponding scvfs
			// if keepSCVF copy them into constrained scvf array
			if (keepSCVF){
				for (size_t j=0;j<2;j++){
					m_vConstrainedSCVF[m_numConstrainedSCVF+j] = m_vSCVF[edgeCo[j]];
					if (m_vConstrainedSCVF[edgeCo[j]].To==side){
						m_vConstrainedSCVF[edgeCo[j]].From=side;
						m_vConstrainedSCVF[edgeCo[j]].Normal*=-1;
					}
				}
				m_numConstrainedSCVF += 2;
			}
			for (size_t j=0;j<2;j++){
				if (m_vSCVF[edgeCo[j]].From==side) m_vSCVF[edgeCo[j]].From=m_numDofs+j;
				if (m_vSCVF[edgeCo[j]].To==side) m_vSCVF[edgeCo[j]].To=m_numDofs+j;
			}
			m_vSCV[ind].Normal = normal;
			m_vSCV[ind].Vol = vol;
			m_vSCV[ind].nodeID = m_numDofs;
			m_vSCV[ind].vGlobIP = m_vGlobUnkCoords[ind];
			m_vSCV[ind].vLocIP = m_vLocUnkCoords[ind];
			m_vSCV[ind].numSH = rTrialSpace.num_sh();
			m_vSCV[ind].vGloPos[0]=vCornerCoords[edgeCo[0]];
			m_vSCV[ind].vGloPos[1]=globalMidpoint;
			m_vSCV[ind].vGloPos[2]=globalBary;
			m_vSCV[ind].vLocPos[0]=m_rRefElem.corner(edgeCo[0]);
			m_vSCV[ind].vLocPos[1]=localMidpoint;
			m_vSCV[ind].vLocPos[2]=localBary;
			m_vSCV[ind].numCorners = 3;
			rTrialSpace.shapes(&(m_vSCV[ind].vShape[0]), m_vSCV[ind].local_ip());
			rTrialSpace.grads(&(m_vSCV[ind].vLocalGrad[0]), m_vSCV[ind].local_ip());
			// second scv inserted at the end
			m_vSCV[m_numSCV].Normal = normal;
			m_vSCV[m_numSCV].Vol = vol;
			m_vSCV[m_numSCV].nodeID = m_numDofs+1;
			m_vSCV[m_numSCV].vGlobIP = m_vGlobUnkCoords[m_numSCV];
			m_vSCV[m_numSCV].vLocIP = m_vLocUnkCoords[m_numSCV];
			m_vSCV[m_numSCV].numSH = rTrialSpace.num_sh();
			m_vSCV[m_numSCV].vGloPos[0]=vCornerCoords[edgeCo[1]];
			m_vSCV[m_numSCV].vGloPos[1]=globalMidpoint;
			m_vSCV[m_numSCV].vGloPos[2]=globalBary;
			m_vSCV[m_numSCV].vLocPos[0]=m_rRefElem.corner(edgeCo[1]);
			m_vSCV[m_numSCV].vLocPos[1]=localMidpoint;
			m_vSCV[m_numSCV].vLocPos[2]=localBary;
			m_vSCV[m_numSCV].numCorners = 3;
			rTrialSpace.shapes(&(m_vSCV[m_numSCV].vShape[0]), m_vSCV[m_numSCV].local_ip());
			rTrialSpace.grads(&(m_vSCV[m_numSCV].vLocalGrad[0]), m_vSCV[m_numSCV].local_ip());
			// insert new scvf
			m_vSCVF[m_numSCVF].From = m_numDofs;
			m_vSCVF[m_numSCVF].To = m_numDofs+1;
			m_vSCVF[m_numSCVF].vLocPos[0] = localMidpoint;
			m_vSCVF[m_numSCVF].vLocPos[1] = localBary;
			m_vSCVF[m_numSCVF].vGloPos[0] = globalMidpoint;
			m_vSCVF[m_numSCVF].vGloPos[1] = globalBary;
			m_vSCVF[m_numSCVF].numSH = rTrialSpace.num_sh();
			for (int d=0;d<dim;d++) m_vSCVF[m_numSCVF].localIP[d] = 0.5*(localMidpoint[d] + localBary[d]);
			for (int d=0;d<worldDim;d++) m_vSCVF[m_numSCVF].globalIP[d] = 0.5*(globalMidpoint[d] + globalBary[d]);
			rTrialSpace.shapes(&(m_vSCVF[m_numSCVF].vShape[0]), m_vSCVF[m_numSCVF].local_ip());
			rTrialSpace.grads(&(m_vSCVF[m_numSCVF].vLocalGrad[0]), m_vSCVF[m_numSCVF].local_ip());
			ElementNormal<face_type0,worldDim>(m_vSCVF[m_numSCVF].Normal,m_vSCVF[m_numSCVF].vGloPos);
			// insert new constrained dof object
			m_vCD[m_numConstrainedDofs].i = side;
			m_vCD[m_numConstrainedDofs].numConstrainingDofs = 2;
			m_vCD[m_numConstrainedDofs].cDofInd[0] = m_numDofs;
			m_vCD[m_numConstrainedDofs].cDofInd[1] = m_numDofs+1;
			m_vCD[m_numConstrainedDofs].cDofWeights[0] = 0.5;
			m_vCD[m_numConstrainedDofs].cDofWeights[1] = 0.5;
			m_numSCV+=1+keepOffset;
			m_numSCVF+=1;
			m_numDofs+=2;
			m_numConstrainedDofs+=1;
			m_roid = ROID_UNKNOWN;
		}
	} else {
		// dim == 3
		std::vector<Face*> vFaces;
		CollectFacesSorted(vFaces, grid, pElem);
		handledEdges.clear();
		for(size_t face = 0; face < vFaces.size(); ++face){
			ConstrainingFace* constrainingObj = dynamic_cast<ConstrainingFace*>(vFaces[face]);
			if(constrainingObj == NULL) continue;
			// found constraining face
			MathVector<worldDim> globalMidpoint = m_vSCV[face].vGlobIP;
			MathVector<dim> localMidpoint = m_vSCV[face].vLocIP;
			number faceVol = m_vSCV[face].Vol;
			MathVector<worldDim> faceNormal = m_vSCV[face].Normal;
			// get face corners and edges
			size_t faceCo[4];
			size_t faceEdge[4];
			// number of corners of face = number of edges
			size_t numFaceCo = m_rRefElem.num(2,face,0);
			for (size_t j=0;j<numFaceCo;j++) faceCo[j] = m_rRefElem.id(2,face,0,j);
			for (size_t j=0;j<numFaceCo;j++) faceEdge[j] = m_rRefElem.id(2,face,1,j);
			number volSum=0;
			size_t keepOffset=0;
			if (keepSCV) keepOffset=1;
			// compute coordinates of each face and fill scv values
			for (size_t i=0;i<numFaceCo;i++){
				size_t co = faceCo[i];
				size_t nOfEdges=0;
				size_t nbEdges[2];
				// find 2 edges in face belonging to node
				for (size_t j=0;j<m_rRefElem.num(0,co,1);j++){
					size_t candidate = m_rRefElem.id(0,co,1,j);
					bool found = false;
					for (size_t k=0;k<numFaceCo;k++){
						if (faceEdge[k]==candidate){
							found = true;
							break;
						}
					}
					if (found==true){
						nbEdges[nOfEdges] = candidate;
						nOfEdges++;
						if (nOfEdges==2) break;
					}
				}
				// in triangular case switch edges if necessary for correct orientation
				if (numFaceCo==3){
					if (faceEdge[i]==nbEdges[1]){
						nbEdges[1] = nbEdges[0];
						nbEdges[0] = faceEdge[i];
					}
				}
				// keepSCV parameter specifies if scv "side" gets replaced by new one
				size_t ind = m_numSCV+i-1+keepOffset;
				if (i==0){
					if (keepSCV) ind=m_numSCV;
					else ind = face;
				}
				// others are inserted at the end
				m_vSCV[ind].vGloPos[0] = vCornerCoords[co];
				m_vSCV[ind].vLocPos[0] = m_rRefElem.corner(co);
				for (int d=0;d<worldDim;d++){
 					 // edge 0 midpoint
  					m_vSCV[ind].vGloPos[1][d] = 0.5 * ( vCornerCoords[m_rRefElem.id(1,nbEdges[0],0,0)][d] + vCornerCoords[m_rRefElem.id(1,nbEdges[0],0,1)][d] );
					m_vSCV[ind].vLocPos[1][d] = 0.5 * ( m_rRefElem.corner(m_rRefElem.id(1,nbEdges[0],0,0))[d] + m_rRefElem.corner(m_rRefElem.id(1,nbEdges[0],0,1))[d] );
					// edge 1 midpoint
					m_vSCV[ind].vGloPos[numFaceCo-1][d] = 0.5 * ( vCornerCoords[m_rRefElem.id(1,nbEdges[1],0,0)][d] + vCornerCoords[m_rRefElem.id(1,nbEdges[1],0,1)][d] );
					m_vSCV[ind].vLocPos[numFaceCo-1][d] = 0.5 * ( m_rRefElem.corner(m_rRefElem.id(1,nbEdges[1],0,0))[d] + m_rRefElem.corner(m_rRefElem.id(1,nbEdges[1],0,1))[d] );
				}
				if (numFaceCo==4) m_vSCV[ind].vGloPos[2] = globalMidpoint;
				m_vSCV[ind].vGloPos[numFaceCo] = globalBary;
				m_vSCV[ind].vLocPos[numFaceCo] = localBary;
				m_vSCV[ind].numCorners = numFaceCo + 1;
				AveragePositions(m_vGlobUnkCoords[ind], m_vSCV[ind].vGloPos, m_vSCV[ind].numCorners-1);
				AveragePositions(m_vLocUnkCoords[ind], m_vSCV[ind].vLocPos, m_vSCV[ind].numCorners-1);
				m_vSCV[ind].vLocIP = m_vLocUnkCoords[ind];
				m_vSCV[ind].vGlobIP = m_vGlobUnkCoords[ind];
				m_vSCV[ind].numSH = rTrialSpace.num_sh();
				if (numFaceCo==3) m_vSCV[ind].Vol = ElementSize<scv_type0,worldDim>(m_vSCV[ind].vGloPos);
				else m_vSCV[ind].Vol = ElementSize<scv_type1,worldDim>(m_vSCV[ind].vGloPos);
				if (m_vSCV[ind].Vol<0) m_vSCV[ind].Vol *= -1;
				volSum+=m_vSCV[ind].Vol;
				m_vSCV[ind].Normal = faceNormal;
				m_vSCV[ind].Normal *= (number) m_vSCV[ind].Vol / faceVol;
				rTrialSpace.shapes(&(m_vSCV[ind].vShape[0]), m_vSCV[ind].local_ip());
				rTrialSpace.grads(&(m_vSCV[ind].vLocalGrad[0]), m_vSCV[ind].local_ip());
				m_vSCV[ind].nodeID = m_numDofs+i;
			}
			// compute inner scv in triangular case
			if (numFaceCo==3){
				size_t ind = m_numSCV+2+keepOffset;
				m_vSCV[ind].Vol = faceVol - volSum;
				m_vSCV[ind].nodeID=m_numDofs+3;
				m_vSCV[ind].Normal = faceNormal;
				m_vSCV[ind].Normal *= (number) m_vSCV[ind].Vol / faceVol;
				m_vSCV[ind].vGlobIP = m_vSCV[face].vGloPos[1];
				m_vSCV[ind].vLocIP = m_vSCV[face].vLocPos[1];
				for (size_t j=0;j<2;j++){
					m_vSCV[ind].vGlobIP += m_vSCV[m_numSCV+j].vGloPos[1];
					m_vSCV[ind].vLocIP += m_vSCV[m_numSCV+j].vLocPos[1];
				}
				m_vSCV[ind].vGlobIP *= (number)1.0/3.0;
				m_vSCV[ind].vLocIP *= (number)1.0/3.0;
				m_vGlobUnkCoords[ind] = m_vSCV[ind].vGlobIP;
				m_vLocUnkCoords[ind] = m_vSCV[ind].vLocIP;
				m_vSCV[ind].numCorners = numFaceCo + 1;
				m_vSCV[ind].numSH = rTrialSpace.num_sh();
				rTrialSpace.shapes(&(m_vSCV[ind].vShape[0]), m_vSCV[ind].local_ip());
				rTrialSpace.grads(&(m_vSCV[ind].vLocalGrad[0]), m_vSCV[ind].local_ip());
			}
			// copy scvfs into constrained scvfs
			if (keepSCVF){
				for (size_t i=0;i<numFaceCo;i++){
					size_t edge = faceEdge[i];
					m_vConstrainedSCVF[m_numConstrainedSCVF+i] = m_vSCVF[edge];
					if (m_vConstrainedSCVF[m_numConstrainedSCVF+i].To==face){
						m_vConstrainedSCVF[m_numConstrainedSCVF+i].From=face;
						m_vConstrainedSCVF[m_numConstrainedSCVF+i].Normal*=-1;
					}
				}
				m_numConstrainedSCVF += numFaceCo;
			}
			// insert new scvfs, first the ones associated to edges of face
			for (size_t i=0;i<numFaceCo;i++){
				size_t edge = faceEdge[i];
				size_t from = m_vSCVF[edge].From;
				size_t to   = m_vSCVF[edge].To;
				MathVector<worldDim> normal = m_vSCVF[edge].Normal;
				normal*=0.5;
				size_t edgeCo[2];
				size_t scvID[2];
				for (size_t j=0;j<2;j++){
					edgeCo[j] = m_rRefElem.id(1,edge,0,j);
					// find corresponding face vertex (= corresponding scv index)
					for (size_t k=0;k<numFaceCo;k++){
						if (faceCo[k]==edgeCo[j]){
							scvID[j] = m_numDofs+k;
							break;
						}
					}
				}
				// look if edge has already been handled
				if (m_numConstrainedDofs>0){
					bool found=false;
					for (size_t j=0;j<handledEdges.size();j++){
						if (handledEdges[j].index==edge){
							HandledEdge& hE=handledEdges[j];
							found=true;
							// set new from/to values
							for (size_t k=0;k<2;k++){
								if (hE.from){
									m_vSCVF[hE.scvfIndex+k].To=scvID[k];	
								} else {
									m_vSCVF[hE.scvfIndex+k].From=scvID[k];
								}
							}
							// mark edge so associated scvf is not overwritten
							faceEdge[i]=deleted;
							break;
						}
					}
					if (found==true) continue;
				}
				HandledEdge hEdge;
				hEdge.index=edge;
				hEdge.scvfIndex=m_numSCVF;
				MathVector<worldDim> edgeMidGlob;
				MathVector<dim> edgeMidLoc;
				for (int d=0;d<worldDim;d++){
					edgeMidGlob[d] = 0.5 * (vCornerCoords[edgeCo[0]][d]  + vCornerCoords[edgeCo[1]][d]);
					edgeMidLoc[d] = 0.5 * (m_rRefElem.corner(edgeCo[0])[d] +  m_rRefElem.corner(edgeCo[1])[d]);
				}
				for (size_t j=0;j<2;j++){
					hEdge.associatedSCV[j] = scvID[j];
					if (from==face){
						m_vSCVF[m_numSCVF].From = scvID[j];
						m_vSCVF[m_numSCVF].To   = to;
						hEdge.from=true;
					} else {
						m_vSCVF[m_numSCVF].From = from;
						m_vSCVF[m_numSCVF].To 	= scvID[j];
						hEdge.from=false;
					}
					m_vSCVF[m_numSCVF].Normal = normal;
					m_vSCVF[m_numSCVF].vGloPos[0] = vCornerCoords[edgeCo[j]];
					m_vSCVF[m_numSCVF].vLocPos[0] = m_rRefElem.corner(edgeCo[j]);
					m_vSCVF[m_numSCVF].vGloPos[1] = edgeMidGlob;
					m_vSCVF[m_numSCVF].vLocPos[1] = edgeMidLoc;
					m_vSCVF[m_numSCVF].vGloPos[2] = globalBary;
					m_vSCVF[m_numSCVF].vLocPos[2] = localBary;
					m_vSCVF[m_numSCVF].numSH = rTrialSpace.num_sh();
					AveragePositions(m_vSCVF[m_numSCVF].localIP, m_vSCVF[m_numSCVF].vLocPos, m_vSCVF[m_numSCVF].numCo);
					AveragePositions(m_vSCVF[m_numSCVF].globalIP, m_vSCVF[m_numSCVF].vGloPos, m_vSCVF[m_numSCVF].numCo);
					rTrialSpace.shapes(&(m_vSCVF[m_numSCVF].vShape[0]), m_vSCVF[m_numSCVF].local_ip());
					rTrialSpace.grads(&(m_vSCVF[m_numSCVF].vLocalGrad[0]), m_vSCVF[m_numSCVF].local_ip());
					m_numSCVF++;
				}
				handledEdges.push_back(hEdge);
			}
			// scvfs inside the face
			// insert remaining inner scvfs into positions of edge associated scvs
			for (size_t j=0;j<numFaceCo;j++){
				// replaces edge associated scvf
				size_t ii = faceEdge[j];
				if (ii==deleted){ 
					ii = m_numSCVF;
					m_numSCVF++;
				}
				if (numFaceCo==3){
					// compute inner scvfs in triangular case
					m_vSCVF[ii].From = m_numDofs+3;
					m_vSCVF[ii].To = m_numDofs+j;
					size_t ind = face;
					if (j>0) ind = m_numSCV+j-1;
					m_vSCVF[ii].vLocPos[0] = m_vSCV[ind].vLocPos[1];
					m_vSCVF[ii].vLocPos[1] = m_vSCV[ind].vLocPos[2];
					m_vSCVF[ii].vGloPos[0] = m_vSCV[ind].vGloPos[1];
					m_vSCVF[ii].vGloPos[1] = m_vSCV[ind].vGloPos[2];
				}else{
					// compute inner scvfs in quadrilateral case
					m_vSCVF[ii].To = m_numDofs+j;
					m_vSCVF[ii].From = m_numDofs + ((j+1) % 4);
					for (int d=0;d<worldDim;d++){
						m_vSCVF[ii].vLocPos[0][d] = 0.5*( m_rRefElem.corner(faceCo[j])[d] + m_rRefElem.corner(faceCo[(j+1) % 4])[d] );
						m_vSCVF[ii].vGloPos[0][d] = 0.5*( vCornerCoords[faceCo[j]][d] + vCornerCoords[faceCo[(j+1) % 4]][d] );
					}
					m_vSCVF[ii].vLocPos[1] = localMidpoint;
					m_vSCVF[ii].vGloPos[1] = globalMidpoint;
				}
				m_vSCVF[ii].vLocPos[2] = localBary;
				m_vSCVF[ii].vGloPos[2] = globalBary;	
				m_vSCVF[ii].numSH = rTrialSpace.num_sh();
				ElementNormal<face_type0,worldDim>(m_vSCVF[ii].Normal,m_vSCVF[ii].vGloPos);
				AveragePositions(m_vSCVF[ii].globalIP,m_vSCVF[ii].vGloPos,3);
				AveragePositions(m_vSCVF[ii].localIP,m_vSCVF[ii].vLocPos,3);
				rTrialSpace.shapes(&(m_vSCVF[ii].vShape[0]), m_vSCVF[ii].local_ip());
				rTrialSpace.grads(&(m_vSCVF[ii].vLocalGrad[0]), m_vSCVF[ii].local_ip());
			}
			// insert new constrained dof object
			m_vCD[m_numConstrainedDofs].i = face;
			m_vCD[m_numConstrainedDofs].numConstrainingDofs = 4;
			for (size_t i=0;i<4;i++){
				m_vCD[m_numConstrainedDofs].cDofInd[i] = m_numDofs+i;
				size_t ind=face;
				if (i>0) ind = m_numSCV+i-1;
				m_vCD[m_numConstrainedDofs].cDofWeights[i] = (number)m_vSCV[ind].Vol / faceVol;
			}
			m_numSCV+=3+keepOffset;
			m_numDofs+=4;
			m_numConstrainedDofs+=1;
			m_roid = ROID_UNKNOWN;
		}
	}// end of hanging node check

	//\todo compute with one virt. call
//	compute jacobian for linear mapping
	if(rMapping.is_linear())
	{
		MathMatrix<worldDim,dim> JtInv;
		rMapping.jacobian_transposed_inverse(JtInv, m_vSCVF[0].local_ip());
		const number detJ = rMapping.sqrt_gram_det(m_vSCVF[0].local_ip());

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
			rMapping.jacobian_transposed_inverse(m_vSCVF[i].JtInv, m_vSCVF[i].local_ip());
			m_vSCVF[i].detj = rMapping.sqrt_gram_det(m_vSCVF[i].local_ip());
		}
		for(size_t i = 0; i < num_scv(); ++i)
		{
			rMapping.jacobian_transposed_inverse(m_vSCV[i].JtInv, m_vSCV[i].local_ip());
			m_vSCV[i].detj = rMapping.sqrt_gram_det(m_vSCV[i].local_ip());
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

//	if no boundary subsets required, return
	if(num_boundary_subsets() == 0 || ish == NULL) return;
	else update_boundary_faces(pElem, vCornerCoords, ish);
}


/// update geometric data for given element (no shapes)
template <int TDim, int TWorldDim>
void DimCRFVGeometry<TDim, TWorldDim>::
update_geometric_data(GridObject* pElem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
// 	If already update for this element, do nothing
	if(m_pElem == pElem) return; else m_pElem = pElem;

//	refresh local data, if different roid given
	if(m_roid != pElem->reference_object_id())
	{
	//	remember new roid
		m_roid = (ReferenceObjectID) pElem->reference_object_id();

	//	update local data
		update_local_data();
	}

//	get reference element
	try{
	const DimReferenceElement<dim>& rRefElem
		= ReferenceElementProvider::get<dim>(m_roid);

	//  compute barycenter coordinates
	globalBary = vCornerCoords[0];
	for (size_t j=1;j<rRefElem.num(0);j++){
	   globalBary+=vCornerCoords[j];
	}
	globalBary*=1./(number)rRefElem.num(0);

// 	compute global informations for scvf
	for(size_t i = 0; i < num_scvf(); ++i)
	{
		for (size_t j=0;j<m_vSCVF[i].numCo-1;j++){
			m_vSCVF[i].vGloPos[j]=vCornerCoords[rRefElem.id(dim-2,i,0,j)];
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
			m_vSCV[i].vGloPos[m_vSCV[i].numCorners-2-j]=vCornerCoords[rRefElem.id(dim-1,i,0,j)];
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
		     ElementNormal<face_type1,worldDim>(m_vSCV[i].Normal, m_vSCV[i].vGloPos);
		};
		// nodes are in reverse order therefore reverse sign to get outward normal
		m_vSCV[i].Normal*=-1;
	}
	
	}
	UG_CATCH_THROW("DimCRFVGeometry: geometric update failed.");
}


template <int TDim, int TWorldDim>
void DimCRFVGeometry<TDim, TWorldDim>::
update_boundary_faces(GridObject* pElem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
//	get grid
	Grid& grid = *(ish->grid());

//	vector of subset indices of side
	std::vector<int> vSubsetIndex;

//	get subset indices for sides (i.e. edge in 2d, faces in 3d)
	if(dim == 1) {
		std::vector<Vertex*> vVertex;
		CollectVertices(vVertex, grid, pElem);
		vSubsetIndex.resize(vVertex.size());
		for(size_t i = 0; i < vVertex.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vVertex[i]);
	}
	if(dim == 2) {
		std::vector<Edge*> vEdges;
		CollectEdgesSorted(vEdges, grid, pElem);
		vSubsetIndex.resize(vEdges.size());
		for(size_t i = 0; i < vEdges.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vEdges[i]);
	}
	if(dim == 3) {
		std::vector<Face*> vFaces;
		CollectFacesSorted(vFaces, grid, pElem);
		vSubsetIndex.resize(vFaces.size());
		for(size_t i = 0; i < vFaces.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vFaces[i]);
	}

	try{
//	const DimReferenceElement<dim>& rRefElem
	//	= ReferenceElementProvider::get<dim>(m_roid);

	DimReferenceMapping<dim, worldDim>& rMapping = ReferenceMappingProvider::get<dim, worldDim>(m_roid);
	rMapping.update(vCornerCoords);

	const LocalShapeFunctionSet<dim>& TrialSpace =
		LocalFiniteElementProvider::get<dim>(m_roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

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

	//	current number of bf
		size_t curr_bf = 0;

	//	loop sides of element
		for(size_t side = 0; side < vSubsetIndex.size(); ++side)
		{
		//	skip non boundary sides
			if(vSubsetIndex[side] != bndIndex) continue;

		//	number of corners of side
			// const int coOfSide = rRefElem.num(dim-1, side, 0); todo use somewhere?

		//	resize vector
			vBF.resize(curr_bf + 1);

		//  fill BF with data from associated SCV
			BF& bf = vBF[curr_bf];
			bf.nodeID =	m_vSCV[side].nodeID;

			bf.localIP =  m_vSCV[side].vLocIP;
			bf.globalIP = m_vSCV[side].vGlobIP;

			bf.Normal = m_vSCV[side].Normal;
			//	compute volume
			bf.Vol = VecTwoNorm(bf.Normal);

			bf.numCo = m_vSCV[side].numCorners-1;

			//	compute shapes and grads
			bf.numSH = TrialSpace.num_sh();
			TrialSpace.shapes(&(bf.vShape[0]), bf.localIP);
			TrialSpace.grads(&(bf.vLocalGrad[0]), bf.localIP);

			//	get reference mapping
			rMapping.jacobian_transposed_inverse(bf.JtInv, bf.localIP);
			bf.detj = rMapping.sqrt_gram_det(bf.localIP);

			//	compute global gradients
			for(size_t sh = 0 ; sh < bf.num_sh(); ++sh)
				MatVecMult(bf.vGlobalGrad[sh], bf.JtInv, bf.vLocalGrad[sh]);

			//	increase curr_bf
			++curr_bf;
		}
	}

	}
	UG_CATCH_THROW("DimCRFVGeometry: update failed.");
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//// Methods for CRFVGeometry class
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <	typename TElem, int TWorldDim>
CRFVGeometry<TElem, TWorldDim>::
CRFVGeometry()
	: m_pElem(NULL), m_rRefElem(Provider<ref_elem_type>::get()),
	  m_rTrialSpace(Provider<local_shape_fct_set_type>::get())
{
	update_local_data();
}

template <	typename TElem, int TWorldDim>
void CRFVGeometry<TElem, TWorldDim>::
update_local_data()
{
//  compute barycenter coordinates
	localBary = m_rRefElem.corner(0);
	for (size_t j=1;j<m_rRefElem.num(0);j++){
	   	localBary+=m_rRefElem.corner(j);
	}
	localBary*=1./(number)m_rRefElem.num(0);

// 	set up local informations for SubControlVolumeFaces (scvf)
// 	each scvf is associated to one vertex (2d) / edge (3d) of the element
	for(size_t i = 0; i < numSCVF; ++i)
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
// 	each scv is associated to one edge (2d) / face (3d) of the element
	for(size_t i = 0; i < numSCV; ++i)
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
	for(size_t i = 0; i < numSCVF; ++i)
	{
		m_rTrialSpace.shapes(&(m_vSCVF[i].vShape[0]), m_vSCVF[i].local_ip());
		m_rTrialSpace.grads(&(m_vSCVF[i].vLocalGrad[0]), m_vSCVF[i].local_ip());
	}

	for(size_t i = 0; i < numSCV; ++i)
	{
		m_rTrialSpace.shapes(&(m_vSCV[i].vShape[0]), m_vSCV[i].local_ip());
		m_rTrialSpace.grads(&(m_vSCV[i].vLocalGrad[0]), m_vSCV[i].local_ip());
	}

// 	copy ip positions in a list for Sub Control Volumes Faces (SCVF)
	for(size_t i = 0; i < num_scvf(); ++i)
		m_vLocSCVF_IP[i] = scvf(i).local_ip();
}


/// update data for given element
template <	typename TElem, int TWorldDim>
void CRFVGeometry<TElem, TWorldDim>::
update(GridObject* elem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
	UG_ASSERT(dynamic_cast<TElem*>(elem) != NULL, "Wrong element type.");
	TElem* pElem = static_cast<TElem*>(elem);

// 	If already update for this element, do nothing
	if(m_pElem == pElem) return; else m_pElem = pElem;

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

//	if no boundary subsets required, return
	if(num_boundary_subsets() == 0 || ish == NULL) return;
	else update_boundary_faces(pElem, vCornerCoords, ish);
}

template <	typename TElem, int TWorldDim>
void CRFVGeometry<TElem, TWorldDim>::
update_boundary_faces(GridObject* pElem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
//	get grid
	Grid& grid = *(ish->grid());

//	vector of subset indices of side
	std::vector<int> vSubsetIndex;

//	get subset indices for sides (i.e. edge in 2d, faces in 3d)
	if(dim == 1) {
		std::vector<Vertex*> vVertex;
		CollectVertices(vVertex, grid, pElem);
		vSubsetIndex.resize(vVertex.size());
		for(size_t i = 0; i < vVertex.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vVertex[i]);
	}
	if(dim == 2) {
		std::vector<Edge*> vEdges;
		CollectEdgesSorted(vEdges, grid, pElem);
		vSubsetIndex.resize(vEdges.size());
		for(size_t i = 0; i < vEdges.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vEdges[i]);
	}
	if(dim == 3) {
		std::vector<Face*> vFaces;
		CollectFacesSorted(vFaces, grid, pElem);
		vSubsetIndex.resize(vFaces.size());
		for(size_t i = 0; i < vFaces.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vFaces[i]);
	}

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

	//	current number of bf
		size_t curr_bf = 0;

	//	loop sides of element
		for(size_t side = 0; side < vSubsetIndex.size(); ++side)
		{
		//	skip non boundary sides
			if(vSubsetIndex[side] != bndIndex) continue;

		//	resize vector
			vBF.resize(curr_bf + 1);

		//  fill BF with data from associated SCV
			BF& bf = vBF[curr_bf];
			bf.nodeID =	m_vSCV[side].nodeID;

			bf.localIP =  m_vSCV[side].vLocIP;
			bf.globalIP = m_vSCV[side].vGlobIP;

			bf.Normal = m_vSCV[side].Normal;
			//	compute volume
			bf.Vol = VecTwoNorm(bf.Normal);

			bf.numCo = m_vSCV[side].numCorners-1;

			//	compute shapes and grads
			m_rTrialSpace.shapes(&(bf.vShape[0]), bf.localIP);
			m_rTrialSpace.grads(&(bf.vLocalGrad[0]), bf.localIP);

			//	get reference mapping
			m_mapping.jacobian_transposed_inverse(bf.JtInv, bf.localIP);
			bf.detj = m_mapping.sqrt_gram_det(bf.localIP);

			//	compute global gradients
			for(size_t sh = 0 ; sh < bf.num_sh(); ++sh)
				MatVecMult(bf.vGlobalGrad[sh], bf.JtInv, bf.vLocalGrad[sh]);

			//	increase curr_bf
			++curr_bf;
		}
	}
}

#ifdef UG_DIM_2
template class DimCRFVGeometry<2, 2>;
template class CRFVGeometry<Triangle, 2>;
template class CRFVGeometry<Quadrilateral, 2>;
#endif

#ifdef UG_DIM_3
template class DimCRFVGeometry<3, 3>;
template class CRFVGeometry<Triangle, 3>;
template class CRFVGeometry<Quadrilateral, 3>;

template class CRFVGeometry<Tetrahedron, 3>;
template class CRFVGeometry<Prism, 3>;
template class CRFVGeometry<Pyramid, 3>;
template class CRFVGeometry<Hexahedron, 3>;
#endif

} // end namespace ug
