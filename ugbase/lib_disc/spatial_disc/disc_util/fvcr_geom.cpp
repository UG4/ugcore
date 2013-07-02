/*
 * fvcr_geom.cpp
 *
 *  Created on: 21.06.2012
 *      Author: Christian Wehner
 *
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
	const DimReferenceElement<dim>& rRefElem
		= ReferenceElementProvider::get<dim>(m_roid);

//	set number of scvf / scv of this roid
	m_numSCV = rRefElem.num(dim-1); // number of faces
	m_numSCVF = rRefElem.num(1); // number of edges

//  compute barycenter coordinates
	localBary = rRefElem.corner(0);
	for (size_t j=1;j<rRefElem.num(0);j++){
	   	localBary+=rRefElem.corner(j);
	}
	localBary*=1./(number)rRefElem.num(0);

// 	set up local informations for SubControlVolumeFaces (scvf)
// 	each scvf is associated to one edge of the element
	for(size_t i = 0; i < m_numSCVF; ++i)
	{
	//	this scvf separates the given edges/faces
		m_vSCVF[i].From = rRefElem.id(dim-2, i, dim-1, 0);// todo handle dim==1
		m_vSCVF[i].To = rRefElem.id(dim-2, i, dim-1, 1);

		for (size_t j=0;j<m_vSCVF[i].numCo-1;j++){
			m_vSCVF[i].vLocPos[j]=rRefElem.corner(rRefElem.id(dim-2,i,0,j));
		}

		m_vSCVF[i].vLocPos[m_vSCVF[i].numCo-1]=localBary;

		AveragePositions(m_vSCVF[i].localIP, m_vSCVF[i].vLocPos, m_vSCVF[i].numCo);
	}

// 	set up local informations for SubControlVolumes (scv)
// 	each scv is associated to one corner of the element
	for(size_t i = 0; i < m_numSCV; ++i)
	{
	//	store associated node
		m_vSCV[i].nodeID = i;

		m_vSCV[i].numCorners = rRefElem.num(dim-1,i,0)+1;
		for (int j=0;j<m_vSCV[i].numCorners-1;j++){
			m_vSCV[i].vLocPos[m_vSCV[i].numCorners-2-j]=rRefElem.corner(rRefElem.id(dim-1,i,0,j));
		}
		AveragePositions(m_vLocUnkCoords[i], m_vSCV[i].vLocPos, m_vSCV[i].numCorners-1);
		m_vSCV[i].vLocIP = m_vLocUnkCoords[i];
		m_vSCV[i].vLocPos[m_vSCV[i].numCorners-1]=localBary;
	}

	/////////////////////////
	// Shapes and Derivatives
	/////////////////////////

	const LocalShapeFunctionSet<dim>& TrialSpace =
		LocalFiniteElementProvider::get<dim>(m_roid, LFEID(LFEID::CROUZEIX_RAVIART, dim, 1));

	m_nsh = TrialSpace.num_sh();

	for(size_t i = 0; i < m_numSCVF; ++i)
	{
		m_vSCVF[i].numSH = TrialSpace.num_sh();
		TrialSpace.shapes(&(m_vSCVF[i].vShape[0]), m_vSCVF[i].local_ip());
		TrialSpace.grads(&(m_vSCVF[i].vLocalGrad[0]), m_vSCVF[i].local_ip());
	}

	for(size_t i = 0; i < m_numSCV; ++i)
	{
		m_vSCV[i].numSH = TrialSpace.num_sh();
		TrialSpace.shapes(&(m_vSCV[i].vShape[0]), m_vSCV[i].local_ip());
		TrialSpace.grads(&(m_vSCV[i].vLocalGrad[0]), m_vSCV[i].local_ip());
	}

	}
	UG_CATCH_THROW("DimCRFVGeometry: update failed.");

// 	copy ip positions in a list for Sub Control Volumes Faces (SCVF)
	for(size_t i = 0; i < m_numSCVF; ++i)
		m_vLocSCVF_IP[i] = scvf(i).local_ip();
}


/// update data for given element
template <int TDim, int TWorldDim>
void DimCRFVGeometry<TDim, TWorldDim>::
update(GeometricObject* pElem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
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

template <int TDim, int TWorldDim>
void DimCRFVGeometry<TDim, TWorldDim>::
update_boundary_faces(GeometricObject* pElem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
//	get grid
	Grid& grid = *(ish->grid());

//	vector of subset indices of side
	std::vector<int> vSubsetIndex;

//	get subset indices for sides (i.e. edge in 2d, faces in 3d)
	if(dim == 1) {
		std::vector<VertexBase*> vVertex;
		CollectVertices(vVertex, grid, pElem);
		vSubsetIndex.resize(vVertex.size());
		for(size_t i = 0; i < vVertex.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vVertex[i]);
	}
	if(dim == 2) {
		std::vector<EdgeBase*> vEdges;
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
// 	each scvf is associated to one edge of the element
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
// 	each scv is associated to one corner of the element
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
update(GeometricObject* elem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
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
update_boundary_faces(GeometricObject* pElem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
//	get grid
	Grid& grid = *(ish->grid());

//	vector of subset indices of side
	std::vector<int> vSubsetIndex;

//	get subset indices for sides (i.e. edge in 2d, faces in 3d)
	if(dim == 1) {
		std::vector<VertexBase*> vVertex;
		CollectVertices(vVertex, grid, pElem);
		vSubsetIndex.resize(vVertex.size());
		for(size_t i = 0; i < vVertex.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vVertex[i]);
	}
	if(dim == 2) {
		std::vector<EdgeBase*> vEdges;
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
