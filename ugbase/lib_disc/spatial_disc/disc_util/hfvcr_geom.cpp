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
	
//  get grid
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
	} else {
		// dim = 3
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
				// first new scv replaces scv nr face
				size_t ind = face;
				// others are inserted at the end
				if (i>0) ind = numSCV+i-1;
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
				if (numFaceCo==3) m_vSCV[ind].Vol = ElementSize<scv_type0,worldDim>(m_vSCV[ind].vGloPos);
				else m_vSCV[ind].Vol = ElementSize<scv_type1,worldDim>(m_vSCV[ind].vGloPos);
				if (m_vSCV[ind].Vol<0) m_vSCV[ind].Vol *= -1;
				m_vSCV[ind].Normal = faceNormal;
				m_vSCV[ind].Normal *= (number) m_vSCV[ind].Vol / faceVol;
				m_rTrialSpace.shapes(&(m_vSCV[ind].vShape[0]), m_vSCV[ind].local_ip());
				m_rTrialSpace.grads(&(m_vSCV[ind].vLocalGrad[0]), m_vSCV[ind].local_ip());
				m_vSCV[ind].nodeID = numDofs+i;
			}
			// compute inner scv in triangular case
			if (numFaceCo==3){
				number volSum=m_vSCV[face].Vol;
				for (size_t j=0;j<2;j++){
					volSum+=m_vSCV[numSCV+j].Vol;
				}
				size_t ind = numSCV+2;
				m_vSCV[ind].Vol = faceVol - volSum;
				m_vSCV[ind].nodeID=numDofs+3;
				m_vSCV[ind].Normal = faceNormal;
				m_vSCV[ind].Normal *= (number) m_vSCV[ind].Vol / faceVol;
				m_vSCV[ind].vGlobIP = m_vSCV[face].vGloPos[1];
				m_vSCV[ind].vLocIP = m_vSCV[face].vLocPos[1];
				for (size_t j=0;j<2;j++){
					m_vSCV[ind].vGlobIP += m_vSCV[numSCV+j].vGloPos[1];
					m_vSCV[ind].vLocIP += m_vSCV[numSCV+j].vLocPos[1];
				}
				m_vSCV[ind].vGlobIP *= (number)1.0/3.0;
				m_vSCV[ind].vLocIP *= (number)1.0/3.0;
				m_vGlobUnkCoords[ind] = m_vSCV[ind].vGlobIP;
				m_vLocUnkCoords[ind] = m_vSCV[ind].vLocIP;
				m_vSCV[ind].numCorners = numFaceCo + 1;
				m_rTrialSpace.shapes(&(m_vSCV[ind].vShape[0]), m_vSCV[ind].local_ip());
				m_rTrialSpace.grads(&(m_vSCV[ind].vLocalGrad[0]), m_vSCV[ind].local_ip());
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
							scvID[j] = numDofs+k;
							break;
						}
					}
				}
				// look if edge has already been handled
				if (numConstrainedDofs>0){
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
				hEdge.scvfIndex=numSCVF;
				MathVector<worldDim> edgeMidGlob;
				MathVector<dim> edgeMidLoc;
				for (int d=0;d<worldDim;d++){
					edgeMidGlob[d] = 0.5 * (vCornerCoords[edgeCo[0]][d]  + vCornerCoords[edgeCo[1]][d]);
					edgeMidLoc[d] = 0.5 * (m_rRefElem.corner(edgeCo[0])[d] +  m_rRefElem.corner(edgeCo[1])[d]);
				}
				for (size_t j=0;j<2;j++){
					hEdge.associatedSCV[j] = scvID[j];
					if (from==face){
						m_vSCVF[numSCVF].From = scvID[j];
						m_vSCVF[numSCVF].To   = to;
						hEdge.from=true;
					} else {
						m_vSCVF[numSCVF].From = from;
						m_vSCVF[numSCVF].To 	= scvID[j];
						hEdge.from=false;
					}
					m_vSCVF[numSCVF].Normal = normal;
					m_vSCVF[numSCVF].vGloPos[0] = vCornerCoords[edgeCo[j]];
					m_vSCVF[numSCVF].vLocPos[0] = m_rRefElem.corner(edgeCo[j]);
					m_vSCVF[numSCVF].vGloPos[1] = edgeMidGlob;
					m_vSCVF[numSCVF].vLocPos[1] = edgeMidLoc;
					m_vSCVF[numSCVF].vGloPos[2] = globalBary;
					m_vSCVF[numSCVF].vLocPos[2] = localBary;
					AveragePositions(m_vSCVF[numSCVF].localIP, m_vSCVF[numSCVF].vLocPos, m_vSCVF[numSCVF].numCo);
					AveragePositions(m_vSCVF[numSCVF].globalIP, m_vSCVF[numSCVF].vGloPos, m_vSCVF[numSCVF].numCo);
					m_rTrialSpace.shapes(&(m_vSCVF[numSCVF].vShape[0]), m_vSCVF[numSCVF].local_ip());
					m_rTrialSpace.grads(&(m_vSCVF[numSCVF].vLocalGrad[0]), m_vSCVF[numSCVF].local_ip());
					numSCVF++;
				}
				handledEdges.push_back(hEdge);
			}
			// scvfs inside the face
			// insert remaining inner scvfs into positions of edge associated scvs
			for (size_t j=0;j<numFaceCo;j++){
				// replaces edge associated scvf
				size_t ii = faceEdge[j];
				if (ii==deleted){ 
					ii = numSCVF;
					numSCVF++;
				}
				if (numFaceCo==3){
					// compute inner scvfs in triangular case
					m_vSCVF[ii].From = numDofs+3;
					m_vSCVF[ii].To = numDofs+j;
					size_t ind = face;
					if (j>0) ind = numSCV+j-1;
					m_vSCVF[ii].vLocPos[0] = m_vSCV[ind].vLocPos[1];
					m_vSCVF[ii].vLocPos[1] = m_vSCV[ind].vLocPos[2];
					m_vSCVF[ii].vGloPos[0] = m_vSCV[ind].vGloPos[1];
					m_vSCVF[ii].vGloPos[1] = m_vSCV[ind].vGloPos[2];
				}else{
					// compute inner scvfs in quadrilateral case
					m_vSCVF[ii].To = numDofs+j;
					m_vSCVF[ii].From = numDofs + ((j+1) % 4);
					for (int d=0;d<worldDim;d++){
						m_vSCVF[ii].vLocPos[0][d] = 0.5*( m_rRefElem.corner(faceCo[j])[d] + m_rRefElem.corner(faceCo[(j+1) % 4])[d] );
						m_vSCVF[ii].vGloPos[0][d] = 0.5*( vCornerCoords[faceCo[j]][d] + vCornerCoords[faceCo[(j+1) % 4]][d] );
					}
					m_vSCVF[ii].vLocPos[1] = localMidpoint;
					m_vSCVF[ii].vGloPos[1] = globalMidpoint;
				}
				m_vSCVF[ii].vLocPos[2] = localBary;
				m_vSCVF[ii].vGloPos[2] = globalBary;	
				ElementNormal<face_type0,worldDim>(m_vSCVF[ii].Normal,m_vSCVF[ii].vGloPos);
				AveragePositions(m_vSCVF[ii].globalIP,m_vSCVF[ii].vGloPos,3);
				AveragePositions(m_vSCVF[ii].localIP,m_vSCVF[ii].vLocPos,3);
				m_rTrialSpace.shapes(&(m_vSCVF[ii].vShape[0]), m_vSCVF[ii].local_ip());
				m_rTrialSpace.grads(&(m_vSCVF[ii].vLocalGrad[0]), m_vSCVF[ii].local_ip());
			}
			// insert new constrained dof object
			m_vCD[numConstrainedDofs].i = face;
			m_vCD[numConstrainedDofs].numConstrainingDofs = 4;
			for (size_t i=0;i<4;i++){
				m_vCD[numConstrainedDofs].cDofInd[i] = numDofs+i;
				size_t ind=face;
				if (i>0) ind = numSCV+i-1;
				m_vCD[numConstrainedDofs].cDofWeights[i] = (number)m_vSCV[ind].Vol / faceVol;
			}
			numSCV+=3;
			numDofs+=4;
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
	UG_LOG("constrained dof indices:\n");
	for (size_t i=0;i<numConstrainedDofs;i++){
		UG_LOG("constrained index = " << m_vCD[i].i << "\n");
		UG_LOG("constraining indices: ");
		for (size_t j=0;j< m_vCD[i].numConstrainingDofs;j++) UG_LOG(m_vCD[i].cDofInd[j] << " ");
		UG_LOG("\n");
		UG_LOG("weights: ");
		for (size_t j=0;j< m_vCD[i].numConstrainingDofs;j++) UG_LOG(m_vCD[i].cDofWeights[j] << " ");
		UG_LOG("\n");
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
