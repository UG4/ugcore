/*
 * hfv1_geom.cpp
 *
 *  Created on: 31.05.2011
 *      Author: andreasvogel
 */

#include "hfv1_geom.h"
#include "common/util/provider.h"

namespace ug{

template <typename TElem, int TWorldDim>
HFV1Geometry<TElem, TWorldDim>::
HFV1Geometry()
	: m_pElem(NULL),  m_rRefElem(Provider<ref_elem_type>::get())
{
	// local corners
	m_vSCV.resize(m_numNaturalSCV);
	m_locMid[0].resize(m_numNaturalSCV);
	for(size_t i = 0; i < m_numNaturalSCV; ++i)
	{
		m_vSCV[i].nodeId = i;
		m_vSCV[i].m_vLocPos[0] = m_rRefElem.corner(i);
		m_locMid[0][i] = m_rRefElem.corner(i);
	}

	// compute center
	m_locMid[dim].resize(1);
	m_gloMid[dim].resize(1);
	m_locMid[dim][0] = 0.0;
	for(size_t i = 0; i < m_locMid[0].size(); ++i)
	{
		m_locMid[dim][0] += m_locMid[0][i];
	}
	m_locMid[dim][0] *= 1./(m_locMid[0].size());
}


template <typename TElem, int TWorldDim>
void HFV1Geometry<TElem, TWorldDim>::
update(GridObject* elem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
	UG_ASSERT(dynamic_cast<TElem*>(elem) != NULL, "Wrong element type.");
	TElem* pElem = static_cast<TElem*>(elem);

	// If already update for this element, do nothing
	if(m_pElem == pElem) return;
	else m_pElem = pElem;

	// get grid
	Grid& grid = *(ish->grid());

	// reset to natural nodes
	m_gloMid[0].resize(m_numNaturalSCV);
	m_locMid[0].resize(m_numNaturalSCV);

	// remember global position of nodes
	for(size_t i = 0; i < m_numNaturalSCV; ++i)
		m_gloMid[0][i] = vCornerCoords[i];

	// compute center
	m_gloMid[dim][0] = 0.0;
	for(size_t i = 0; i < m_gloMid[0].size(); ++i)
	{
		m_gloMid[dim][0] += m_gloMid[0][i];
	}
	m_gloMid[dim][0] *= 1./(m_gloMid[0].size());

	// get natural edges (and faces if in 3d)
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, grid, pElem);

	// compute Nodes
	m_vSCVF.clear();
	m_vNewEdgeInfo.clear();
	m_vNatEdgeInfo.clear(); m_vNatEdgeInfo.resize(m_numNaturalSCVF);
    UG_ASSERT(vEdges.size() == m_numNaturalSCVF, "Not correct number of edges found, only " << vEdges.size() << "Edges");
	for(size_t i = 0; i < vEdges.size(); ++i)
	{
		// natural ids of end of edge
		const size_t from = m_rRefElem.id(1, i, 0, 0);
		const size_t to = m_rRefElem.id(1, i, 0, 1);

		// choose whether to insert two or one new edge
		switch(vEdges[i]->container_section())
		{
		case CSEDGE_CONSTRAINED_EDGE:
		case CSEDGE_REGULAR_EDGE:
			// classic case: Just set corner ids
			if(dim == 2)
			{
				const size_t numSCVF = m_vSCVF.size();
				m_vSCVF.resize(numSCVF + 1);
				m_vSCVF[numSCVF].m_from = from;
				m_vSCVF[numSCVF].m_to = to;
			}
			if(dim == 3)
			{
				const size_t numNewEdgeInfo = m_vNewEdgeInfo.size();
				m_vNatEdgeInfo[i].numChildEdges = 1;
				m_vNatEdgeInfo[i].childEdge[0] = numNewEdgeInfo;

				m_vNewEdgeInfo.resize(numNewEdgeInfo + 1);
				m_vNewEdgeInfo[numNewEdgeInfo].m_from = from;
				m_vNewEdgeInfo[numNewEdgeInfo].m_to = to;
			}
			break;

		case CSEDGE_CONSTRAINING_EDGE:
			{
				// insert hanging node in list of nodes
				const size_t newNodeId = m_gloMid[0].size();
				m_gloMid[0].resize(newNodeId + 1);
				m_locMid[0].resize(newNodeId + 1);
				VecInterpolateLinear(	m_gloMid[0].back(),
										m_gloMid[0][from],
										m_gloMid[0][to],
										0.5);
				VecInterpolateLinear(	m_locMid[0].back(),
										m_locMid[0][from],
										m_locMid[0][to],
										0.5);

				if(dim == 2)
				{
					// insert two edges with nodeIds
					const size_t numSCVF = m_vSCVF.size();
					m_vSCVF.resize(numSCVF + 2);

					m_vSCVF[numSCVF].m_from = from;
					m_vSCVF[numSCVF].m_to = newNodeId;

					m_vSCVF[numSCVF+1].m_from = newNodeId;
					m_vSCVF[numSCVF+1].m_to = to;
				}
				if(dim == 3)
				{
					// Mapping NaturalEdges -> New Edges
					const size_t numNewEdgeInfo = m_vNewEdgeInfo.size();
					m_vNatEdgeInfo[i].nodeId = newNodeId;
					m_vNatEdgeInfo[i].numChildEdges = 2;
					m_vNatEdgeInfo[i].childEdge[0] = numNewEdgeInfo;
					m_vNatEdgeInfo[i].childEdge[1] = numNewEdgeInfo + 1;

					m_vNewEdgeInfo.resize(numNewEdgeInfo + 2);

					m_vNewEdgeInfo[numNewEdgeInfo].m_from = from;
					m_vNewEdgeInfo[numNewEdgeInfo].m_to = newNodeId;

					m_vNewEdgeInfo[numNewEdgeInfo+1].m_from = newNodeId;
					m_vNewEdgeInfo[numNewEdgeInfo+1].m_to = to;
				}
			}
			break;

		default: UG_THROW("Cannot detect type of edge.");
		}
	}

	// for 3d case also check faces for hanging nodes
	if(dim == 3)
	{
		std::vector<Face*> vFaces;
		CollectFacesSorted(vFaces, grid, pElem);

		// compute Nodes
		MathVector<dim> locSideMid;
		MathVector<worldDim> gloSideMid;
		for(size_t i = 0; i < vFaces.size(); ++i)
		{
			///////////
			// case QUADRILATERAL with all edges hanging and hanging node in middle
			///////////
			if(vFaces[i]->container_section() == CSFACE_CONSTRAINING_QUADRILATERAL)
			{
				// insert hanging node in list of nodes
				const size_t newNodeId = m_gloMid[0].size();
				m_gloMid[0].resize(newNodeId + 1);
				m_locMid[0].resize(newNodeId + 1);

				// compute position of new (hanging) node
				compute_side_midpoints(i, m_locMid[0][newNodeId], m_gloMid[0][newNodeId]);

				// loop constrained faces
				for(size_t j = 0; j < m_rRefElem.num(2, i, 1); ++j)
				{
					const size_t jplus1 = (j+1)%4;

					// natural edges
					const size_t natEdId1 = m_rRefElem.id(2, i, 1, j);
					const size_t natEdId2 = m_rRefElem.id(2, i, 1, jplus1);

					// corner of the face
					const size_t cornerId = m_rRefElem.id(2,i, 0, jplus1);

					// refined edges that belong to this face
					const size_t edId1 = get_child_edge_of_corner(natEdId1, cornerId);
					const size_t edId2 = get_child_edge_of_corner(natEdId2, cornerId);

					// nodes of hanging edges
					const size_t hangEdNodeId1 = m_vNatEdgeInfo[natEdId1].node_id();
					const size_t hangEdNodeId2 = m_vNatEdgeInfo[natEdId2].node_id();

					// mid point of hanging side
					compute_side_midpoints(	cornerId, newNodeId,
											hangEdNodeId1, hangEdNodeId2,
											locSideMid, gloSideMid);

					// add side midpoint to already existing scvf of this side
					const size_t numSCVF = m_vSCVF.size();
					m_vSCVF.resize(numSCVF + 4);

					m_vSCVF[numSCVF].m_from = m_vNewEdgeInfo[edId1].from();
					m_vSCVF[numSCVF].m_to = m_vNewEdgeInfo[edId1].to();
					m_vSCVF[numSCVF].m_vLocPos[2] = locSideMid;
					m_vSCVF[numSCVF].m_vGloPos[2] = gloSideMid;

					m_vSCVF[numSCVF+1].m_from = m_vNewEdgeInfo[edId2].from();
					m_vSCVF[numSCVF+1].m_to = m_vNewEdgeInfo[edId2].to();
					m_vSCVF[numSCVF+1].m_vLocPos[2] = locSideMid;
					m_vSCVF[numSCVF+1].m_vGloPos[2] = gloSideMid;

					m_vSCVF[numSCVF+2].m_from = hangEdNodeId1;
					m_vSCVF[numSCVF+2].m_to = newNodeId;
					m_vSCVF[numSCVF+2].m_vLocPos[2] = locSideMid;
					m_vSCVF[numSCVF+2].m_vGloPos[2] = gloSideMid;

					m_vSCVF[numSCVF+3].m_from = hangEdNodeId2;
					m_vSCVF[numSCVF+3].m_to = newNodeId;
					m_vSCVF[numSCVF+3].m_vLocPos[2] = locSideMid;
					m_vSCVF[numSCVF+3].m_vGloPos[2] = gloSideMid;
				}
			}
			///////////
			// case TRIANGLE with all edges hanging, that matches a refined element on other side
			///////////
			else if (vFaces[i]->container_section() == CSFACE_CONSTRAINING_TRIANGLE)
			{
				bool bAllConstraining = true;
				for(size_t j = 0; j < m_rRefElem.num(2, i, 1); ++j)
					if(vEdges[m_rRefElem.id(2, i, 1, j)]->container_section() != CSEDGE_CONSTRAINING_EDGE)
						bAllConstraining = false;

				if(!bAllConstraining) continue;

				// compute position of new (hanging) node
				compute_side_midpoints(i, locSideMid, gloSideMid);

				// loop constrained faces
				for(size_t j = 0; j < m_rRefElem.num(2, i, 1); ++j)
				{
					const size_t jplus1 = (j+1)%3;

					// natural edges
					const size_t natEdId1 = m_rRefElem.id(2, i, 1, j);
					const size_t natEdId2 = m_rRefElem.id(2, i, 1, jplus1);

					// corner of the face
					const size_t cornerId = m_rRefElem.id(2,i, 0, jplus1);

					// nodes of hanging edges
					const size_t hangEdNodeId1 = m_vNatEdgeInfo[natEdId1].node_id();
					const size_t hangEdNodeId2 = m_vNatEdgeInfo[natEdId2].node_id();

					MathVector<dim> locSmallSideMid;
					MathVector<worldDim> gloSmallSideMid;

					// mid point of hanging side
					compute_side_midpoints(	cornerId, hangEdNodeId1, hangEdNodeId2,
											locSmallSideMid, gloSmallSideMid);

					// add side midpoint to already existing scvf of this side
					const size_t numSCVF = m_vSCVF.size();
					m_vSCVF.resize(numSCVF + 4);

					m_vSCVF[numSCVF].m_from = hangEdNodeId1;
					m_vSCVF[numSCVF].m_to = hangEdNodeId2;
					m_vSCVF[numSCVF].m_vLocPos[2] = locSideMid;
					m_vSCVF[numSCVF].m_vGloPos[2] = gloSideMid;

					m_vSCVF[numSCVF+1].m_from = hangEdNodeId1;
					m_vSCVF[numSCVF+1].m_to = hangEdNodeId2;
					m_vSCVF[numSCVF+1].m_vLocPos[2] = locSmallSideMid;
					m_vSCVF[numSCVF+1].m_vGloPos[2] = gloSmallSideMid;

					m_vSCVF[numSCVF+2].m_from = hangEdNodeId1;
					m_vSCVF[numSCVF+2].m_to = cornerId;
					m_vSCVF[numSCVF+2].m_vLocPos[2] = locSmallSideMid;
					m_vSCVF[numSCVF+2].m_vGloPos[2] = gloSmallSideMid;

					m_vSCVF[numSCVF+3].m_from = cornerId;
					m_vSCVF[numSCVF+3].m_to = hangEdNodeId2;
					m_vSCVF[numSCVF+3].m_vLocPos[2] = locSmallSideMid;
					m_vSCVF[numSCVF+3].m_vGloPos[2] = gloSmallSideMid;
				}
			}
			//////////
			// other cases: Not all edges hanging (i.e. neighbor not refined)
			///////////
			else
			{
				// compute side midpoint
				compute_side_midpoints(i, locSideMid, gloSideMid);

				// connect all edges with side midpoint
				for(size_t j = 0; j < m_rRefElem.num(2, i, 1); ++j)
				{
					const size_t natEdgeId = m_rRefElem.id(2, i, 1, j);
					for(size_t e = 0; e < m_vNatEdgeInfo[natEdgeId].num_child_edges(); ++e)
					{
						const size_t edgeId = m_vNatEdgeInfo[natEdgeId].child_edge(e);

						const size_t numSCVF = m_vSCVF.size();
						m_vSCVF.resize(numSCVF + 1);

						m_vSCVF[numSCVF].m_from = m_vNewEdgeInfo[edgeId].from();
						m_vSCVF[numSCVF].m_to = m_vNewEdgeInfo[edgeId].to();
						m_vSCVF[numSCVF].m_vGloPos[2] = gloSideMid;
						m_vSCVF[numSCVF].m_vLocPos[2] = locSideMid;
					}
				}
			}
		}
	}

	// resize number of scv (== number of nodes)
	if(dim <= 2)
	{
		m_vSCV.resize(m_gloMid[0].size());
	}
	else
	{
		m_vSCV.clear();
	}

	// compute scvf
	for(size_t i = 0; i < m_vSCVF.size(); ++i)
	{
		// compute midpoints of edges
		{
			const size_t from = m_vSCVF[i].m_from;
			const size_t to = m_vSCVF[i].m_to;

			VecInterpolateLinear(	m_vSCVF[i].m_vLocPos[0],
									m_locMid[0][from],
									m_locMid[0][to],
									0.5);
			VecInterpolateLinear(	m_vSCVF[i].m_vGloPos[0],
									m_gloMid[0][from],
									m_gloMid[0][to],
									0.5);
		}

		// set center of elem as part of scvf
		m_vSCVF[i].m_vGloPos[1] = m_gloMid[dim][0];
		m_vSCVF[i].m_vLocPos[1] = m_locMid[dim][0];

		// integration point
		AveragePositions(m_vSCVF[i].localIP, m_vSCVF[i].m_vLocPos, SCVF::m_numCorners);
		AveragePositions(m_vSCVF[i].globalIP, m_vSCVF[i].m_vGloPos, SCVF::m_numCorners);

		// normal
		HangingNormalOnSCVF<ref_elem_type, worldDim>(m_vSCVF[i].Normal, &(m_vSCVF[i].m_vGloPos[0]));

		// TODO: In 3D check orientation
		if(dim == 3)
		{
			const size_t from = m_vSCVF[i].m_from;
			const size_t to = m_vSCVF[i].m_to;

			MathVector<worldDim> diff;
			VecSubtract(diff, m_gloMid[0][to], m_gloMid[0][from]);

			if(VecDot(diff, m_vSCVF[i].Normal) < 0)
			{
				m_vSCVF[i].m_from = to;
				m_vSCVF[i].m_to = from;
			}
		}

		// write edge midpoints to as corners of scv
		if(dim == 2)
		{
			const size_t from = m_vSCVF[i].m_from;
			const size_t to = m_vSCVF[i].m_to;

			m_vSCV[from].m_vLocPos[1] = m_vSCVF[i].m_vLocPos[0];
			m_vSCV[from].m_vGloPos[1] = m_vSCVF[i].m_vGloPos[0];

			m_vSCV[to].m_vLocPos[3] = m_vSCVF[i].m_vLocPos[0];
			m_vSCV[to].m_vGloPos[3] = m_vSCVF[i].m_vGloPos[0];
		}
	}

	// compute size of scv
	if(dim <= 2)
	{
		for(size_t i = 0; i < m_vSCV.size(); ++i)
		{
			// set node id
			m_vSCV[i].nodeId = i;

			// start at node
			m_vSCV[i].m_vLocPos[0] = m_locMid[0][i];
			m_vSCV[i].m_vGloPos[0] = m_gloMid[0][i];

			if(dim == 1)
			{
				// center of element
				m_vSCV[i].m_vLocPos[1] = m_locMid[dim][0];
				m_vSCV[i].m_vGloPos[1] = m_gloMid[dim][0];
			}
			else if (dim == 2)
			{
				// center of element
				m_vSCV[i].m_vLocPos[2] = m_locMid[dim][0];
				m_vSCV[i].m_vGloPos[2] = m_gloMid[dim][0];
			}

           m_vSCV[i].vol = ElementSize<typename SCV::scv_type, worldDim>(&(m_vSCV[i].m_vGloPos[0]));;
		}
	}

	if(dim == 3)
	{
		m_vSCV.resize(m_vSCVF.size() * 2);

		for(size_t i = 0; i < m_vSCVF.size(); ++i)
        {
			for(size_t n = 0; n < 3; ++n)
			{
				m_vSCV[2*i + 0].m_vLocPos[n+1] = m_vSCVF[i].m_vLocPos[n];
				m_vSCV[2*i + 0].m_vGloPos[n+1] = m_vSCVF[i].m_vGloPos[n];
				m_vSCV[2*i + 1].m_vLocPos[n+1] = m_vSCVF[i].m_vLocPos[n];
				m_vSCV[2*i + 1].m_vGloPos[n+1] = m_vSCVF[i].m_vGloPos[n];
			}
			const size_t from = m_vSCVF[i].m_from;
			const size_t to = m_vSCVF[i].m_to;

			m_vSCV[2*i + 0].m_vLocPos[0] = m_locMid[0][from];
			m_vSCV[2*i + 0].m_vGloPos[0] = m_gloMid[0][from];
			m_vSCV[2*i + 0].nodeId = from;
			m_vSCV[2*i + 0].m_numCorners = 4;

			m_vSCV[2*i + 1].m_vLocPos[0] = m_locMid[0][to];
			m_vSCV[2*i + 1].m_vGloPos[0] = m_gloMid[0][to];
			m_vSCV[2*i + 1].nodeId = to;
			m_vSCV[2*i + 1].m_numCorners = 4;

			m_vSCV[2*i + 0].vol = ElementSize<typename SCV::scv_type, worldDim>(&(m_vSCV[2*i + 0].m_vGloPos[0]));;
			m_vSCV[2*i + 1].vol = ElementSize<typename SCV::scv_type, worldDim>(&(m_vSCV[2*i + 1].m_vGloPos[0]));;
		}
	}

	/////////////////////////
	// Shapes and Derivatives
	/////////////////////////
	m_rMapping.update(vCornerCoords);

	const size_t num_sh = ref_elem_type::numCorners;
	m_numSh = num_sh;

	for(size_t i = 0; i < num_scvf(); ++i)
	{
		m_rMapping.jacobian_transposed_inverse(m_vSCVF[i].JtInv, m_vSCVF[i].localIP);
		m_vSCVF[i].detj = m_rMapping.sqrt_gram_det(m_vSCVF[i].localIP);

		const LocalShapeFunctionSet<ref_elem_type::dim>& TrialSpace =
				LocalFiniteElementProvider::
					get<ref_elem_type::dim>
						(ref_elem_type::REFERENCE_OBJECT_ID,
						 LFEID(LFEID::LAGRANGE, ref_elem_type::dim, 1));

		m_vSCVF[i].vShape.resize(num_sh);
		m_vSCVF[i].localGrad.resize(num_sh);
		m_vSCVF[i].globalGrad.resize(num_sh);

		TrialSpace.shapes(&(m_vSCVF[i].vShape[0]), m_vSCVF[i].localIP);
		TrialSpace.grads(&(m_vSCVF[i].localGrad[0]), m_vSCVF[i].localIP);

		for(size_t sh = 0 ; sh < num_sh; ++sh)
		{
			MatVecMult((m_vSCVF[i].globalGrad)[sh], m_vSCVF[i].JtInv, (m_vSCVF[i].localGrad)[sh]);
		}

	}


	///////////////////////////
	// Copy ip pos in list
	///////////////////////////

// 	loop Sub Control Volumes (SCV)
	m_vGlobSCVIP.clear();
	m_vLocSCVIP.clear();
	for(size_t i = 0; i < num_scv(); ++i)
	{
	//	get current SCV
		const SCV& rSCV = scv(i);

	//	copy
		m_vGlobSCVIP.push_back(rSCV.global_ip());
		m_vLocSCVIP.push_back(rSCV.local_ip());
	}

// 	loop Sub Control Volumes Faces (SCVF)
	m_vGlobSCVFIP.clear();
	m_vLocSCVFIP.clear();
	for(size_t i = 0; i < num_scvf(); ++i)
	{
	//	get current SCVF
		const SCVF& rSCVF = scvf(i);

	//  copy
		m_vGlobSCVFIP.push_back(rSCVF.global_ip());
		m_vLocSCVFIP.push_back(rSCVF.local_ip());
	}

//	print();
}

// debug output
template <typename TElem, int TWorldDim>
void HFV1Geometry<TElem, TWorldDim>::print()
{
	UG_LOG("\nFVG hanging debug output\n");
	UG_LOG("LocalCenter=" << m_locMid << ", globalCenter="<<m_gloMid<<"\n");
	for(size_t i = 0; i < m_vSCV.size(); ++i)
	{
		UG_LOG(i<<" SCV: ");
		UG_LOG("node_id=" << m_vSCV[i].node_id());
		UG_LOG(", local_pos="<< m_vSCV[i].local_ip());
		UG_LOG(", global_pos="<< m_vSCV[i].global_ip());
		UG_LOG(", vol=" << m_vSCV[i].volume());
/*		UG_LOG("\n    localCorner=" << m_vSCV[i].m_vLocPos[0]);
		UG_LOG(", localSide1=" << m_vSCV[i].m_vLocPos[1]);
		UG_LOG(", localCenter=" << m_vSCV[i].m_vLocPos[2]);
		UG_LOG(", localSide2=" << m_vSCV[i].m_vLocPos[3]);
		UG_LOG("\n    globalCorner=" << m_vSCV[i].m_vGloPos[0]);
		UG_LOG(", globalSide1=" << m_vSCV[i].m_vGloPos[1]);
		UG_LOG(", globalCenter=" << m_vSCV[i].m_vGloPos[2]);
		UG_LOG(", globalSide2=" << m_vSCV[i].m_vGloPos[3]);
*/
		UG_LOG("\n");
	}
	UG_LOG("\n");
	for(size_t i = 0; i < m_vSCVF.size(); ++i)
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

template <int TDim, int TWorldDim>
void DimHFV1Geometry<TDim, TWorldDim>::
update_local_data()
{
	//	remember new roid
	m_roid = (ReferenceObjectID) m_pElem->reference_object_id();
	
	// get new reference element
	m_rRefElem = ReferenceElementProvider::get<dim>(m_roid);
	
	// get new reference mapping
	m_rMapping = &ReferenceMappingProvider::get<dim, worldDim>(m_roid);

	m_numNaturalSCV = (m_roid != ROID_PYRAMID) ? m_rRefElem.num(0) : 8; // number of corners
	m_numNaturalSCVF = (m_roid != ROID_PYRAMID) ? m_rRefElem.num(1) : 12; // number of edges
	
	// local corners
	m_vSCV.resize(m_numNaturalSCV);
	m_locMid[0].resize(m_numNaturalSCV);
	for(size_t i = 0; i < m_numNaturalSCV; ++i)
	{
		m_vSCV[i].nodeId = i;
		m_vSCV[i].m_vLocPos[0] = m_rRefElem.corner(i);
		m_locMid[0][i] = m_rRefElem.corner(i);
	}

	// compute center
	m_locMid[dim].resize(1);
	m_gloMid[dim].resize(1);
	m_locMid[dim][0] = 0.0;
	for(size_t i = 0; i < m_locMid[0].size(); ++i)
	{
		m_locMid[dim][0] += m_locMid[0][i];
	}
	m_locMid[dim][0] *= 1./(m_locMid[0].size());
}


template <int TDim, int TWorldDim>
void DimHFV1Geometry<TDim, TWorldDim>::
update(GridObject* pElem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
	// 	If already update for this element, do nothing
	if(m_pElem == pElem) return; else m_pElem = pElem;
	
	//	refresh local data, if different roid given
	if (m_roid != m_pElem->reference_object_id()) update_local_data();
	
	// get grid
	Grid& grid = *(ish->grid());

	// reset to natural nodes
	m_gloMid[0].resize(m_numNaturalSCV);
	m_locMid[0].resize(m_numNaturalSCV);

	// remember global position of nodes
	for(size_t i = 0; i < m_numNaturalSCV; ++i)
		m_gloMid[0][i] = vCornerCoords[i];

	// compute center
	m_gloMid[dim][0] = 0.0;
	for(size_t i = 0; i < m_gloMid[0].size(); ++i)
	{
		m_gloMid[dim][0] += m_gloMid[0][i];
	}
	m_gloMid[dim][0] *= 1./(m_gloMid[0].size());

	// get natural edges (and faces if in 3d)
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, grid, pElem);

	// compute Nodes
	m_vSCVF.clear();
	m_vNewEdgeInfo.clear();
	m_vNatEdgeInfo.clear(); m_vNatEdgeInfo.resize(m_numNaturalSCVF);
    UG_ASSERT(vEdges.size() == m_numNaturalSCVF, "Not correct number of edges found, only " << vEdges.size() << "Edges");
	for(size_t i = 0; i < vEdges.size(); ++i)
	{
		// natural ids of end of edge
		const size_t from = m_rRefElem.id(1, i, 0, 0);
		const size_t to = m_rRefElem.id(1, i, 0, 1);

		// choose weather to insert two or one new edge
		switch(vEdges[i]->container_section())
		{
		case CSEDGE_CONSTRAINED_EDGE:
		case CSEDGE_REGULAR_EDGE:
			// classic case: Just set corner ids
			if(dim == 2)
			{
				const size_t numSCVF = m_vSCVF.size();
				m_vSCVF.resize(numSCVF + 1);
				m_vSCVF[numSCVF].m_from = from;
				m_vSCVF[numSCVF].m_to = to;
			}
			if(dim == 3)
			{
				const size_t numNewEdgeInfo = m_vNewEdgeInfo.size();
				m_vNatEdgeInfo[i].numChildEdges = 1;
				m_vNatEdgeInfo[i].childEdge[0] = numNewEdgeInfo;

				m_vNewEdgeInfo.resize(numNewEdgeInfo + 1);
				m_vNewEdgeInfo[numNewEdgeInfo].m_from = from;
				m_vNewEdgeInfo[numNewEdgeInfo].m_to = to;
			}
			break;

		case CSEDGE_CONSTRAINING_EDGE:
			{
				// insert hanging node in list of nodes
				const size_t newNodeId = m_gloMid[0].size();
				m_gloMid[0].resize(newNodeId + 1);
				m_locMid[0].resize(newNodeId + 1);
				VecInterpolateLinear(	m_gloMid[0].back(),
										m_gloMid[0][from],
										m_gloMid[0][to],
										0.5);
				VecInterpolateLinear(	m_locMid[0].back(),
										m_locMid[0][from],
										m_locMid[0][to],
										0.5);

				if(dim == 2)
				{
					// insert two edges with nodeIds
					const size_t numSCVF = m_vSCVF.size();
					m_vSCVF.resize(numSCVF + 2);

					m_vSCVF[numSCVF].m_from = from;
					m_vSCVF[numSCVF].m_to = newNodeId;

					m_vSCVF[numSCVF+1].m_from = newNodeId;
					m_vSCVF[numSCVF+1].m_to = to;
				}
				if(dim == 3)
				{
					// Mapping NaturalEdges -> New Edges
					const size_t numNewEdgeInfo = m_vNewEdgeInfo.size();
					m_vNatEdgeInfo[i].nodeId = newNodeId;
					m_vNatEdgeInfo[i].numChildEdges = 2;
					m_vNatEdgeInfo[i].childEdge[0] = numNewEdgeInfo;
					m_vNatEdgeInfo[i].childEdge[1] = numNewEdgeInfo + 1;

					m_vNewEdgeInfo.resize(numNewEdgeInfo + 2);

					m_vNewEdgeInfo[numNewEdgeInfo].m_from = from;
					m_vNewEdgeInfo[numNewEdgeInfo].m_to = newNodeId;

					m_vNewEdgeInfo[numNewEdgeInfo+1].m_from = newNodeId;
					m_vNewEdgeInfo[numNewEdgeInfo+1].m_to = to;
				}
			}
			break;

		default: UG_THROW("Cannot detect type of edge.");
		}
	}

	// for 3d case also check faces for hanging nodes
	if(dim == 3)
	{
		std::vector<Face*> vFaces;
		CollectFacesSorted(vFaces, grid, pElem);

		// compute Nodes
		MathVector<dim> locSideMid;
		MathVector<worldDim> gloSideMid;
		for(size_t i = 0; i < vFaces.size(); ++i)
		{
			///////////
			// case QUADRILATERAL with all edges hanging and hanging node in middle
			///////////
			if(vFaces[i]->container_section() == CSFACE_CONSTRAINING_QUADRILATERAL)
			{
				// insert hanging node in list of nodes
				const size_t newNodeId = m_gloMid[0].size();
				m_gloMid[0].resize(newNodeId + 1);
				m_locMid[0].resize(newNodeId + 1);

				// compute position of new (hanging) node
				compute_side_midpoints(i, m_locMid[0][newNodeId], m_gloMid[0][newNodeId]);

				// loop constrained faces
				for(size_t j = 0; j < m_rRefElem.num(2, i, 1); ++j)
				{
					const size_t jplus1 = (j+1)%4;

					// natural edges
					const size_t natEdId1 = m_rRefElem.id(2, i, 1, j);
					const size_t natEdId2 = m_rRefElem.id(2, i, 1, jplus1);

					// corner of the face
					const size_t cornerId = m_rRefElem.id(2,i, 0, jplus1);

					// refined edges that belong to this face
					const size_t edId1 = get_child_edge_of_corner(natEdId1, cornerId);
					const size_t edId2 = get_child_edge_of_corner(natEdId2, cornerId);

					// nodes of hanging edges
					const size_t hangEdNodeId1 = m_vNatEdgeInfo[natEdId1].node_id();
					const size_t hangEdNodeId2 = m_vNatEdgeInfo[natEdId2].node_id();

					// mid point of hanging side
					compute_side_midpoints(	cornerId, newNodeId,
											hangEdNodeId1, hangEdNodeId2,
											locSideMid, gloSideMid);

					// add side midpoint to already existing scvf of this side
					const size_t numSCVF = m_vSCVF.size();
					m_vSCVF.resize(numSCVF + 4);

					m_vSCVF[numSCVF].m_from = m_vNewEdgeInfo[edId1].from();
					m_vSCVF[numSCVF].m_to = m_vNewEdgeInfo[edId1].to();
					m_vSCVF[numSCVF].m_vLocPos[2] = locSideMid;
					m_vSCVF[numSCVF].m_vGloPos[2] = gloSideMid;

					m_vSCVF[numSCVF+1].m_from = m_vNewEdgeInfo[edId2].from();
					m_vSCVF[numSCVF+1].m_to = m_vNewEdgeInfo[edId2].to();
					m_vSCVF[numSCVF+1].m_vLocPos[2] = locSideMid;
					m_vSCVF[numSCVF+1].m_vGloPos[2] = gloSideMid;

					m_vSCVF[numSCVF+2].m_from = hangEdNodeId1;
					m_vSCVF[numSCVF+2].m_to = newNodeId;
					m_vSCVF[numSCVF+2].m_vLocPos[2] = locSideMid;
					m_vSCVF[numSCVF+2].m_vGloPos[2] = gloSideMid;

					m_vSCVF[numSCVF+3].m_from = hangEdNodeId2;
					m_vSCVF[numSCVF+3].m_to = newNodeId;
					m_vSCVF[numSCVF+3].m_vLocPos[2] = locSideMid;
					m_vSCVF[numSCVF+3].m_vGloPos[2] = gloSideMid;
				}
			}
			///////////
			// case TRIANGLE with all edges hanging, that matches a refined element on other side
			///////////
			else if (vFaces[i]->container_section() == CSFACE_CONSTRAINING_TRIANGLE)
			{
				bool bAllConstraining = true;
				for(size_t j = 0; j < m_rRefElem.num(2, i, 1); ++j)
					if(vEdges[m_rRefElem.id(2, i, 1, j)]->container_section() != CSEDGE_CONSTRAINING_EDGE)
						bAllConstraining = false;

				if(!bAllConstraining) continue;

				// compute position of new (hanging) node
				compute_side_midpoints(i, locSideMid, gloSideMid);

				// loop constrained faces
				for(size_t j = 0; j < m_rRefElem.num(2, i, 1); ++j)
				{
					const size_t jplus1 = (j+1)%3;

					// natural edges
					const size_t natEdId1 = m_rRefElem.id(2, i, 1, j);
					const size_t natEdId2 = m_rRefElem.id(2, i, 1, jplus1);

					// corner of the face
					const size_t cornerId = m_rRefElem.id(2,i, 0, jplus1);

					// nodes of hanging edges
					const size_t hangEdNodeId1 = m_vNatEdgeInfo[natEdId1].node_id();
					const size_t hangEdNodeId2 = m_vNatEdgeInfo[natEdId2].node_id();

					MathVector<dim> locSmallSideMid;
					MathVector<worldDim> gloSmallSideMid;

					// mid point of hanging side
					compute_side_midpoints(	cornerId, hangEdNodeId1, hangEdNodeId2,
											locSmallSideMid, gloSmallSideMid);

					// add side midpoint to already existing scvf of this side
					const size_t numSCVF = m_vSCVF.size();
					m_vSCVF.resize(numSCVF + 4);

					m_vSCVF[numSCVF].m_from = hangEdNodeId1;
					m_vSCVF[numSCVF].m_to = hangEdNodeId2;
					m_vSCVF[numSCVF].m_vLocPos[2] = locSideMid;
					m_vSCVF[numSCVF].m_vGloPos[2] = gloSideMid;

					m_vSCVF[numSCVF+1].m_from = hangEdNodeId1;
					m_vSCVF[numSCVF+1].m_to = hangEdNodeId2;
					m_vSCVF[numSCVF+1].m_vLocPos[2] = locSmallSideMid;
					m_vSCVF[numSCVF+1].m_vGloPos[2] = gloSmallSideMid;

					m_vSCVF[numSCVF+2].m_from = hangEdNodeId1;
					m_vSCVF[numSCVF+2].m_to = cornerId;
					m_vSCVF[numSCVF+2].m_vLocPos[2] = locSmallSideMid;
					m_vSCVF[numSCVF+2].m_vGloPos[2] = gloSmallSideMid;

					m_vSCVF[numSCVF+3].m_from = cornerId;
					m_vSCVF[numSCVF+3].m_to = hangEdNodeId2;
					m_vSCVF[numSCVF+3].m_vLocPos[2] = locSmallSideMid;
					m_vSCVF[numSCVF+3].m_vGloPos[2] = gloSmallSideMid;
				}
			}
			//////////
			// other cases: Not all edges hanging (i.e. neighbor not refined)
			///////////
			else
			{
				// compute side midpoint
				compute_side_midpoints(i, locSideMid, gloSideMid);

				// connect all edges with side midpoint
				for(size_t j = 0; j < m_rRefElem.num(2, i, 1); ++j)
				{
					const size_t natEdgeId = m_rRefElem.id(2, i, 1, j);
					for(size_t e = 0; e < m_vNatEdgeInfo[natEdgeId].num_child_edges(); ++e)
					{
						const size_t edgeId = m_vNatEdgeInfo[natEdgeId].child_edge(e);

						const size_t numSCVF = m_vSCVF.size();
						m_vSCVF.resize(numSCVF + 1);

						m_vSCVF[numSCVF].m_from = m_vNewEdgeInfo[edgeId].from();
						m_vSCVF[numSCVF].m_to = m_vNewEdgeInfo[edgeId].to();
						m_vSCVF[numSCVF].m_vGloPos[2] = gloSideMid;
						m_vSCVF[numSCVF].m_vLocPos[2] = locSideMid;
					}
				}
			}
		}
	}

	// resize number of scv (== number of nodes)
	if(dim <= 2)
	{
		m_vSCV.resize(m_gloMid[0].size());
	}
	else
	{
		m_vSCV.clear();
	}

	// compute scvf
	for(size_t i = 0; i < m_vSCVF.size(); ++i)
	{
		// compute midpoints of edges
		{
			const size_t from = m_vSCVF[i].m_from;
			const size_t to = m_vSCVF[i].m_to;

			VecInterpolateLinear(	m_vSCVF[i].m_vLocPos[0],
									m_locMid[0][from],
									m_locMid[0][to],
									0.5);
			VecInterpolateLinear(	m_vSCVF[i].m_vGloPos[0],
									m_gloMid[0][from],
									m_gloMid[0][to],
									0.5);
		}

		// set center of elem as part of scvf
		m_vSCVF[i].m_vGloPos[1] = m_gloMid[dim][0];
		m_vSCVF[i].m_vLocPos[1] = m_locMid[dim][0];

		// integration point
		AveragePositions(m_vSCVF[i].localIP, m_vSCVF[i].m_vLocPos, SCVF::m_numCorners);
		AveragePositions(m_vSCVF[i].globalIP, m_vSCVF[i].m_vGloPos, SCVF::m_numCorners);

		// normal
		switch (m_roid){
			case ROID_EDGE: HangingNormalOnSCVF<elem_type_0, worldDim>(m_vSCVF[i].Normal, &(m_vSCVF[i].m_vGloPos[0])); break;
			case ROID_TRIANGLE: HangingNormalOnSCVF<elem_type_0, worldDim>(m_vSCVF[i].Normal, &(m_vSCVF[i].m_vGloPos[0])); break;
			case ROID_QUADRILATERAL: HangingNormalOnSCVF<elem_type_1, worldDim>(m_vSCVF[i].Normal, &(m_vSCVF[i].m_vGloPos[0])); break;
			case ROID_TETRAHEDRON: HangingNormalOnSCVF<elem_type_0, worldDim>(m_vSCVF[i].Normal, &(m_vSCVF[i].m_vGloPos[0])); break;
			case ROID_PYRAMID: HangingNormalOnSCVF<elem_type_1, worldDim>(m_vSCVF[i].Normal, &(m_vSCVF[i].m_vGloPos[0])); break;
			case ROID_PRISM: HangingNormalOnSCVF<elem_type_2, worldDim>(m_vSCVF[i].Normal, &(m_vSCVF[i].m_vGloPos[0])); break;
			case ROID_HEXAHEDRON: HangingNormalOnSCVF<elem_type_3, worldDim>(m_vSCVF[i].Normal, &(m_vSCVF[i].m_vGloPos[0])); break;
			default: UG_THROW("unsupported element type"); break;
		}

		// TODO: In 3D check orientation
		if(dim == 3)
		{
			const size_t from = m_vSCVF[i].m_from;
			const size_t to = m_vSCVF[i].m_to;

			MathVector<worldDim> diff;
			VecSubtract(diff, m_gloMid[0][to], m_gloMid[0][from]);

			if(VecDot(diff, m_vSCVF[i].Normal) < 0)
			{
				m_vSCVF[i].m_from = to;
				m_vSCVF[i].m_to = from;
			}
		}

		// write edge midpoints to as corners of scv
		if(dim == 2)
		{
			const size_t from = m_vSCVF[i].m_from;
			const size_t to = m_vSCVF[i].m_to;

			m_vSCV[from].m_vLocPos[1] = m_vSCVF[i].m_vLocPos[0];
			m_vSCV[from].m_vGloPos[1] = m_vSCVF[i].m_vGloPos[0];

			m_vSCV[to].m_vLocPos[3] = m_vSCVF[i].m_vLocPos[0];
			m_vSCV[to].m_vGloPos[3] = m_vSCVF[i].m_vGloPos[0];
		}
	}

	// compute size of scv
	if(dim <= 2)
	{
		for(size_t i = 0; i < m_vSCV.size(); ++i)
		{
			// set node id
			m_vSCV[i].nodeId = i;

			// start at node
			m_vSCV[i].m_vLocPos[0] = m_locMid[0][i];
			m_vSCV[i].m_vGloPos[0] = m_gloMid[0][i];

			if(dim == 1)
			{
				// center of element
				m_vSCV[i].m_vLocPos[1] = m_locMid[dim][0];
				m_vSCV[i].m_vGloPos[1] = m_gloMid[dim][0];
			}
			else if (dim == 2)
			{
				// center of element
				m_vSCV[i].m_vLocPos[2] = m_locMid[dim][0];
				m_vSCV[i].m_vGloPos[2] = m_gloMid[dim][0];
			}

           m_vSCV[i].vol = ElementSize<typename SCV::scv_type, worldDim>(&(m_vSCV[i].m_vGloPos[0]));;
		}
	}

	if(dim == 3)
	{
		m_vSCV.resize(m_vSCVF.size() * 2);

		for(size_t i = 0; i < m_vSCVF.size(); ++i)
        {
			for(size_t n = 0; n < 3; ++n)
			{
				m_vSCV[2*i + 0].m_vLocPos[n+1] = m_vSCVF[i].m_vLocPos[n];
				m_vSCV[2*i + 0].m_vGloPos[n+1] = m_vSCVF[i].m_vGloPos[n];
				m_vSCV[2*i + 1].m_vLocPos[n+1] = m_vSCVF[i].m_vLocPos[n];
				m_vSCV[2*i + 1].m_vGloPos[n+1] = m_vSCVF[i].m_vGloPos[n];
			}
			const size_t from = m_vSCVF[i].m_from;
			const size_t to = m_vSCVF[i].m_to;

			m_vSCV[2*i + 0].m_vLocPos[0] = m_locMid[0][from];
			m_vSCV[2*i + 0].m_vGloPos[0] = m_gloMid[0][from];
			m_vSCV[2*i + 0].nodeId = from;
			m_vSCV[2*i + 0].m_numCorners = 4;

			m_vSCV[2*i + 1].m_vLocPos[0] = m_locMid[0][to];
			m_vSCV[2*i + 1].m_vGloPos[0] = m_gloMid[0][to];
			m_vSCV[2*i + 1].nodeId = to;
			m_vSCV[2*i + 1].m_numCorners = 4;

			m_vSCV[2*i + 0].vol = ElementSize<typename SCV::scv_type, worldDim>(&(m_vSCV[2*i + 0].m_vGloPos[0]));;
			m_vSCV[2*i + 1].vol = ElementSize<typename SCV::scv_type, worldDim>(&(m_vSCV[2*i + 1].m_vGloPos[0]));;
		}
	}

	/////////////////////////
	// Shapes and Derivatives
	/////////////////////////
	m_rMapping->update(vCornerCoords);
	
	const LocalShapeFunctionSet<dim>& TrialSpace =
		LocalFiniteElementProvider::get<dim>(m_roid, LFEID(LFEID::LAGRANGE, dim, 1));

	const size_t num_sh = TrialSpace.num_sh();
	m_numSh = num_sh;

	for(size_t i = 0; i < num_scvf(); ++i)
	{
		m_rMapping->jacobian_transposed_inverse(m_vSCVF[i].JtInv, m_vSCVF[i].localIP);
		m_vSCVF[i].detj = m_rMapping->sqrt_gram_det(m_vSCVF[i].localIP);

		m_vSCVF[i].vShape.resize(num_sh);
		m_vSCVF[i].localGrad.resize(num_sh);
		m_vSCVF[i].globalGrad.resize(num_sh);

		TrialSpace.shapes(&(m_vSCVF[i].vShape[0]), m_vSCVF[i].localIP);
		TrialSpace.grads(&(m_vSCVF[i].localGrad[0]), m_vSCVF[i].localIP);

		for(size_t sh = 0 ; sh < num_sh; ++sh)
		{
			MatVecMult((m_vSCVF[i].globalGrad)[sh], m_vSCVF[i].JtInv, (m_vSCVF[i].localGrad)[sh]);
		}

	}
	
	for(size_t i = 0; i < num_scv(); ++i)
	{
		m_rMapping->jacobian_transposed_inverse(m_vSCV[i].JtInv, m_vSCVF[i].m_vLocPos[0]);
		m_vSCVF[i].detj = m_rMapping->sqrt_gram_det(m_vSCV[i].m_vLocPos[0]);

		TrialSpace.shapes(&(m_vSCV[i].vShape[0]), m_vSCV[i].m_vLocPos[0]);
		TrialSpace.grads(&(m_vSCV[i].localGrad[0]), m_vSCV[i].m_vLocPos[0]);

		for(size_t sh = 0 ; sh < num_sh; ++sh)
		{
			MatVecMult((m_vSCV[i].globalGrad)[sh], m_vSCV[i].JtInv, (m_vSCV[i].localGrad)[sh]);
		}
	}

	///////////////////////////
	// Copy ip pos in list
	///////////////////////////

// 	loop Sub Control Volumes (SCV)
	m_vGlobSCVIP.clear();
	m_vLocSCVIP.clear();
	for(size_t i = 0; i < num_scv(); ++i)
	{
	//	get current SCV
		const SCV& rSCV = scv(i);

	//	copy
		m_vGlobSCVIP.push_back(rSCV.global_ip());
		m_vLocSCVIP.push_back(rSCV.local_ip());
	}

// 	loop Sub Control Volumes Faces (SCVF)
	m_vGlobSCVFIP.clear();
	m_vLocSCVFIP.clear();
	for(size_t i = 0; i < num_scvf(); ++i)
	{
	//	get current SCVF
		const SCVF& rSCVF = scvf(i);

	//  copy
		m_vGlobSCVFIP.push_back(rSCVF.global_ip());
		m_vLocSCVFIP.push_back(rSCVF.local_ip());
	}
	
}



////////////////////////////////////////////////////////////////////////////////
// HFV1ManifoldGeometry
////////////////////////////////////////////////////////////////////////////////

template <typename TElem, int TWorldDim>
HFV1ManifoldGeometry<TElem, TWorldDim>::
HFV1ManifoldGeometry() : m_pElem(NULL), m_rRefElem(Provider<ref_elem_type>::get()), m_ssi(-1)
{
	// set corners of element as local centers of nodes
	m_vBF.resize(m_numNaturalBF);
	m_locMid[0].resize(m_numNaturalBF);
	for (size_t i = 0; i < m_rRefElem.num(0); ++i)
		m_locMid[0][i] = m_rRefElem.corner(i);

	// compute local elem center
	m_locMid[dim].resize(1);
	m_gloMid[dim].resize(1);
	m_locMid[dim][0] = 0.0;
	for (size_t i = 0; i < m_locMid[0].size(); ++i)
	{
		m_locMid[dim][0] += m_locMid[0][i];
	}
	m_locMid[dim][0] *= 1.0/(m_locMid[0].size());
}


/// update data for given element
template <typename TElem, int TWorldDim>
void HFV1ManifoldGeometry<TElem, TWorldDim>::
update(GridObject* elem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
	/* FIXME: This is ugly. We need this distinction, because elem might be of type
	 * ConstrainingTriangle or ConstrainingQuadrilateral. Unfortunately, TElem can only be
	 * of type Triangle or Quadrilateral, even though the Constraining* variants work
	 * with the exact same code.
	 * Best would be for Constraining* to inherit from * (as is the case for ConstrainingEdge).
	 */
	if (dim == 1)
	{
		UG_ASSERT(dynamic_cast<Edge*>(elem), "Wrong element type.");
	}
	else if (dim == 2)
	{
		UG_ASSERT(dynamic_cast<Face*>(elem), "Wrong element type.");
	}

	// 	if already update for this element, do nothing
	// we can do this only because we do not use m_pElem elsewhere,
	// otherwise this cast would be very dangerous!
	if (m_pElem == reinterpret_cast<TElem*>(elem)) return;
	else m_pElem = reinterpret_cast<TElem*>(elem);

	//	store subset index for geometry
	Grid& grid = *(ish->grid());
	if (dim == 1)
	{
		std::vector<Edge*> vEdge;
		CollectEdges(vEdge, grid, static_cast<Edge*>(elem));
		UG_ASSERT(vEdge.size(),"No edge contained in 1D manifold element!");
		m_ssi = ish->get_subset_index(vEdge[0]);
	}
	if (dim == 2)
	{
		std::vector<Face*> vFace;
		CollectFaces(vFace, grid, static_cast<Face*>(elem));
		UG_ASSERT(vFace.size(),"No face contained in 2D manifold element!");
		m_ssi = ish->get_subset_index(vFace[0]);
	}

	// reset to natural nodes
	m_gloMid[0].resize(m_numNaturalBF);
	m_locMid[0].resize(m_numNaturalBF);

	// remember global position of nodes
	for (size_t i = 0; i < m_numNaturalBF; i++)
		m_gloMid[0][i] = vCornerCoords[i];

	// compute global center
	m_gloMid[dim][0] = 0.0;
	for (size_t i = 0; i < m_gloMid[0].size(); i++)
		m_gloMid[dim][0] += m_gloMid[0][i];

	m_gloMid[dim][0] *= 1./(m_gloMid[0].size());

	// get natural edges
	std::vector<Edge*> vEdges;
	CollectEdgesSorted(vEdges, grid, elem);

	// compute nodes
	m_vBF.clear();
	m_vNewEdgeInfo.clear();
	m_vNatEdgeInfo.clear(); m_vNatEdgeInfo.resize(m_numNaturalBFS);
	UG_ASSERT(vEdges.size() == m_numNaturalBFS, "Incorrect number of edges found, only " << vEdges.size() << "Edges");

	bool bAllConstraining = true;
	for (size_t i = 0; i < vEdges.size(); i++)
	{
		// natural ids of end of edge
		const size_t from = m_rRefElem.id(1, i, 0, 0);
		const size_t to = m_rRefElem.id(1, i, 0, 1);

		// choose whether to insert two or one new edge
		switch (vEdges[i]->container_section())
		{
			case CSEDGE_CONSTRAINED_EDGE:
			// classic case (no hanging node on edge)
			case CSEDGE_REGULAR_EDGE:
				if (dim == 1)
				{
					// set loc & glo corner coords
					const size_t numBF = m_vBF.size();
					m_vBF.resize(numBF + 2);
					MathVector<dim> locMidPt;
					MathVector<worldDim> gloMidPt;
					VecInterpolateLinear(locMidPt, m_locMid[0][from], m_locMid[0][to], 0.5);
					VecInterpolateLinear(gloMidPt, m_gloMid[0][from], m_gloMid[0][to], 0.5);

					m_vBF[numBF].m_vLocPos[0] = m_locMid[0][from];
					m_vBF[numBF].m_vLocPos[1] = locMidPt;
					m_vBF[numBF].m_vGloPos[0] = m_gloMid[0][from];
					m_vBF[numBF].m_vGloPos[1] = gloMidPt;
					m_vBF[numBF].nodeId = from;

					m_vBF[numBF+1].m_vLocPos[0] = m_locMid[0][to];
					m_vBF[numBF+1].m_vLocPos[1] = locMidPt;
					m_vBF[numBF+1].m_vGloPos[0] = m_gloMid[0][to];
					m_vBF[numBF+1].m_vGloPos[1] = gloMidPt;
					m_vBF[numBF+1].nodeId = to;
				}
				if (dim == 2)
				{
					// remember this edge (and that it is not constraining)
					const size_t numNewEdgeInfo = m_vNewEdgeInfo.size();
					m_vNatEdgeInfo[i].numChildEdges = 1;
					m_vNatEdgeInfo[i].childEdge[0] = numNewEdgeInfo;

					m_vNewEdgeInfo.resize(numNewEdgeInfo + 1);
					m_vNewEdgeInfo[numNewEdgeInfo].m_from = from;
					m_vNewEdgeInfo[numNewEdgeInfo].m_to = to;

					bAllConstraining = false;
				}
				break;

			// hanging node case
			case CSEDGE_CONSTRAINING_EDGE:
				{
					// insert hanging node in list of nodes
					const size_t newNodeId = m_gloMid[0].size();
					m_gloMid[0].resize(newNodeId + 1);
					m_locMid[0].resize(newNodeId + 1);
					VecInterpolateLinear(	m_gloMid[0].back(),
											m_gloMid[0][from],
											m_gloMid[0][to],
											0.5);
					VecInterpolateLinear(	m_locMid[0].back(),
											m_locMid[0][from],
											m_locMid[0][to],
											0.5);

					if (dim == 1)	// cannot happen! (loop is only over RegularEdge for dim==1!)
					{
						// insert two edges with nodeIds and set loc & glo corner coords
						const size_t numBF = m_vBF.size();
						m_vBF.resize(numBF + 4);

						MathVector<dim> locMidPt;
						MathVector<worldDim> gloMidPt;
						VecInterpolateLinear(locMidPt, m_locMid[0][from], m_locMid[0][newNodeId], 0.5);
						VecInterpolateLinear(gloMidPt, m_gloMid[0][from], m_gloMid[0][newNodeId], 0.5);

						m_vBF[numBF].m_vLocPos[0] = m_locMid[0][from];
						m_vBF[numBF].m_vLocPos[1] = locMidPt;
						m_vBF[numBF].m_vGloPos[0] = m_gloMid[0][from];
						m_vBF[numBF].m_vGloPos[1] = gloMidPt;
						m_vBF[numBF].nodeId = from;

						m_vBF[numBF+1].m_vLocPos[0] = m_locMid[0][newNodeId];
						m_vBF[numBF+1].m_vLocPos[1] = locMidPt;
						m_vBF[numBF+1].m_vGloPos[0] = m_gloMid[0][newNodeId];
						m_vBF[numBF+1].m_vGloPos[1] = gloMidPt;
						m_vBF[numBF+1].nodeId = newNodeId;

						VecInterpolateLinear(locMidPt, m_locMid[0][to], m_locMid[0][newNodeId], 0.5);
						VecInterpolateLinear(gloMidPt, m_gloMid[0][to], m_gloMid[0][newNodeId], 0.5);

						m_vBF[numBF+2].m_vLocPos[0] = m_locMid[0][to];
						m_vBF[numBF+2].m_vLocPos[1] = locMidPt;
						m_vBF[numBF+2].m_vGloPos[0] = m_gloMid[0][to];
						m_vBF[numBF+2].m_vGloPos[1] = gloMidPt;
						m_vBF[numBF+2].nodeId = to;

						m_vBF[numBF+3].m_vLocPos[0] = m_locMid[0][newNodeId];
						m_vBF[numBF+3].m_vLocPos[1] = locMidPt;
						m_vBF[numBF+3].m_vGloPos[0] = m_gloMid[0][newNodeId];
						m_vBF[numBF+3].m_vGloPos[1] = gloMidPt;
						m_vBF[numBF+3].nodeId = newNodeId;
					}

					if (dim == 2)
					{
						// mapping naturalEdges -> new edges
						const size_t numNewEdgeInfo = m_vNewEdgeInfo.size();
						m_vNatEdgeInfo[i].nodeId = newNodeId;
						m_vNatEdgeInfo[i].numChildEdges = 2;
						m_vNatEdgeInfo[i].childEdge[0] = numNewEdgeInfo;
						m_vNatEdgeInfo[i].childEdge[1] = numNewEdgeInfo + 1;

						m_vNewEdgeInfo.resize(numNewEdgeInfo + 2);

						m_vNewEdgeInfo[numNewEdgeInfo].m_from = from;
						m_vNewEdgeInfo[numNewEdgeInfo].m_to = newNodeId;

						m_vNewEdgeInfo[numNewEdgeInfo+1].m_from = newNodeId;
						m_vNewEdgeInfo[numNewEdgeInfo+1].m_to = to;
					}
				}
				break;

			default: UG_THROW("Cannot detect type of edge.");
		}
	}

	// for 3d case hanging fv1 manifold geometry depends on whether or not it is induced by the
	// hanging fv1 geometry of a refined element or not
	if (dim == 2)
	{
		std::vector<Face*> vFaces;
		CollectFacesSorted(vFaces, grid, elem);

		if (vFaces.size() != 1)
			UG_THROW("Expected number of faces to be exactly 1, but number is " << vFaces.size() << ".");

		Face* face = vFaces[0];

		// compute Nodes
		MathVector<dim> locSideMid;
		MathVector<worldDim> gloSideMid;

		///////////
		// case QUADRILATERAL with all edges hanging and hanging node in middle
		//
		// cannot happen, since only Quadrilaterals are looped, not ConstrainingQuadrilaterals!
		///////////
		if (face->container_section() == CSFACE_CONSTRAINING_QUADRILATERAL)
		{
			// insert hanging node in list of nodes
			const size_t newNodeId = m_gloMid[0].size();
			m_gloMid[0].resize(newNodeId + 1);
			m_locMid[0].resize(newNodeId + 1);

			// compute position of new (hanging) node
			compute_side_midpoints(m_locMid[0][newNodeId], m_gloMid[0][newNodeId]);

			// loop constrained edges
			for (size_t j = 0; j < m_rRefElem.num(2, 0, 1); ++j)
			{
				const size_t jplus1 = (j+1)%4;

				// natural edges
				const size_t natEdId1 = m_rRefElem.id(2, 0, 1, j);
				const size_t natEdId2 = m_rRefElem.id(2, 0, 1, jplus1);

				// corner between the two edges
				const size_t cornerId = m_rRefElem.id(2, 0, 0, jplus1);

				// nodes of hanging edges
				const size_t hangEdNodeId1 = m_vNatEdgeInfo[natEdId1].node_id();
				const size_t hangEdNodeId2 = m_vNatEdgeInfo[natEdId2].node_id();

				// mid point of hanging side
				compute_side_midpoints(	cornerId, newNodeId,
										hangEdNodeId1, hangEdNodeId2,
										locSideMid, gloSideMid);

				// create corresponding boundary faces
				const size_t numBF = m_vBF.size();
				m_vBF.resize(numBF + 4);

				MathVector<dim> locEdgeMid0, locEdgeMid1, locEdgeMid2, locEdgeMid3;
				MathVector<worldDim> gloEdgeMid0, gloEdgeMid1, gloEdgeMid2, gloEdgeMid3;
				VecInterpolateLinear(locEdgeMid0, m_locMid[0][cornerId], m_locMid[0][hangEdNodeId2], 0.5);
				VecInterpolateLinear(gloEdgeMid0, m_gloMid[0][cornerId], m_gloMid[0][hangEdNodeId2], 0.5);
				VecInterpolateLinear(locEdgeMid1, m_locMid[0][hangEdNodeId2], m_locMid[0][newNodeId], 0.5);
				VecInterpolateLinear(gloEdgeMid1, m_gloMid[0][hangEdNodeId2], m_gloMid[0][newNodeId], 0.5);
				VecInterpolateLinear(locEdgeMid2, m_locMid[0][newNodeId], m_locMid[0][hangEdNodeId1], 0.5);
				VecInterpolateLinear(gloEdgeMid2, m_gloMid[0][newNodeId], m_gloMid[0][hangEdNodeId1], 0.5);
				VecInterpolateLinear(locEdgeMid3, m_locMid[0][hangEdNodeId1], m_locMid[0][cornerId], 0.5);
				VecInterpolateLinear(gloEdgeMid3, m_gloMid[0][hangEdNodeId1], m_gloMid[0][cornerId], 0.5);

				m_vBF[numBF].m_vLocPos[0] = m_locMid[0][cornerId];
				m_vBF[numBF].m_vLocPos[1] = locEdgeMid0;
				m_vBF[numBF].m_vLocPos[2] = locSideMid;
				m_vBF[numBF].m_vLocPos[3] = locEdgeMid3;
				m_vBF[numBF].m_vGloPos[0] = m_gloMid[0][cornerId];
				m_vBF[numBF].m_vGloPos[1] = gloEdgeMid0;
				m_vBF[numBF].m_vGloPos[2] = gloSideMid;
				m_vBF[numBF].m_vGloPos[3] = gloEdgeMid3;
				m_vBF[numBF].nodeId = cornerId;

				m_vBF[numBF+1].m_vLocPos[0] = m_locMid[0][hangEdNodeId2];
				m_vBF[numBF+1].m_vLocPos[1] = locEdgeMid1;
				m_vBF[numBF+1].m_vLocPos[2] = locSideMid;
				m_vBF[numBF+1].m_vLocPos[3] = locEdgeMid0;
				m_vBF[numBF+1].m_vGloPos[0] = m_gloMid[0][hangEdNodeId2];
				m_vBF[numBF+1].m_vGloPos[1] = gloEdgeMid1;
				m_vBF[numBF+1].m_vGloPos[2] = gloSideMid;
				m_vBF[numBF+1].m_vGloPos[3] = gloEdgeMid0;
				m_vBF[numBF+1].nodeId = hangEdNodeId2;

				m_vBF[numBF+2].m_vLocPos[0] = m_locMid[0][newNodeId];
				m_vBF[numBF+2].m_vLocPos[1] = locEdgeMid2;
				m_vBF[numBF+2].m_vLocPos[2] = locSideMid;
				m_vBF[numBF+2].m_vLocPos[3] = locEdgeMid1;
				m_vBF[numBF+2].m_vGloPos[0] = m_gloMid[0][newNodeId];
				m_vBF[numBF+2].m_vGloPos[1] = gloEdgeMid2;
				m_vBF[numBF+2].m_vGloPos[2] = gloSideMid;
				m_vBF[numBF+2].m_vGloPos[3] = gloEdgeMid1;
				m_vBF[numBF+2].nodeId = newNodeId;

				m_vBF[numBF+3].m_vLocPos[0] = m_locMid[0][hangEdNodeId1];
				m_vBF[numBF+3].m_vLocPos[1] = locEdgeMid3;
				m_vBF[numBF+3].m_vLocPos[2] = locSideMid;
				m_vBF[numBF+3].m_vLocPos[3] = locEdgeMid2;
				m_vBF[numBF+3].m_vGloPos[0] = m_gloMid[0][hangEdNodeId1];
				m_vBF[numBF+3].m_vGloPos[1] = gloEdgeMid3;
				m_vBF[numBF+3].m_vGloPos[2] = gloSideMid;
				m_vBF[numBF+3].m_vGloPos[3] = gloEdgeMid2;
				m_vBF[numBF+3].nodeId = hangEdNodeId1;
			}
		}
		///////////
		// case TRIANGLE with all edges hanging, that matches a refined element on other side
		//
		// cannot happen, since only Triangles are looped, not ConstrainingTriangles!
		///////////
		else if (bAllConstraining && face->container_section() == CSFACE_CONSTRAINING_TRIANGLE)
		{
			// compute position of side mid
			compute_side_midpoints(locSideMid, gloSideMid);

			// loop constrained edges
			for (size_t j = 0; j < m_rRefElem.num(2, 0, 1); ++j)
			{
				const size_t jplus1 = (j+1)%3;

				// natural edges
				const size_t natEdId1 = m_rRefElem.id(2, 0, 1, j);
				const size_t natEdId2 = m_rRefElem.id(2, 0, 1, jplus1);

				// corner between the two edges
				const size_t cornerId = m_rRefElem.id(2, 0, 0, jplus1);

				// nodes of hanging edges
				const size_t hangEdNodeId1 = m_vNatEdgeInfo[natEdId1].node_id();
				const size_t hangEdNodeId2 = m_vNatEdgeInfo[natEdId2].node_id();

				MathVector<dim> locSmallSideMid;
				MathVector<worldDim> gloSmallSideMid;

				// mid point of hanging side
				compute_side_midpoints(	cornerId, hangEdNodeId1, hangEdNodeId2,
										locSmallSideMid, gloSmallSideMid);

				MathVector<dim> locEdgeMid0, locEdgeMid1;
				MathVector<worldDim> gloEdgeMid0, gloEdgeMid1;
				VecInterpolateLinear(locEdgeMid0, m_locMid[0][cornerId], m_locMid[0][hangEdNodeId2], 0.5);
				VecInterpolateLinear(gloEdgeMid0, m_gloMid[0][cornerId], m_gloMid[0][hangEdNodeId2], 0.5);
				VecInterpolateLinear(locEdgeMid1, m_locMid[0][hangEdNodeId1], m_locMid[0][cornerId], 0.5);
				VecInterpolateLinear(gloEdgeMid1, m_gloMid[0][hangEdNodeId1], m_gloMid[0][cornerId], 0.5);

				// create corresponding boundary faces
				const size_t numBF = m_vBF.size();
				m_vBF.resize(numBF + 3);

				m_vBF[numBF].m_vLocPos[0] = m_locMid[0][cornerId];
				m_vBF[numBF].m_vLocPos[1] = locEdgeMid0;
				m_vBF[numBF].m_vLocPos[2] = locSmallSideMid;
				m_vBF[numBF].m_vLocPos[3] = locEdgeMid1;
				m_vBF[numBF].m_vGloPos[0] = m_gloMid[0][cornerId];
				m_vBF[numBF].m_vGloPos[1] = gloEdgeMid0;
				m_vBF[numBF].m_vGloPos[2] = gloSmallSideMid;
				m_vBF[numBF].m_vGloPos[3] = gloEdgeMid1;
				m_vBF[numBF].nodeId = cornerId;

				m_vBF[numBF+1].m_vLocPos[0] = m_locMid[0][hangEdNodeId2];
				m_vBF[numBF+1].m_vLocPos[1] = locSideMid;
				m_vBF[numBF+1].m_vLocPos[2] = locSmallSideMid;
				m_vBF[numBF+1].m_vLocPos[3] = locEdgeMid0;
				m_vBF[numBF+1].m_vGloPos[0] = m_gloMid[0][hangEdNodeId2];
				m_vBF[numBF+1].m_vGloPos[1] = gloSideMid;
				m_vBF[numBF+1].m_vGloPos[2] = gloSmallSideMid;
				m_vBF[numBF+1].m_vGloPos[3] = gloEdgeMid0;
				m_vBF[numBF+1].nodeId = hangEdNodeId2;

				m_vBF[numBF+2].m_vLocPos[0] = m_locMid[0][hangEdNodeId1];
				m_vBF[numBF+2].m_vLocPos[1] = locEdgeMid0;
				m_vBF[numBF+2].m_vLocPos[2] = locSmallSideMid;
				m_vBF[numBF+2].m_vLocPos[3] = locSideMid;
				m_vBF[numBF+2].m_vGloPos[0] = m_gloMid[0][hangEdNodeId1];
				m_vBF[numBF+2].m_vGloPos[1] = gloEdgeMid0;
				m_vBF[numBF+2].m_vGloPos[2] = gloSmallSideMid;
				m_vBF[numBF+2].m_vGloPos[3] = gloSideMid;
				m_vBF[numBF+2].nodeId = hangEdNodeId1;
			}
		}
		//////////
		// other cases: neighbor not refined, but still hanging nodes
		///////////
		else
		{
			// compute position of side mid
			compute_side_midpoints(locSideMid, gloSideMid);

			// loop edges
			for (size_t j = 0; j < m_rRefElem.num(2, 0, 1); ++j)
			{
				const size_t jplus1 = (j+1) % m_rRefElem.num(2, 0, 1);
				const size_t jplus2 = (j+2) % m_rRefElem.num(2, 0, 1);

				// natural edges
				const size_t natEdId1 = m_rRefElem.id(2, 0, 1, j);
				const size_t natEdId2 = m_rRefElem.id(2, 0, 1, jplus1);

				// corner between the two edges
				const size_t cornerId = m_rRefElem.id(2, 0, 0, jplus1);
				const size_t prevCornerId = m_rRefElem.id(2, 0, 0, j);
				const size_t nextCornerId = m_rRefElem.id(2, 0, 0, jplus2);

				// mid points of edges
				MathVector<dim> locEdgeMidNext, locEdgeMidCurr;
				MathVector<worldDim> gloEdgeMidNext, gloEdgeMidCurr;

				// next side is hanging
				if (m_vNatEdgeInfo[natEdId2].num_child_edges() == 2)
				{
					const size_t hangEdNodeId2 = m_vNatEdgeInfo[natEdId2].node_id();
					VecInterpolateLinear(locEdgeMidNext, m_locMid[0][cornerId], m_locMid[0][hangEdNodeId2], 0.5);
					VecInterpolateLinear(gloEdgeMidNext, m_gloMid[0][cornerId], m_gloMid[0][hangEdNodeId2], 0.5);
				}
				// next side is not hanging
				else
				{
					VecInterpolateLinear(locEdgeMidNext, m_locMid[0][cornerId], m_locMid[0][nextCornerId], 0.5);
					VecInterpolateLinear(gloEdgeMidNext, m_gloMid[0][cornerId], m_gloMid[0][nextCornerId], 0.5);
				}

				// current side is hanging
				if (m_vNatEdgeInfo[natEdId1].num_child_edges() == 2)
				{
					const size_t hangEdNodeId1 = m_vNatEdgeInfo[natEdId1].node_id();

					MathVector<dim> locEdgeMidHang;
					MathVector<worldDim> gloEdgeMidHang;
					VecInterpolateLinear(locEdgeMidCurr, m_locMid[0][cornerId], m_locMid[0][hangEdNodeId1], 0.5);
					VecInterpolateLinear(gloEdgeMidCurr, m_gloMid[0][cornerId], m_gloMid[0][hangEdNodeId1], 0.5);
					VecInterpolateLinear(locEdgeMidHang, m_locMid[0][prevCornerId], m_locMid[0][hangEdNodeId1], 0.5);
					VecInterpolateLinear(gloEdgeMidHang, m_gloMid[0][prevCornerId], m_gloMid[0][hangEdNodeId1], 0.5);

					const size_t numBF = m_vBF.size();
					m_vBF.resize(numBF + 1);

					// create BF for hanging node
					m_vBF[numBF].m_vLocPos[0] = m_locMid[0][hangEdNodeId1];
					m_vBF[numBF].m_vLocPos[1] = locEdgeMidCurr;
					m_vBF[numBF].m_vLocPos[2] = locSideMid;
					m_vBF[numBF].m_vLocPos[3] = locEdgeMidHang;
					m_vBF[numBF].m_vGloPos[0] = m_gloMid[0][hangEdNodeId1];
					m_vBF[numBF].m_vGloPos[1] = gloEdgeMidCurr;
					m_vBF[numBF].m_vGloPos[2] = gloSideMid;
					m_vBF[numBF].m_vGloPos[3] = gloEdgeMidHang;
					m_vBF[numBF].nodeId = hangEdNodeId1;
				}
				// current side is not hanging
				else
				{
					VecInterpolateLinear(locEdgeMidCurr, m_locMid[0][cornerId], m_locMid[0][prevCornerId], 0.5);
					VecInterpolateLinear(gloEdgeMidCurr, m_gloMid[0][cornerId], m_gloMid[0][prevCornerId], 0.5);
				}

				const size_t numBF = m_vBF.size();
				m_vBF.resize(numBF + 1);

				// create BF for corner
				m_vBF[numBF].m_vLocPos[0] = m_locMid[0][cornerId];
				m_vBF[numBF].m_vLocPos[1] = locEdgeMidNext;
				m_vBF[numBF].m_vLocPos[2] = locSideMid;
				m_vBF[numBF].m_vLocPos[3] = locEdgeMidCurr;
				m_vBF[numBF].m_vGloPos[0] = m_gloMid[0][cornerId];
				m_vBF[numBF].m_vGloPos[1] = gloEdgeMidNext;
				m_vBF[numBF].m_vGloPos[2] = gloSideMid;
				m_vBF[numBF].m_vGloPos[3] = gloEdgeMidCurr;
				m_vBF[numBF].nodeId = cornerId;
			}
		}
	}

	// compute size of bf
	for (size_t i = 0; i < m_vBF.size(); ++i)
		m_vBF[i].vol = ElementSize<bf_type, worldDim>(&(m_vBF[i].m_vGloPos[0]));

	// set integration points
	// TODO: The calculated IPs (in the bf midpoints) and shapes are never used
	// as local_ip() and global_ip() always returns the corresponding node coords.
	// Either change local_ip() and global_ip() or remove the members localIP and globalIP
	// and calculate accordingly.
	for (size_t i = 0; i < num_bf(); ++i)
	{
		if (dim == 1) // RegularEdge
		{
			AveragePositions(m_vBF[i].localIP, m_vBF[i].m_vLocPos, BF::numCorners);
			AveragePositions(m_vBF[i].globalIP, m_vBF[i].m_vGloPos, BF::numCorners);
		}
		else if (dim == 2)	// Quadrilateral
		{
			AveragePositions(m_vBF[i].localIP, m_vBF[i].m_vLocPos, BF::numCorners);
			AveragePositions(m_vBF[i].globalIP, m_vBF[i].m_vGloPos, BF::numCorners);
		}
		else {UG_ASSERT(0, "Dimension higher than 2 not implemented.");}
	}

	// shapes
	// A word of warning: This is only meaningful,
	// if the trial space is piecewise linear on tetrahedra/triangles!
	for (size_t i = 0; i < num_bf(); ++i)
	{
		const LocalShapeFunctionSet<ref_elem_type::dim>& TrialSpace =
				LocalFiniteElementProvider::
					get<ref_elem_type::dim>
					(ref_elem_type::REFERENCE_OBJECT_ID,
					 LFEID(LFEID::LAGRANGE, ref_elem_type::dim, 1));

		const size_t num_sh = num_bf();
		m_vBF[i].vShape.resize(num_sh);

		TrialSpace.shapes(&(m_vBF[i].vShape[0]), m_vBF[i].localIP);
	}

	// copy ip positions in list
	m_vGlobBFIP.clear();
	for (size_t i = 0; i < num_bf(); i++)
	{
		const BF& rBF = bf(i);

		m_vLocBFIP.push_back(rBF.local_ip());
		m_vGlobBFIP.push_back(rBF.global_ip());
	}

	//print();
}


/// debug output
template <typename TElem, int TWorldDim>
void HFV1ManifoldGeometry<TElem, TWorldDim>::print()
{
	UG_LOG("\nHanging FVG debug output:\n");
	for (size_t i = 0; i < m_vBF.size(); ++i)
	{
		UG_LOG("SCV " << i << ": ");
		UG_LOG("node_id=" << m_vBF[i].node_id());
		UG_LOG(", local_ip="<< m_vBF[i].local_ip());
		UG_LOG(", global_ip="<< m_vBF[i].global_ip());
		UG_LOG(", vol=" << m_vBF[i].volume());
		UG_LOG("\n");
		for (size_t j = 0; j < m_vBF[i].num_corners(); j++)
			UG_LOG("    localCorner " << j << "=" << m_vBF[i].m_vLocPos[j]);
		UG_LOG("\n");
		for (size_t j = 0; j < m_vBF[i].num_corners(); j++)
			UG_LOG("    globalCorner " << j << "=" << m_vBF[i].m_vGloPos[j]);
		UG_LOG("\n");
	}
	UG_LOG("\n");
	/*
	for(size_t i = 0; i < m_vSCVF.size(); ++i)
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
	*/
}


template class HFV1Geometry<RegularEdge, 1>;
template class HFV1Geometry<RegularEdge, 2>;
template class HFV1Geometry<RegularEdge, 3>;

template class HFV1Geometry<Triangle, 2>;
template class HFV1Geometry<Triangle, 3>;

template class HFV1Geometry<Quadrilateral, 2>;
template class HFV1Geometry<Quadrilateral, 3>;

template class HFV1Geometry<Tetrahedron, 3>;
template class HFV1Geometry<Prism, 3>;
template class HFV1Geometry<Pyramid, 3>;
template class HFV1Geometry<Hexahedron, 3>;

//////////////////////
// DimHFV1Geometry
//////////////////////
template class DimHFV1Geometry<1, 1>;
template class DimHFV1Geometry<1, 2>;
template class DimHFV1Geometry<1, 3>;

template class DimHFV1Geometry<2, 2>;
template class DimHFV1Geometry<2, 3>;

template class DimHFV1Geometry<3, 3>;

//////////////////////
// Manifold
//////////////////////
template class HFV1ManifoldGeometry<Edge, 2>;
template class HFV1ManifoldGeometry<Triangle, 3>;
template class HFV1ManifoldGeometry<Quadrilateral, 3>;


} // end namespace ug

