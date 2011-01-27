/*
 * hanging_finite_volume_geometry.h
 *
 *  Created on: 08.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__HANGING_FINITE_VOLUME_GEOMETRY__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__HANGING_FINITE_VOLUME_GEOMETRY__

// extern libraries
#include <cmath>
#include <vector>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lib_grid.h"

// library intern includes
#include "../../reference_element/reference_element.h"
#include "../../local_shape_function_set/local_shape_function_set_provider.h"
#include "./finite_volume_util.h"

namespace ug{

template <	typename TElem,
			int TWorldDim>
class HFV1Geometry {
	private:
		// type of reference element
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

		// number of SubControlVolumes
		static const size_t m_numNaturalSCV = ref_elem_type::num_corners;

		// number of SubControlVolumeFaces
		static const size_t m_numNaturalSCVF = ref_elem_type::num_edges;

	public:
		// dimension of reference element
		static const int dim = ref_elem_type::dim;

		// dimension of world
		static const int world_dim = TWorldDim;

		// Hanging node flag: this Geometry does support hanging nodes
		static const bool usesHangingNodes = true;

	protected:
		struct MidID
		{
				MidID() : dim(0), id(0) {};
				MidID(size_t dim_, size_t id_) : dim(dim_), id(id_) {};
				size_t dim;
				size_t id;
		};

	public:
		class SCVF
		{
			public:
				// dimension of reference element
				static const int dim = ref_elem_type::dim;

				// dimension of world
				static const int world_dim = TWorldDim;

			private:
				// let outer class access private members
				friend class HFV1Geometry<TElem, TWorldDim>;

				// number of integration points
				static const size_t m_numIP = 1;

				// Number of corners of scvf
				static const size_t m_numCorners = dim;

			public:
				SCVF() {};

				/// index of SubControlVolume on one side of the scvf
				inline size_t from() const {return m_from;}

				/// index of SubControlVolume on one side of the scvf
				inline size_t to() const {return m_to;}

				/// number of integration points on scvf
				inline size_t num_ip() const {return m_numIP;}

				/// local integration point of scvf
				inline const MathVector<dim>& local_ip(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return localIP;}

				/// global integration point of scvf
				inline const MathVector<world_dim>& global_ip(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return globalIP;}

				/// normal on scvf (points direction "from"->"to"). Norm is equal to area
				inline const MathVector<world_dim>& normal() const {return Normal;} // includes area

				/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<dim,world_dim>& JTInv(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return JtInv;}

				/// Determinante of Jacobian in integration point
				inline number detJ(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return detj;}

				/// number of shape functions
				inline size_t num_sh() const {return vShape.size();}

				/// value of shape function i in integration point
				inline number shape(size_t i, size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return vShape[i];}

				/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t i, size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return localGrad[i];}

				/// value of global gradient of shape function i in integration point
				inline const MathVector<world_dim>& global_grad(size_t i, size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return globalGrad[i];}

				/// number of corners, that bound the scvf
				inline size_t num_corners() const {return m_numCorners;}

				/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t i) const
					{UG_ASSERT(i < num_corners(), "Invalid corner index."); return m_vLocPos[i];}

				/// return glbal corner number i
				inline const MathVector<world_dim>& global_corner(size_t i) const
					{UG_ASSERT(i < num_corners(), "Invalid corner index."); return m_vGloPos[i];}

			private:
				// This scvf separates the scv with the ids given in "from" and "to"
				// The computed normal points in direction from->to
				size_t m_from, m_to;

				// ordering is:
				// 2D: edgeMidPoint, CenterOfElement
				// 3D: edgeMidPoint, Side one, CenterOfElement, Side two
				MathVector<dim> m_vLocPos[m_numCorners]; // local corners of scvf
				MathVector<world_dim> m_vGloPos[m_numCorners]; // global corners of scvf
				MidID m_midId[m_numCorners]; // dimension and id of object, that's midpoint bounds the scvf

				// scvf part
				MathVector<dim> localIP; // local integration point
				MathVector<world_dim> globalIP; // global intergration point
				MathVector<world_dim> Normal; // normal (incl. area)

				// shapes and derivatives
				std::vector<number> vShape; // shapes at ip
				std::vector<MathVector<dim> > localGrad; // local grad at ip
				std::vector<MathVector<world_dim> > globalGrad; // global grad at ip
				MathMatrix<world_dim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

		class SCV
		{
			private:
				// let outer class access private members
				friend class HFV1Geometry<TElem, TWorldDim>;

				// number of integration points
				static const size_t m_numIP = 1;

				// Number of corners of scvf
				static const size_t m_maxNumCorners = hanging_finite_volume_traits<ref_elem_type, TWorldDim>::MaxNumCornersOfSCV;

				// type of element the subcontrol volume represents
				typedef typename hanging_finite_volume_traits<ref_elem_type, TWorldDim>::scv_type scv_type;

			public:
				SCV() : m_numCorners(m_maxNumCorners) {};

				/// node id that this scv is associated to
				inline size_t node_id() const {return nodeId;}

				/// number of integration points
				inline size_t num_ip() const {return m_numIP;}

				/// local integration point of scv
				inline const MathVector<dim>& local_ip(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return m_vLocPos[0];}

				/// global integration point
				inline const MathVector<world_dim>& global_ip(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return m_vGloPos[0];}

				/// volume of scv
				inline number volume() const {return vol;}

				/// number of corners, that bound the scvf
				inline size_t num_corners() const {return m_numCorners;}

				/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t i) const
					{UG_ASSERT(i < num_corners(), "Invalid corner index."); return m_vLocPos[i];}

				/// return glbal corner number i
				inline const MathVector<world_dim>& global_corner(size_t i) const
					{UG_ASSERT(i < num_corners(), "Invalid corner index."); return m_vGloPos[i];}

			private:
				size_t nodeId; // node id of associated node
				size_t m_numCorners;

				// The ordering is: Corner, ...
				MathVector<dim> m_vLocPos[m_maxNumCorners]; // local position of node
				MathVector<world_dim> m_vGloPos[m_maxNumCorners]; // global position of node
				MidID m_midId[m_maxNumCorners]; // dimension and id of object, that's midpoint bounds the scv
				number vol;
		};


	public:
		HFV1Geometry() : m_pElem(NULL)
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

		bool update(TElem* elem, const ISubsetHandler& ish, const MathVector<world_dim>* vCornerCoords)
		{
			// If already update for this element, do nothing
			if(m_pElem == elem) return true;
			else m_pElem = elem;

			// get grid
			Grid& grid = *ish.get_assigned_grid();

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
			std::vector<EdgeBase*> vEdges;
			CollectEdgesSorted(vEdges, grid, elem);

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
				switch(vEdges[i]->shared_pipe_section())
				{
				case SPSEDGE_CONSTRAINED_EDGE:
				case SPSEDGE_EDGE:
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

				case SPSEDGE_CONSTRAINING_EDGE:
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

				default: UG_LOG("Cannot detect type of edge.\n"); return false;
				}
			}

			// for 3d case also check faces for hanging nodes
			if(dim == 3)
			{
				std::vector<Face*> vFaces;
				CollectFacesSorted(vFaces, grid, elem);

				// compute Nodes
				MathVector<dim> locSideMid;
				MathVector<world_dim> gloSideMid;
				for(size_t i = 0; i < vFaces.size(); ++i)
				{
					///////////
					// case QUADRILATERAL with all edges hanging and hanging node in middle
					///////////
					if(vFaces[i]->shared_pipe_section() == SPSFACE_CONSTRAINING_QUADRILATERAL)
					{
						// insert hanging node in list of nodes
						const size_t newNodeId = m_gloMid[0].size();
						m_gloMid[0].resize(newNodeId + 1);
						m_locMid[0].resize(newNodeId + 1);

						// compute position of new (hanging) node
						compute_side_midpoints(i, m_locMid[0][newNodeId], m_gloMid[0][newNodeId]);

						// loop constrained faces
						for(size_t j = 0; j < m_rRefElem.num_obj_of_obj(2, i, 1); ++j)
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
					else if (vFaces[i]->shared_pipe_section() == SPSFACE_CONSTRAINING_TRIANGLE)
					{
						// compute position of new (hanging) node
						compute_side_midpoints(i, locSideMid, gloSideMid);

						// loop constrained faces
						for(size_t j = 0; j < m_rRefElem.num_obj_of_obj(2, i, 1); ++j)
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
							MathVector<world_dim> gloSmallSideMid;

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
						for(size_t j = 0; j < m_rRefElem.num_obj_of_obj(2, i, 1); ++j)
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
				HangingNormalOnSCVF<ref_elem_type, world_dim>(m_vSCVF[i].Normal, &(m_vSCVF[i].m_vGloPos[0]));

				// TODO: In 3D check orientation
				if(dim == 3)
				{
					const size_t from = m_vSCVF[i].m_from;
					const size_t to = m_vSCVF[i].m_to;

					MathVector<world_dim> diff;
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

                   m_vSCV[i].vol = ElementSize<typename SCV::scv_type, world_dim>(&(m_vSCV[i].m_vGloPos[0]));;
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

					m_vSCV[2*i + 0].vol = ElementSize<typename SCV::scv_type, world_dim>(&(m_vSCV[2*i + 0].m_vGloPos[0]));;
					m_vSCV[2*i + 1].vol = ElementSize<typename SCV::scv_type, world_dim>(&(m_vSCV[2*i + 1].m_vGloPos[0]));;
				}
			}

			/////////////////////////
			// Shapes and Derivatives
			/////////////////////////
			m_rMapping.update(vCornerCoords);

			for(size_t i = 0; i < num_scvf(); ++i)
			{
				if(!m_rMapping.jacobian_transposed_inverse(m_vSCVF[i].localIP, m_vSCVF[i].JtInv))
					{UG_LOG("Cannot compute jacobian transposed.\n"); return false;}
				if(!m_rMapping.jacobian_det(m_vSCVF[i].localIP, m_vSCVF[i].detj))
					{UG_LOG("Cannot compute jacobian determinate.\n"); return false;}

				const LocalShapeFunctionSet<ref_elem_type>& TrialSpace =
						LocalShapeFunctionSetProvider::
							get_local_shape_function_set<ref_elem_type>
								(LocalShapeFunctionSetID(LocalShapeFunctionSetID::LAGRANGE, 1));

				const size_t num_sh = ref_elem_type::num_corners;
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

			// 	loop ips
				for(size_t ip = 0; ip < rSCV.num_ip(); ++ip)
				{
					m_vGlobSCVIP.push_back(rSCV.global_ip(ip));
					m_vLocSCVIP.push_back(rSCV.local_ip(ip));
				}
			}

		// 	loop Sub Control Volumes Faces (SCVF)
			m_vGlobSCVFIP.clear();
			m_vLocSCVFIP.clear();
			for(size_t i = 0; i < num_scvf(); ++i)
			{
			//	get current SCVF
				const SCVF& rSCVF = scvf(i);

			// 	loop ips
				for(size_t ip = 0; ip < rSCVF.num_ip(); ++ip)
				{
					m_vGlobSCVFIP.push_back(rSCVF.global_ip(ip));
					m_vLocSCVFIP.push_back(rSCVF.local_ip(ip));
				}
			}


			//print();
			return true;
		}

		// debug output
		void print()
		{
			UG_LOG("\nFVG hanging debug output\n");
			UG_LOG("LocalCenter=" << m_locMid << ", globalCenter="<<m_gloMid<<"\n");
			for(size_t i = 0; i < m_vSCV.size(); ++i)
			{
				UG_LOG(i<<" SCV: ");
				UG_LOG("node_id=" << m_vSCV[i].node_id());
				UG_LOG(", local_pos="<< m_vSCV[i].local_ip(0));
				UG_LOG(", global_pos="<< m_vSCV[i].global_ip(0));
				UG_LOG(", vol=" << m_vSCV[i].volume());
				UG_LOG("\n    localCorner=" << m_vSCV[i].m_vLocPos[0]);
				UG_LOG(", localSide1=" << m_vSCV[i].m_vLocPos[1]);
				UG_LOG(", localCenter=" << m_vSCV[i].m_vLocPos[2]);
				UG_LOG(", localSide2=" << m_vSCV[i].m_vLocPos[3]);
				UG_LOG("\n    globalCorner=" << m_vSCV[i].m_vGloPos[0]);
				UG_LOG(", globalSide1=" << m_vSCV[i].m_vGloPos[1]);
				UG_LOG(", globalCenter=" << m_vSCV[i].m_vGloPos[2]);
				UG_LOG(", globalSide2=" << m_vSCV[i].m_vGloPos[3]);

				UG_LOG("\n");
			}
			UG_LOG("\n");
			for(size_t i = 0; i < m_vSCVF.size(); ++i)
			{
				UG_LOG(i<<" SCVF: ");
				UG_LOG("from=" << m_vSCVF[i].from()<<", to="<<m_vSCVF[i].to());
				UG_LOG(", local_pos="<< m_vSCVF[i].local_ip(0));
				UG_LOG(", global_pos="<< m_vSCVF[i].global_ip(0));
				UG_LOG(", normal=" << m_vSCVF[i].normal());
				UG_LOG("\n    Shapes:\n");
				for(size_t sh=0; sh < m_vSCVF[i].num_sh(); ++sh)
				{
					UG_LOG("         " <<sh << ": shape="<<m_vSCVF[i].shape(sh, 0));
					UG_LOG(", global_grad="<<m_vSCVF[i].global_grad(sh, 0));
					UG_LOG(", local_grad="<<m_vSCVF[i].local_grad(sh, 0));
					UG_LOG("\n");
				}
			}
			UG_LOG("\n");
		}

	public:
		/// number of SubControlVolumeFaces
		inline size_t num_scvf() const {return m_vSCVF.size();}

		/// const access to SubControlVolumeFace number i
		inline const SCVF& scvf(size_t i) const
			{UG_ASSERT(i < num_scvf(), "Invalid Index."); return m_vSCVF[i];}

		/// number of SubControlVolumes
		size_t num_scv() const {return m_vSCV.size();}

		/// const access to SubControlVolume number i
		inline const SCV& scv(size_t i) const
			{UG_ASSERT(i < num_scv(), "Invalid Index."); return m_vSCV[i];}

	public:
		/// returns all ips of scv as they appear in scv loop
		const MathVector<world_dim>* scvf_global_ips() const {return &m_vGlobSCVFIP[0];}

		/// returns number of all scv ips
		size_t num_scvf_global_ips() const {return m_vGlobSCVFIP.size();}

		/// returns all ips of scv as they appear in scv loop
		const MathVector<dim>* scvf_local_ips() const {return &m_vLocSCVFIP[0];}

		/// returns number of all scv ips
		size_t num_scvf_local_ips() const {return m_vLocSCVFIP.size();}

		/// returns all ips of scv as they appear in scv loop
		const MathVector<world_dim>* scv_global_ips() const {return &m_vGlobSCVIP[0];}

		/// returns number of all scv ips
		size_t num_scv_global_ips() const {return m_vGlobSCVIP.size();}

		/// returns all ips of scv as they appear in scv loop
		const MathVector<dim>* scv_local_ips() const {return &m_vLocSCVIP[0];}

		/// returns number of all scv ips
		size_t num_scv_local_ips() const {return m_vLocSCVIP.size();}

	protected:
		std::vector<MathVector<world_dim> > m_vGlobSCVFIP;
		std::vector<MathVector<dim> > m_vLocSCVFIP;
		std::vector<MathVector<world_dim> > m_vGlobSCVIP;
		std::vector<MathVector<dim> > m_vLocSCVIP;

	protected:
		void copy_local_corners(SCVF& scvf)
		{
			for(size_t i = 0; i < scvf.num_corners(); ++i)
			{
				const size_t dim = scvf.m_midId[i].dim;
				const size_t id = scvf.m_midId[i].id;
				scvf.m_vLocPos[i] = m_locMid[dim][id];
			}
		}

		void copy_global_corners(SCVF& scvf)
		{
			for(size_t i = 0; i < scvf.num_corners(); ++i)
			{
				const size_t dim = scvf.m_midId[i].dim;
				const size_t id = scvf.m_midId[i].id;
				scvf.m_vGloPos[i] = m_gloMid[dim][id];
			}
		}

		void copy_local_corners(SCV& scv)
		{
			for(size_t i = 0; i < scv.num_corners(); ++i)
			{
				const size_t dim_ = scv.m_midId[i].dim;
				const size_t id = scv.m_midId[i].id;
				UG_ASSERT(dim_ >= 0 && dim_ <= (size_t)dim, "Dimension wrong");
				UG_ASSERT(id < m_locMid[dim_].size(), "id " << id << " in dim="<<dim_<<" wrong. (size is "<< m_locMid[dim_].size()<<")\n");
				scv.m_vLocPos[i] = m_locMid[dim_][id];
			}
		}

		void copy_global_corners(SCV& scv)
		{
			for(size_t i = 0; i < scv.num_corners(); ++i)
			{
				const size_t dim_ = scv.m_midId[i].dim;
				const size_t id = scv.m_midId[i].id;
				UG_ASSERT(dim_ >= 0 && dim_ <= (size_t)dim, "Dimension wrong");
				UG_ASSERT(id < m_gloMid[dim_].size(), "id " << id << " in dim="<<dim_<<" wrong. (size is "<< m_gloMid[dim_].size()<<")\n");
				scv.m_vGloPos[i] = m_gloMid[dim_][id];
			}
		}

		// i = number of side
		void compute_side_midpoints(size_t i, MathVector<dim>& locSideMid, MathVector<world_dim>& gloSideMid)
		{
			const size_t coID0 = m_rRefElem.id(2, i, 0, 0);
			locSideMid = m_locMid[0][coID0];
			gloSideMid = m_gloMid[0][coID0];

			// add corner coordinates of the corners of the geometric object
			for(size_t j = 1; j < m_rRefElem.num_obj_of_obj(2, i, 0); ++j)
			{
				const size_t coID = m_rRefElem.id(2, i, 0, j);
				locSideMid += m_locMid[0][coID];
				gloSideMid += m_gloMid[0][coID];
			}

			// scale for correct averaging
			locSideMid *= 1./(m_rRefElem.num_obj_of_obj(2, i, 0));
			gloSideMid *= 1./(m_rRefElem.num_obj_of_obj(2, i, 0));
		}

		// i, j, k, l = number nodes
		void compute_side_midpoints(size_t i, size_t j, size_t k, size_t l,
									MathVector<dim>& locSideMid, MathVector<world_dim>& gloSideMid)
		{
			VecAdd(locSideMid, m_locMid[0][i], m_locMid[0][j], m_locMid[0][k], m_locMid[0][l]);
			VecAdd(gloSideMid, m_gloMid[0][i], m_gloMid[0][j], m_gloMid[0][k], m_gloMid[0][l]);

			// scale for correct averaging
			locSideMid *= 0.25;
			gloSideMid *= 0.25;
		}

		// i, j, k = number nodes
		void compute_side_midpoints(size_t i, size_t j, size_t k,
									MathVector<dim>& locSideMid, MathVector<world_dim>& gloSideMid)
		{
			VecAdd(locSideMid, m_locMid[0][i], m_locMid[0][j], m_locMid[0][k]);
			VecAdd(gloSideMid, m_gloMid[0][i], m_gloMid[0][j], m_gloMid[0][k]);

			// scale for correct averaging
			locSideMid *= 1./3.;
			gloSideMid *= 1./3.;
		}

		// returns edgeID of child edge, that has corner co. i is the natural edge index
		size_t get_child_edge_of_corner(size_t i, size_t co)
        {
         	for(size_t e = 0; e < m_vNatEdgeInfo[i].num_child_edges(); ++e)
         	{
         		const size_t childId = m_vNatEdgeInfo[i].child_edge(e);
         		if(m_vNewEdgeInfo[childId].from() == co)
         			return childId;
         		if(m_vNewEdgeInfo[childId].to() == co)
         			return childId;
         	}
         	return -1;
        }

	private:
		struct NewEdgeInfo
		{
			friend class HFV1Geometry<TElem, TWorldDim>;
		public:
			NewEdgeInfo() : m_from(-1), m_to(-1){}
            size_t from() const {return m_from;}
            size_t to() const {return m_to;}
		private:
			size_t m_from;
			size_t m_to;
		};

        struct NatEdgeInfo
        {
                friend class HFV1Geometry<TElem, TWorldDim>;
            public:
            	NatEdgeInfo() : nodeId(-1), numChildEdges(0) {for(size_t i=0; i<2; ++i) childEdge[i] = -1;}
                bool is_hanging() const {return numChildEdges == 2;}
                size_t node_id() const {UG_ASSERT(is_hanging(), "Should only be called, if edge is hanging."); return nodeId;}
                size_t num_child_edges() const {return numChildEdges;}
                size_t child_edge(size_t i) const {UG_ASSERT(i < num_child_edges(), "Wrong id."); return childEdge[i];}

            private:
                size_t nodeId;
                size_t numChildEdges;
                size_t childEdge[2];
        };

         // help array: maps NaturalEdge -> NodeOnEdge (-1 if no node on edge)
        std::vector<NatEdgeInfo> m_vNatEdgeInfo;
        std::vector<NewEdgeInfo> m_vNewEdgeInfo;

	private:
		// pointer to current element
		TElem* m_pElem;

		std::vector<MathVector<dim> > m_locMid[dim+1];
		std::vector<MathVector<world_dim> > m_gloMid[dim+1];

		// SubControlVolumeFaces
		std::vector<SCVF> m_vSCVF;

		// SubControlVolumes
		std::vector<SCV> m_vSCV;

		// Reference Mapping
		ReferenceMapping<ref_elem_type, world_dim> m_rMapping;

		// Reference Element
		ref_elem_type m_rRefElem;

};

}

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__HANGING_FINITE_VOLUME_GEOMETRY__ */
