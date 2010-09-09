/*
 * finite_volume_geometry.h
 *
 *  Created on: 04.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_GEOMETRY__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_GEOMETRY__

// extern libraries
#include <cmath>
#include <vector>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lib_grid.h"

// library intern includes
#include "lib_discretization/common/common.h"
#include "../../reference_element/reference_element.h"
#include "../../local_shape_function_set/local_shape_function_set_factory.h"
#include "./finite_volume_util.h"

namespace ug{


template <	typename TElem,
			int TWorldDim>
class FV1Geometry {
	private:
		// type of reference element
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

		// number of SubControlVolumes
		static const size_t m_numSCV = ref_elem_type::num_corners;

		// type of SubControlVolume
		typedef typename finite_volume_traits<ref_elem_type, TWorldDim>::scv_type scv_type;

		// number of SubControlVolumeFaces
		static const size_t m_numSCVF = ref_elem_type::num_edges;

	public:
		// dimension of reference element
		static const int dim = ref_elem_type::dim;

		// dimension of world
		static const int world_dim = TWorldDim;

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
			private:
				// let outer class access private members
				friend class FV1Geometry<TElem, TWorldDim>;

				// number of integration points
				static const size_t m_numIP = 1;

				// Number of corners of scvf
				static const size_t m_numCorners = finite_volume_traits<ref_elem_type, TWorldDim>::NumCornersOfSCVF;

			public:
				SCVF() {}

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
				// 1D: edgeMidPoint
				// 2D: edgeMidPoint, CenterOfElement
				// 3D: edgeMidPoint, Side one, CenterOfElement, Side two
				MathVector<dim> m_vLocPos[m_numCorners]; // local corners of scvf
				MathVector<world_dim> m_vGloPos[m_numCorners]; // global corners of scvf
				MidID midId[m_numCorners]; // dimension and id of object, that's midpoint bounds the scvf

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
				friend class FV1Geometry<TElem, TWorldDim>;

				// number of integration points
				static const size_t m_numIP = 1;

				// Number of corners of scvf
				static const size_t m_maxNumCorners = finite_volume_traits<ref_elem_type, TWorldDim>::MaxNumCornersOfSCV;

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
				MathVector<dim> m_vLocPos[m_maxNumCorners]; // local position of node
				MathVector<world_dim> m_vGloPos[m_maxNumCorners]; // global position of node
				MidID midId[m_maxNumCorners]; // dimension and id of object, that's midpoint bounds the scv
				number vol;
		};

	public:
		FV1Geometry()
		{
			// set corners of element as local centers of nodes
			for(size_t i = 0; i < m_refElem.num_obj(0); ++i)
				m_locMid[0][i] = m_refElem.corner(i);

			// compute local midpoints for all geometric objects with  0 < d <= dim
			for(int d = 1; d <= dim; ++d)
			{
				// loop geometric objects of dimension d
				for(size_t i = 0; i < m_refElem.num_obj(d); ++i)
				{
					// set first node
					const size_t coID0 = m_refElem.id(d, i, 0, 0);
					m_locMid[d][i] = m_locMid[0][coID0];

					// add corner coordinates of the corners of the geometric object
					for(size_t j = 1; j < m_refElem.num_obj_of_obj(d, i, 0); ++j)
					{
						const size_t coID = m_refElem.id(d, i, 0, j);
						m_locMid[d][i] += m_locMid[0][coID];
					}

					// scale for correct averaging
					m_locMid[d][i] *= 1./(m_refElem.num_obj_of_obj(d, i, 0));
				}
			}

			// set up local informations for SubControlVolumeFaces (scvf)
			// each scvf is associated to one edge of the element
			for(size_t i = 0; i < num_scvf(); ++i)
			{
				m_vSCVF[i].m_from = m_refElem.id(1, i, 0, 0);
				m_vSCVF[i].m_to = m_refElem.id(1, i, 0, 1);

				// set mid ids
				{
					// start at edge midpoint
					m_vSCVF[i].midId[0] = MidID(1,i);

					// loop up dimension
					if(dim == 2)
					{
						m_vSCVF[i].midId[1] = MidID(dim, 0); // center of element
					}
					else if (dim == 3)
					{
						m_vSCVF[i].midId[1] = MidID(2, m_refElem.id(1, i, 2, 0)); // side 0
						m_vSCVF[i].midId[2] = MidID(dim, 0); // center of element
						m_vSCVF[i].midId[3] = MidID(2, m_refElem.id(1, i, 2, 1)); // side 1
					}
				}

				// copy local corners of scvf
				copy_local_corners(m_vSCVF[i]);

				// integration point
				average_corners<dim>(m_vSCVF[i].localIP, m_vSCVF[i].m_vLocPos, SCVF::m_numCorners);
/*				VecInterpolateLinear(	m_vSCVF[i].localIP,
										m_vSCVF[i].m_vLocPos[0],
										m_vSCVF[i].m_vLocPos[1],
										0.5);
*/			}

			// set up local informations for SubControlVolumes (scv)
			// each scv is associated to one corner of the element
			for(size_t i = 0; i < num_scv(); ++i)
			{
				m_vSCV[i].nodeId = i;

				if(dim == 1)
				{
					m_vSCV[i].midId[0] = MidID(0, i); // set node as corner of scv
					m_vSCV[i].midId[1] = MidID(dim, 0);	// center of element
				}
				else if(dim == 2)
				{
					m_vSCV[i].midId[0] = MidID(0, i); // set node as corner of scv
					m_vSCV[i].midId[1] = MidID(1, m_refElem.id(0, i, 1, 0)); // edge 1
					m_vSCV[i].midId[2] = MidID(dim, 0);	// center of element
					m_vSCV[i].midId[3] = MidID(1, m_refElem.id(0, i, 1, 1)); // edge 2
				}
				else if(dim == 3 && (ref_elem_type::REFERENCE_OBJECT_ID != ROID_PYRAMID || i != num_scv()-1))
				{
					m_vSCV[i].midId[0] = MidID(0, i); // set node as corner of scv
					m_vSCV[i].midId[1] = MidID(1, m_refElem.id(0, i, 1, 1)); // edge 1
					m_vSCV[i].midId[2] = MidID(2, m_refElem.id(0, i, 2, 0)); // face 0
					m_vSCV[i].midId[3] = MidID(1, m_refElem.id(0, i, 1, 0)); // edge 0
					m_vSCV[i].midId[4] = MidID(1, m_refElem.id(0, i, 1, 2)); // edge 2
					m_vSCV[i].midId[5] = MidID(2, m_refElem.id(0, i, 2, 2)); // face 2
					m_vSCV[i].midId[6] = MidID(dim, 0);	// center of element
					m_vSCV[i].midId[7] = MidID(2, m_refElem.id(0, i, 2, 1)); // face 1
				}
				// TODO: Implement last ControlVolume for Pyramid
				else if(dim == 3 && ref_elem_type::REFERENCE_OBJECT_ID == ROID_PYRAMID && i == num_scv()-1)
				{
					// this scv has 10 corners
					m_vSCV[i].m_numCorners = 10;
					//UG_ASSERT(0, "Last SCV for Pyramid must be implemented");
				}
				else {UG_ASSERT(0, "Dimension higher that 3 not implemented.");}

				// copy local corners of scv
				copy_local_corners(m_vSCV[i]);
			}

			/////////////////////////
			// Shapes and Derivatives
			/////////////////////////
			for(size_t i = 0; i < num_scvf(); ++i)
			{
				const LocalShapeFunctionSet<ref_elem_type>& TrialSpace =
						LocalShapeFunctionSetFactory::inst().get_local_shape_function_set<ref_elem_type>(LSFS_LAGRANGEP1);

				const size_t num_sh = m_numSCV;
				m_vSCVF[i].vShape.resize(num_sh);
				m_vSCVF[i].localGrad.resize(num_sh);
				m_vSCVF[i].globalGrad.resize(num_sh);
				for(size_t sh = 0 ; sh < num_sh; ++sh)
				{
					if(!TrialSpace.evaluate(sh, m_vSCVF[i].localIP, (m_vSCVF[i].vShape)[sh]))
						{UG_LOG("Cannot evaluate local shape.\n"); UG_ASSERT(0, "Error in Constructor.");}
					if(!TrialSpace.evaluate_grad(sh, m_vSCVF[i].localIP, (m_vSCVF[i].localGrad)[sh]))
						{UG_LOG("Cannot evaluate local grad.\n"); UG_ASSERT(0, "Error in Constructor.");}
				}
			}
		}

		bool update(TElem* elem, const Grid& mg, MathVector<world_dim> vCornerCoords[])
		{
			// remember global position of nodes
			for(size_t i = 0; i < m_refElem.num_obj(0); ++i)
				m_gloMid[0][i] = vCornerCoords[i];

			// compute global midpoints for all geometric objects with  0 < d <= dim
			for(int d = 1; d <= dim; ++d)
			{
				// loop geometric objects of dimension d
				for(size_t i = 0; i < m_refElem.num_obj(d); ++i)
				{
					// set first node
					const size_t coID0 = m_refElem.id(d, i, 0, 0);
					m_gloMid[d][i] = m_gloMid[0][coID0];

					// add corner coordinates of the corners of the geometric object
					for(size_t j = 1; j < m_refElem.num_obj_of_obj(d, i, 0); ++j)
					{
						const size_t coID = m_refElem.id(d, i, 0, j);
						m_gloMid[d][i] += m_gloMid[0][coID];
					}

					// scale for correct averaging
					m_gloMid[d][i] *= 1./(m_refElem.num_obj_of_obj(d, i, 0));
				}
			}

			// compute global informations for scvf
			for(size_t i = 0; i < num_scvf(); ++i)
			{
				// copy local corners of scvf
				copy_global_corners(m_vSCVF[i]);

				// integration point
				average_corners<world_dim>(m_vSCVF[i].globalIP, m_vSCVF[i].m_vGloPos, SCVF::m_numCorners);
/*				VecInterpolateLinear(	m_vSCVF[i].globalIP,
										m_vSCVF[i].m_vGloPos[0],
										m_vSCVF[i].m_vGloPos[1],
										0.5);
*/
				// normal on scvf
				NormalOnSCVF<ref_elem_type, world_dim>(m_vSCVF[i].Normal, m_vSCVF[i].m_vGloPos);
			}

			// compute size of scv
			for(size_t i = 0; i < num_scv(); ++i)
			{
				// copy global corners
				copy_global_corners(m_vSCV[i]);

				// compute volume of scv
				if(m_vSCV[i].m_numCorners != 10)
				{
					m_vSCV[i].vol = ElementSize<scv_type, world_dim>(m_vSCV[i].m_vGloPos);
				}
				else
				{
					// special case for pyramid, last scv
				}
			}

			// compute shapes and derivatives
			m_mapping.update(vCornerCoords);

			for(size_t i = 0; i < num_scvf(); ++i)
			{
				if(!m_mapping.jacobian_transposed_inverse(m_vSCVF[i].localIP, m_vSCVF[i].JtInv))
					{UG_LOG("Cannot compute jacobian transposed.\n"); return false;}
				if(!m_mapping.jacobian_det(m_vSCVF[i].localIP, m_vSCVF[i].detj))
					{UG_LOG("Cannot compute jacobian determinate.\n"); return false;}

				for(size_t sh = 0 ; sh < num_scv(); ++sh)
					MatVecMult((m_vSCVF[i].globalGrad)[sh], m_vSCVF[i].JtInv, (m_vSCVF[i].localGrad)[sh]);
			}

			//print();
			return true;
		}

	protected:
		template <int TDim>
		void average_corners(MathVector<TDim>& vOut, const MathVector<TDim>* vCornerCoords, size_t num)
		{
			vOut = vCornerCoords[0];
			for(size_t j = 1; j < num; ++j)
			{
				vOut += vCornerCoords[j];
			}
			vOut *= 1./(number)num;
		}

		void copy_local_corners(SCVF& scvf)
		{
			for(size_t i = 0; i < scvf.m_numCorners; ++i)
			{
				const size_t dim = scvf.midId[i].dim;
				const size_t id = scvf.midId[i].id;
				scvf.m_vLocPos[i] = m_locMid[dim][id];
			}
		}

		void copy_global_corners(SCVF& scvf)
		{
			for(size_t i = 0; i < scvf.m_numCorners; ++i)
			{
				const size_t dim = scvf.midId[i].dim;
				const size_t id = scvf.midId[i].id;
				scvf.m_vGloPos[i] = m_gloMid[dim][id];
			}
		}

		void copy_local_corners(SCV& scv)
		{
			for(size_t i = 0; i < scv.m_numCorners; ++i)
			{
				const size_t dim = scv.midId[i].dim;
				const size_t id = scv.midId[i].id;
				scv.m_vLocPos[i] = m_locMid[dim][id];
			}
		}

		void copy_global_corners(SCV& scv)
		{
			for(size_t i = 0; i < scv.m_numCorners; ++i)
			{
				const size_t dim = scv.midId[i].dim;
				const size_t id = scv.midId[i].id;
				scv.m_vGloPos[i] = m_gloMid[dim][id];
			}
		}

	protected:
		// debug output
		void print() const
		{
			UG_LOG("\nFV debug output\n");
			for(size_t i = 0; i < num_scv(); ++i)
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
			for(size_t i = 0; i < num_scvf(); ++i)
			{
				UG_LOG(i<<" SCVF: ");
				UG_LOG("from=" << m_vSCVF[i].from()<<", to="<<m_vSCVF[i].to());
				UG_LOG(", local_pos="<< m_vSCVF[i].local_ip(0));
				UG_LOG(", global_pos="<< m_vSCVF[i].global_ip(0));
				UG_LOG(", normal=" << m_vSCVF[i].normal());
				UG_LOG("\n    localEdgeMid=" << m_vSCVF[i].m_vLocPos[0]);
				UG_LOG(", localCenter=" << m_vSCVF[i].m_vLocPos[1]);
				UG_LOG(", globalEdgeMid=" << m_vSCVF[i].m_vGloPos[0]);
				UG_LOG(", globalCenter=" << m_vSCVF[i].m_vGloPos[1]);
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
		inline size_t num_scvf() const {return m_numSCVF;};

		/// const access to SubControlVolumeFace number i
		inline const SCVF& scvf(size_t i) const
		{
			UG_ASSERT(i < num_scvf(), "Invalid Index.");
			return m_vSCVF[i];
		}

		/// number of SubControlVolumes
		size_t num_scv() const {return m_numSCV;}

		/// const access to SubControlVolume number i
		inline const SCV& scv(size_t i) const
		{
			UG_ASSERT(i < num_scv(), "Invalid Index.");
			return m_vSCV[i];
		}


	private:
		// local and global geom object midpoints for each dimension (most objects in 1 dim, i.e. number of edges)
		MathVector<dim> m_locMid[dim+1][m_numSCVF];
		MathVector<world_dim> m_gloMid[dim+1][m_numSCVF];

		SCVF m_vSCVF[m_numSCVF];
		SCV m_vSCV[m_numSCV];

		ReferenceMapping<ref_elem_type, world_dim> m_mapping;
		ref_elem_type m_refElem;

};

}

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_GEOMETRY__ */
