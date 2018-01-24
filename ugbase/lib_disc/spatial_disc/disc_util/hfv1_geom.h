/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__HANGING_FINITE_VOLUME_GEOMETRY__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__HANGING_FINITE_VOLUME_GEOMETRY__

// extern libraries
#include <cmath>
#include <vector>

// other ug4 modules
#include "common/common.h"

// library intern includes
#include "lib_disc/reference_element/reference_element.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "fv_util.h"
#include "fv1_geom.h"

namespace ug{

template <	typename TElem, int TWorldDim>
class HFV1Geometry : public FVGeometryBase{
	private:
	/// type of reference element
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	/// number of SubControlVolumes
		static const size_t m_numNaturalSCV = ref_elem_type::numCorners;

	/// number of SubControlVolumeFaces
		static const size_t m_numNaturalSCVF = ref_elem_type::numEdges;

	public:
	/// dimension of reference element
		static const int dim = ref_elem_type::dim;

	/// dimension of world
		static const int worldDim = TWorldDim;

	/// Hanging node flag: this Geometry does support hanging nodes
		static const bool usesHangingNodes = true;

	/// flag indicating if local data may change
		static const bool staticLocalData = true;

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
			/// dimension of reference element
				static const int dim = ref_elem_type::dim;

			/// dimension of world
				static const int worldDim = TWorldDim;

			private:
			/// let outer class access private members
				friend class HFV1Geometry<TElem, TWorldDim>;

			/// number of integration points
				static const size_t m_numIP = 1;

			/// Number of corners of scvf
				static const size_t m_numCorners = hfv1_traits<ref_elem_type, TWorldDim>::NumCornersOfSCVF;

			public:
				SCVF() {};

			/// index of SubControlVolume on one side of the scvf
			/// NO! return value is the associated NODE_ID; this is not the same!
				inline size_t from() const {return m_from;}

			/// index of SubControlVolume on one side of the scvf
			/// NO! return value is the associated NODE_ID; this is not the same!
				inline size_t to() const {return m_to;}

			/// number of integration points on scvf
				inline size_t num_ip() const {return m_numIP;}

			/// local integration point of scvf
				inline const MathVector<dim>& local_ip() const {return localIP;}

			/// global integration point of scvf
				inline const MathVector<worldDim>& global_ip() const {return globalIP;}

			/// normal on scvf (points direction "from"->"to"). Norm is equal to area
				inline const MathVector<worldDim>& normal() const {return Normal;} // includes area

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim,dim>& JTInv() const {return JtInv;}

			/// Determinante of Jacobian in integration point
				inline number detJ() const {return detj;}

			/// number of shape functions
				inline size_t num_sh() const {return vShape.size();}

			/// value of shape function i in integration point
				inline number shape(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid ip index"); return vShape[sh];}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid ip index"); return localGrad[sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector() const {return &localGrad[0];}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid ip index"); return globalGrad[sh];}

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector() const {return &globalGrad[0];}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return m_numCorners;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid corner index."); return m_vLocPos[co];}

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid corner index."); return m_vGloPos[co];}

			private:
				// This scvf separates the scv with the ids given in "from" and "to"
				// The computed normal points in direction from->to
				size_t m_from, m_to;

				// ordering is:
				// 2D: edgeMidPoint, CenterOfElement
				// 3D: edgeMidPoint, Side one, CenterOfElement, Side two
				MathVector<dim> m_vLocPos[m_numCorners]; // local corners of scvf
				MathVector<worldDim> m_vGloPos[m_numCorners]; // global corners of scvf
				MidID m_midId[m_numCorners]; // dimension and id of object, that's midpoint bounds the scvf

				// scvf part
				MathVector<dim> localIP; // local integration point
				MathVector<worldDim> globalIP; // global intergration point
				MathVector<worldDim> Normal; // normal (incl. area)

				// shapes and derivatives
				std::vector<number> vShape; // shapes at ip
				std::vector<MathVector<dim> > localGrad; // local grad at ip
				std::vector<MathVector<worldDim> > globalGrad; // global grad at ip
				MathMatrix<worldDim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

		class SCV
		{
			private:
			//  let outer class access private members
				friend class HFV1Geometry<TElem, TWorldDim>;

			/// number of integration points
				static const size_t m_numIP = 1;

			/// Number of corners of scvf
				static const size_t m_maxNumCorners = hfv1_traits<ref_elem_type, TWorldDim>::MaxNumCornersOfSCV;

			/// type of element the subcontrol volume represents
				typedef typename hfv1_traits<ref_elem_type, TWorldDim>::scv_type scv_type;

			public:
				SCV() : m_numCorners(m_maxNumCorners) {};

			/// node id that this scv is associated to
				inline size_t node_id() const {return nodeId;}

			/// number of integration points
				inline size_t num_ip() const {return m_numIP;}

			/// local integration point of scv
				inline const MathVector<dim>& local_ip() const {return m_vLocPos[0];}

			/// global integration point
				inline const MathVector<worldDim>& global_ip() const {return m_vGloPos[0];}

			/// volume of scv
				inline number volume() const {return vol;}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return m_numCorners;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid corner index."); return m_vLocPos[co];}

			/// return global corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid corner index."); return m_vGloPos[co];}

			private:
				size_t nodeId; // node id of associated node
				size_t m_numCorners;

				// The ordering is: Corner, ...
				MathVector<dim> m_vLocPos[m_maxNumCorners]; // local position of node
				MathVector<worldDim> m_vGloPos[m_maxNumCorners]; // global position of node
				MidID m_midId[m_maxNumCorners]; // dimension and id of object, that's midpoint bounds the scv
				number vol;
		};


	public:
	///	constructor
		HFV1Geometry();

	///	update values for an element
		void update(GridObject* pElem, const MathVector<worldDim>* vCornerCoords,
		            			 const ISubsetHandler* ish = NULL);

	// 	debug output
		void print();

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

	/// number of shape functions
		inline size_t num_sh() const
		{
			return m_numSh;
		};

	public:
	/// returns number of all scv ips
		size_t num_scvf_ips() const {return m_vGlobSCVFIP.size();}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<worldDim>* scvf_global_ips() const {return &m_vGlobSCVFIP[0];}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<dim>* scvf_local_ips() const {return &m_vLocSCVFIP[0];}

	/// returns number of all scv ips
		size_t num_scv_ips() const {return m_vGlobSCVIP.size();}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<worldDim>* scv_global_ips() const {return &m_vGlobSCVIP[0];}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<dim>* scv_local_ips() const {return &m_vLocSCVIP[0];}

	/// return local coords for node ID
		const MathVector<dim>& local_node_position(size_t nodeID) const
		{
			UG_ASSERT(nodeID < m_locMid[0].size(), "Invalid node id.");
			return m_locMid[0][nodeID];
		}

	/// return global coords for node ID
		const MathVector<worldDim>& global_node_position(size_t nodeID) const
		{
			UG_ASSERT(nodeID < m_gloMid[0].size(), "Invalid node id.");
			return m_gloMid[0][nodeID];
		}

	protected:
		std::vector<MathVector<worldDim> > m_vGlobSCVFIP;
		std::vector<MathVector<dim> > m_vLocSCVFIP;
		std::vector<MathVector<worldDim> > m_vGlobSCVIP;
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
				UG_ASSERT(dim_ <= (size_t)dim, "Dimension wrong");
				UG_ASSERT(id < m_locMid[dim_].size(), "id " << id << " in dim="
				          <<dim_<<" wrong. (size is "<< m_locMid[dim_].size()<<")\n");
				scv.m_vLocPos[i] = m_locMid[dim_][id];
			}
		}

		void copy_global_corners(SCV& scv)
		{
			for(size_t i = 0; i < scv.num_corners(); ++i)
			{
				const size_t dim_ = scv.m_midId[i].dim;
				const size_t id = scv.m_midId[i].id;
				UG_ASSERT(dim_ <= (size_t)dim, "Dimension wrong");
				UG_ASSERT(id < m_gloMid[dim_].size(), "id " << id << " in dim="
				          <<dim_<<" wrong. (size is "<< m_gloMid[dim_].size()<<")\n");
				scv.m_vGloPos[i] = m_gloMid[dim_][id];
			}
		}

		// i = number of side
		void compute_side_midpoints(size_t i,
		                            MathVector<dim>& locSideMid,
		                            MathVector<worldDim>& gloSideMid)
		{
			if(dim>1)
			{
				const size_t coID0 = m_rRefElem.id(2, i, 0, 0);
				locSideMid = m_locMid[0][coID0];
				gloSideMid = m_gloMid[0][coID0];

				// add corner coordinates of the corners of the geometric object
				for(size_t j = 1; j < m_rRefElem.num(2, i, 0); ++j)
				{
					const size_t coID = m_rRefElem.id(2, i, 0, j);
					locSideMid += m_locMid[0][coID];
					gloSideMid += m_gloMid[0][coID];
				}

				// scale for correct averaging
				locSideMid *= 1./(m_rRefElem.num(2, i, 0));
				gloSideMid *= 1./(m_rRefElem.num(2, i, 0));
			}
		}

		// i, j, k, l = number nodes
		void compute_side_midpoints(size_t i, size_t j, size_t k, size_t l,
									MathVector<dim>& locSideMid,
									MathVector<worldDim>& gloSideMid)
		{
			VecAdd(locSideMid, m_locMid[0][i], m_locMid[0][j], m_locMid[0][k], m_locMid[0][l]);
			VecAdd(gloSideMid, m_gloMid[0][i], m_gloMid[0][j], m_gloMid[0][k], m_gloMid[0][l]);

			// scale for correct averaging
			locSideMid *= 0.25;
			gloSideMid *= 0.25;
		}

		// i, j, k = number nodes
		void compute_side_midpoints(size_t i, size_t j, size_t k,
									MathVector<dim>& locSideMid, MathVector<worldDim>& gloSideMid)
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

    //  help array: maps NaturalEdge -> NodeOnEdge (-1 if no node on edge)
        std::vector<NatEdgeInfo> m_vNatEdgeInfo;
        std::vector<NewEdgeInfo> m_vNewEdgeInfo;

	private:
	// 	pointer to current element
		TElem* m_pElem;

		std::vector<MathVector<dim> > m_locMid[dim+1];
		std::vector<MathVector<worldDim> > m_gloMid[dim+1];

	// 	SubControlVolumeFaces
		std::vector<SCVF> m_vSCVF;

	// 	SubControlVolumes
		std::vector<SCV> m_vSCV;

	//	number of shape functions
		size_t m_numSh;

	// 	Reference Mapping
		ReferenceMapping<ref_elem_type, worldDim> m_rMapping;

	// 	Reference Element
		const ref_elem_type& m_rRefElem;
};

template <	int TDim, int TWorldDim = TDim>
class DimHFV1Geometry : public FVGeometryBase{
	private:
	/// number of SubControlVolumes
		size_t m_numNaturalSCV;

	/// number of SubControlVolumeFaces
		size_t m_numNaturalSCVF;
		
	/// max number of shapes
		static const size_t m_maxNSH = fv1_dim_traits<TDim, TWorldDim>::maxNSH;

	public:
	/// dimension of reference element
		static const int dim = TDim;

	/// dimension of world
		static const int worldDim = TWorldDim;

	/// Hanging node flag: this Geometry does support hanging nodes
		static const bool usesHangingNodes = true;

	/// flag indicating if local data may change
		static const bool staticLocalData = true;
		
	/// traits	
		typedef hdimfv1_traits<TDim> traits;
		
	/// element type names
		typedef typename traits::elem_type_0 elem_type_0;
		typedef typename traits::elem_type_1 elem_type_1;
		typedef typename traits::elem_type_2 elem_type_2;
		typedef typename traits::elem_type_3 elem_type_3;
		typedef typename traits::elem_type_4 elem_type_4;

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
			/// dimension of reference element
				static const int dim = TDim;

			/// dimension of world
				static const int worldDim = TWorldDim;

			private:
			/// let outer class access private members
				friend class DimHFV1Geometry<TDim, TWorldDim>;

			/// number of integration points
				static const size_t m_numIP = 1;

			/// Number of corners of scvf
				static const size_t m_numCorners = traits::NumCornersOfSCVF;

			public:
				SCVF() {};

			/// index of SubControlVolume on one side of the scvf
				inline size_t from() const {return m_from;}

			/// index of SubControlVolume on one side of the scvf
				inline size_t to() const {return m_to;}

			/// number of integration points on scvf
				inline size_t num_ip() const {return m_numIP;}

			/// local integration point of scvf
				inline const MathVector<dim>& local_ip() const {return localIP;}

			/// global integration point of scvf
				inline const MathVector<worldDim>& global_ip() const {return globalIP;}

			/// normal on scvf (points direction "from"->"to"). Norm is equal to area
				inline const MathVector<worldDim>& normal() const {return Normal;} // includes area

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim,dim>& JTInv() const {return JtInv;}

			/// Determinante of Jacobian in integration point
				inline number detJ() const {return detj;}

			/// number of shape functions
				inline size_t num_sh() const {return vShape.size();}

			/// value of shape function i in integration point
				inline number shape(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid ip index"); return vShape[sh];}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid ip index"); return localGrad[sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector() const {return &localGrad[0];}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid ip index"); return globalGrad[sh];}

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector() const {return &globalGrad[0];}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return m_numCorners;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid corner index."); return m_vLocPos[co];}

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid corner index."); return m_vGloPos[co];}

			private:
				// This scvf separates the scv with the ids given in "from" and "to"
				// The computed normal points in direction from->to
				size_t m_from, m_to;

				// ordering is:
				// 2D: edgeMidPoint, CenterOfElement
				// 3D: edgeMidPoint, Side one, CenterOfElement, Side two
				MathVector<dim> m_vLocPos[m_numCorners]; // local corners of scvf
				MathVector<worldDim> m_vGloPos[m_numCorners]; // global corners of scvf
				MidID m_midId[m_numCorners]; // dimension and id of object, that's midpoint bounds the scvf

				// scvf part
				MathVector<dim> localIP; // local integration point
				MathVector<worldDim> globalIP; // global intergration point
				MathVector<worldDim> Normal; // normal (incl. area)

				// shapes and derivatives
				std::vector<number> vShape; // shapes at ip
				std::vector<MathVector<dim> > localGrad; // local grad at ip
				std::vector<MathVector<worldDim> > globalGrad; // global grad at ip
				MathMatrix<worldDim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

		class SCV
		{
			private:
			//  let outer class access private members
				friend class DimHFV1Geometry<TDim, TWorldDim>;

			/// number of integration points
				static const size_t m_numIP = 1;

			/// Number of corners of scvf
				static const size_t m_maxNumCorners = traits::MaxNumCornersOfSCV;

			/// type of element the subcontrol volume represents
				typedef typename traits::scv_type scv_type;

			public:
				SCV() : m_numCorners(m_maxNumCorners) {};

			/// node id that this scv is associated to
				inline size_t node_id() const {return nodeId;}

			/// number of integration points
				inline size_t num_ip() const {return m_numIP;}

			/// local integration point of scv
				inline const MathVector<dim>& local_ip() const {return m_vLocPos[0];}

			/// global integration point
				inline const MathVector<worldDim>& global_ip() const {return m_vGloPos[0];}

			/// volume of scv
				inline number volume() const {return vol;}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return m_numCorners;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid corner index."); return m_vLocPos[co];}

			/// return global corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid corner index."); return m_vGloPos[co];}
					
			/// number of shape functions
				inline size_t num_sh() const {return numSH;}

			/// value of shape function i in integration point
				inline number shape(size_t sh) const {return vShape[sh];}

			/// vector of shape functions in ip point
				inline const number* shape_vector() const {return vShape;}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return localGrad[sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector() const {return localGrad;}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return globalGrad[sh];}

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector() const {return globalGrad;}
				
			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim,dim>& JTInv() const {return JtInv;}

			/// Determinant of Jacobian in integration point
				inline number detJ() const {return detj;}

			private:
				size_t nodeId; // node id of associated node
				size_t m_numCorners;

			// The ordering is: Corner, ...
				MathVector<dim> m_vLocPos[m_maxNumCorners]; // local position of node
				MathVector<worldDim> m_vGloPos[m_maxNumCorners]; // global position of node
				MidID m_midId[m_maxNumCorners]; // dimension and id of object, that's midpoint bounds the scv
				number vol;
				
			// shapes and derivatives
				size_t numSH;
				number vShape[m_maxNSH]; // shapes at ip
				MathVector<dim> localGrad[m_maxNSH]; // local grad at ip
				MathVector<worldDim> globalGrad[m_maxNSH]; // global grad at ip
				MathMatrix<worldDim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};


	public:
	///	constructor
		DimHFV1Geometry() : m_pElem(NULL), m_roid(ROID_UNKNOWN) {};
		
	/// update local data for an element
		void update_local_data();

	///	update values for an element
		void update(GridObject* pElem, const MathVector<worldDim>* vCornerCoords,
		            			 const ISubsetHandler* ish = NULL);

	// 	debug output
	//	void print();

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

	/// number of shape functions
		inline size_t num_sh() const
		{
			return m_numSh;
		};

	public:
	/// returns number of all scv ips
		size_t num_scvf_ips() const {return m_vGlobSCVFIP.size();}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<worldDim>* scvf_global_ips() const {return &m_vGlobSCVFIP[0];}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<dim>* scvf_local_ips() const {return &m_vLocSCVFIP[0];}

	/// returns number of all scv ips
		size_t num_scv_ips() const {return m_vGlobSCVIP.size();}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<worldDim>* scv_global_ips() const {return &m_vGlobSCVIP[0];}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<dim>* scv_local_ips() const {return &m_vLocSCVIP[0];}

	protected:
		std::vector<MathVector<worldDim> > m_vGlobSCVFIP;
		std::vector<MathVector<dim> > m_vLocSCVFIP;
		std::vector<MathVector<worldDim> > m_vGlobSCVIP;
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
				UG_ASSERT(dim_ <= (size_t)dim, "Dimension wrong");
				UG_ASSERT(id < m_locMid[dim_].size(), "id " << id << " in dim="
				          <<dim_<<" wrong. (size is "<< m_locMid[dim_].size()<<")\n");
				scv.m_vLocPos[i] = m_locMid[dim_][id];
			}
		}

		void copy_global_corners(SCV& scv)
		{
			for(size_t i = 0; i < scv.num_corners(); ++i)
			{
				const size_t dim_ = scv.m_midId[i].dim;
				const size_t id = scv.m_midId[i].id;
				UG_ASSERT(dim_ <= (size_t)dim, "Dimension wrong");
				UG_ASSERT(id < m_gloMid[dim_].size(), "id " << id << " in dim="
				          <<dim_<<" wrong. (size is "<< m_gloMid[dim_].size()<<")\n");
				scv.m_vGloPos[i] = m_gloMid[dim_][id];
			}
		}

		// i = number of side
		void compute_side_midpoints(size_t i,
		                            MathVector<dim>& locSideMid,
		                            MathVector<worldDim>& gloSideMid)
		{
			if(dim>1)
			{
				const size_t coID0 = m_rRefElem.id(2, i, 0, 0);
				locSideMid = m_locMid[0][coID0];
				gloSideMid = m_gloMid[0][coID0];

				// add corner coordinates of the corners of the geometric object
				for(size_t j = 1; j < m_rRefElem.num(2, i, 0); ++j)
				{
					const size_t coID = m_rRefElem.id(2, i, 0, j);
					locSideMid += m_locMid[0][coID];
					gloSideMid += m_gloMid[0][coID];
				}

				// scale for correct averaging
				locSideMid *= 1./(m_rRefElem.num(2, i, 0));
				gloSideMid *= 1./(m_rRefElem.num(2, i, 0));
			}
		}

		// i, j, k, l = number nodes
		void compute_side_midpoints(size_t i, size_t j, size_t k, size_t l,
									MathVector<dim>& locSideMid,
									MathVector<worldDim>& gloSideMid)
		{
			VecAdd(locSideMid, m_locMid[0][i], m_locMid[0][j], m_locMid[0][k], m_locMid[0][l]);
			VecAdd(gloSideMid, m_gloMid[0][i], m_gloMid[0][j], m_gloMid[0][k], m_gloMid[0][l]);

			// scale for correct averaging
			locSideMid *= 0.25;
			gloSideMid *= 0.25;
		}

		// i, j, k = number nodes
		void compute_side_midpoints(size_t i, size_t j, size_t k,
									MathVector<dim>& locSideMid, MathVector<worldDim>& gloSideMid)
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
			friend class DimHFV1Geometry<TDim, TWorldDim>;
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
                friend class DimHFV1Geometry<TDim, TWorldDim>;
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

    //  help array: maps NaturalEdge -> NodeOnEdge (-1 if no node on edge)
        std::vector<NatEdgeInfo> m_vNatEdgeInfo;
        std::vector<NewEdgeInfo> m_vNewEdgeInfo;

	private:
	///	pointer to current element
		GridObject* m_pElem;
		
	///	current reference object id
		ReferenceObjectID m_roid;

		std::vector<MathVector<dim> > m_locMid[dim+1];
		std::vector<MathVector<worldDim> > m_gloMid[dim+1];

	// 	SubControlVolumeFaces
		std::vector<SCVF> m_vSCVF;

	// 	SubControlVolumes
		std::vector<SCV> m_vSCV;

	//	number of shape functions
		size_t m_numSh;

	/// reference element
		DimReferenceElement<dim> m_rRefElem;

	/// reference mapping
		DimReferenceMapping<dim, worldDim>* m_rMapping;
		
};


////////////////////////////////////////////////////////////////////////////////
// Hanging FV1 Manifold Geometry
////////////////////////////////////////////////////////////////////////////////

template <typename TElem, int TWorldDim>
class HFV1ManifoldGeometry
{
	public:
	// 	type of element
		typedef TElem elem_type;

	// 	type of reference element
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	public:
	// 	order
		static const int order = 1;

	// number of natural bfs (eq. of SCV)
		static const size_t m_numNaturalBF = ref_elem_type::numCorners;

	// number of natural boundary face sides (eq. of SCVF)
		static const size_t m_numNaturalBFS = ref_elem_type::numEdges;

	// 	dimension of reference element
		static const int dim = ref_elem_type::dim;

	// 	dimension of world
		static const int worldDim = TWorldDim;

	// 	type of BoundaryFaces
		typedef typename hfv1_traits<ref_elem_type, dim>::scv_type bf_type;

	// 	Hanging node flag: this geometry supports hanging nodes
		static const bool usesHangingNodes = true;

	/// flag indicating if local data may change
		static const bool staticLocalData = true;

	protected:
		struct MidID
		{
				MidID() : dim(0), id(0) {};
				MidID(size_t dim_, size_t id_) : dim(dim_), id(id_) {};
				size_t dim;
				size_t id;
		};

	public:
		class BF
		{
			private:
			// 	let outer class access private members
				friend class HFV1ManifoldGeometry<TElem, TWorldDim>;

			// 	number of integration points
				static const size_t m_numIP = 1;

			// 	max number of corners of bf
				static const size_t numCorners = hfv1_traits<ref_elem_type, dim>::MaxNumCornersOfSCV;

			public:
				BF() {};

			/// node id that this bf is associated to
				inline size_t node_id() const {return nodeId;}

			/// number of integration points
				inline size_t num_ip() const {return m_numIP;}

			/// local integration point of scvf
				inline const MathVector<dim>& local_ip() const {return localIP;}

			/// global integration point of scvf
				inline const MathVector<worldDim>& global_ip() const {return globalIP;}

			/// volume of bf
				inline number volume() const {return vol;}

			/// number of corners, that bound the bf
				inline size_t num_corners() const {return numCorners;}

			/// return local position of corner number i
				inline const MathVector<dim>& local_corner(size_t i) const
					{UG_ASSERT(i < num_corners(), "Invalid index."); return m_vLocPos[i];}

			/// return global position of corner number i
				inline const MathVector<worldDim>& global_corner(size_t i) const
					{UG_ASSERT(i < num_corners(), "Invalid index."); return m_vGloPos[i];}

			/// number of shape functions
				inline size_t num_sh() const {return vShape.size();}

			/// value of shape function i in integration point
				inline number shape(size_t i, size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid index"); return vShape[i];}

			private:
				size_t nodeId;										// id of associated node

				// CORNERS: ordering is:
				// 1D: corner, centerOfElement
				// 2D: corner, side one, centerOfElement, side two
				MathVector<dim> m_vLocPos[numCorners];			// local position of node
				MathVector<worldDim> m_vGloPos[numCorners];	// global position of node

				//IPs & shapes
				MathVector<dim> localIP; // local integration point
				MathVector<worldDim> globalIP; // global integration point
				std::vector<number> vShape; // shapes at ip

				number vol;
		};

	protected:
		std::vector<MathVector<dim> > m_vLocBFIP;
		std::vector<MathVector<worldDim> > m_vGlobBFIP;

	public:
	/// constructor
		HFV1ManifoldGeometry();

	///	update data for given element
		void update(GridObject* elem, const MathVector<worldDim>* vCornerCoords,
		            const ISubsetHandler* ish = NULL);

	/// print information about hfvg to log
		void print();
/*
	/// get vector of corners for current element
		const MathVector<worldDim>* corners() const {return m_gloMid[0];}
*/
	/// number of BoundaryFaces
		inline size_t num_bf() const {return m_vBF.size();}

	/// const access to Boundary Face number i
		inline const BF& bf(size_t i) const
			{UG_ASSERT(i < num_bf(), "Invalid Index."); return m_vBF[i];}

	/// returns all ips of scvf as they appear in scv loop
		const MathVector<worldDim>* bf_global_ips() const {return &m_vGlobBFIP[0];}

	/// returns number of all scvf ips
		size_t num_bf_global_ips() const {return m_vGlobBFIP.size();}

	/// returns all ips of scvf as they appear in scv loop
		const MathVector<dim>* bf_local_ips() const {return &m_vLocBFIP[0];}

	/// returns number of all scvf ips
		size_t num_bf_local_ips() const {return m_vLocBFIP.size();}

	/// returns subset index
		int subset_index() const
		{
			if (m_ssi != -1) return m_ssi;
			UG_THROW("Subset index of geometry unknown.")
		}

	protected:
		void compute_side_midpoints(MathVector<dim>& locSideMid,
								   MathVector<worldDim>& gloSideMid)
		{
			if (worldDim > 1)
			{
				const size_t coID0 = m_rRefElem.id(2, 0, 0, 0);
				locSideMid = m_locMid[0][coID0];
				gloSideMid = m_gloMid[0][coID0];

				// add corner coordinates of the corners of the geometric object
				for (size_t j = 1; j < m_rRefElem.num(2, 0, 0); ++j)
				{
					const size_t coID = m_rRefElem.id(2, 0, 0, j);
					locSideMid += m_locMid[0][coID];
					gloSideMid += m_gloMid[0][coID];
				}

				// scale for correct averaging
				locSideMid *= 1./(m_rRefElem.num(2, 0, 0));
				gloSideMid *= 1./(m_rRefElem.num(2, 0, 0));
			}
		}

		// i, j, k, l = number nodes
		void compute_side_midpoints(size_t i, size_t j, size_t k, size_t l,
									MathVector<dim>& locSideMid,
									MathVector<worldDim>& gloSideMid)
		{
			VecAdd(locSideMid, m_locMid[0][i], m_locMid[0][j], m_locMid[0][k], m_locMid[0][l]);
			VecAdd(gloSideMid, m_gloMid[0][i], m_gloMid[0][j], m_gloMid[0][k], m_gloMid[0][l]);

			// scale for correct averaging
			locSideMid *= 0.25;
			gloSideMid *= 0.25;
		}

		// i, j, k = number nodes
		void compute_side_midpoints(size_t i, size_t j, size_t k,
									MathVector<dim>& locSideMid, MathVector<worldDim>& gloSideMid)
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
			for (size_t e = 0; e < m_vNatEdgeInfo[i].num_child_edges(); ++e)
			{
				const size_t childId = m_vNatEdgeInfo[i].child_edge(e);
				if (m_vNewEdgeInfo[childId].from() == co)
					return childId;
				if (m_vNewEdgeInfo[childId].to() == co)
					return childId;
			}
			return -1;
		}

	private:
		struct NewEdgeInfo
		{
			friend class HFV1ManifoldGeometry<TElem, TWorldDim>;

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
			friend class HFV1ManifoldGeometry<TElem, TWorldDim>;
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

	//  help array: maps NaturalEdge -> NodeOnEdge (-1 if no node on edge)
		std::vector<NatEdgeInfo> m_vNatEdgeInfo;
		std::vector<NewEdgeInfo> m_vNewEdgeInfo;

	private:
	// 	pointer to current element
		GridObject* m_pElem;

	// 	local and global geom object midpoints for each dimension
		std::vector<MathVector<dim> > m_locMid[dim+1];
		std::vector<MathVector<worldDim> > m_gloMid[dim+1];

	// 	BndFaces
		std::vector<BF> m_vBF;

	// 	Reference Mapping
		ReferenceMapping<ref_elem_type, worldDim> m_rMapping;

	// 	Reference Element
		const ref_elem_type& m_rRefElem;

	//	subset index of the element represented
		int m_ssi;
};


}

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__HANGING_FINITE_VOLUME_GEOMETRY__ */
