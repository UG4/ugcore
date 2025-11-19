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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FV1_GEOMETRY__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FV1_GEOMETRY__

// extern libraries
#include <cmath>
#include <map>
#include <vector>

// other ug4 modules
#include "common/common.h"

// library intern includes
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_disc/reference_element/reference_element.h"
#include "lib_disc/reference_element/reference_element_traits.h"
#include "lib_disc/reference_element/reference_mapping.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/local_finite_element/lagrange/lagrangep1.h"
#include "lib_disc/quadrature/gauss/gauss_quad.h"
#include "fv_util.h"
#include "fv_geom_base.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// FV1 Geometry for Reference Element Type
////////////////////////////////////////////////////////////////////////////////

// forward declaration
template <typename TElem, int TWorldDim, bool TCondensed> class FV1Geometry_gen;

/// Geometry and shape functions for 1st order Vertex-Centered Finite Volume
/**
 * The class provides the geometry and shape functions for 1st order Vertex-Centered
 * Finite Element Finite Volume method based on the Donald Diagrams.
 *
 * Cf. class FV1Geometry_gen for the implementation.
 * 
 * \tparam	TElem		Element type
 * \tparam	TWorldDim	(physical) world dimension
 */
template <typename TElem, int TWorldDim>
class FV1Geometry : public FV1Geometry_gen<TElem, TWorldDim, false> {};

/// Geometry and shape functions for 1st order Vertex-Centered Finite Volume
/**
 * The class provides the geometry and shape functions for 1st order Vertex-Centered
 * Finite Element Finite Volume method based on the Donald Diagrams.
 * 
 * This class shifts the subcontrol volume face integration
 * points to the edges. This allows to reduce the matrix pattern and to avoid positive
 * off-diagonal entries in some cases. (For ex., the discretization of the Laplacian on
 * a grid of rectangles retains results in the 5-point stencil.) However note that, in many
 * cases, this leads to the discretization order reduction.
 *
 * Cf. class FV1Geometry_gen for the implementation.
 * 
 * \tparam	TElem		Element type
 * \tparam	TWorldDim	(physical) world dimension
 */
template <typename TElem, int TWorldDim>
class FV1CondensedGeometry : public FV1Geometry_gen<TElem, TWorldDim, true> {};

/// Geometry and shape functions for 1st order Vertex-Centered Finite Volume
/**
 * The class provides the geometry and shape functions for 1st order Vertex-Centered
 * Finite Element Finite Volume method based on the Donald Diagrams.
 * 
 * The class provides an option (TCondensed) to shift the subcontrol volume face integration
 * points to the edges. This allows to reduce the matrix pattern and to avoid positive
 * off-diagonal entries in some cases. (For ex., the discretization of the Laplacian on
 * a grid of rectangles retains results in the 5-point stencil.) However note that, in many
 * cases, this leads to the discretization order reduction.
 * 
 * \tparam	TElem		Element type
 * \tparam	TWorldDim	(physical) world dimension
 * \tparam	TCondensed	if to shift the scvf ip's to midpoints of the edges
 */
template <typename TElem, int TWorldDim, bool TCondensed>
class FV1Geometry_gen : public FVGeometryBase
{
	public:
	///	type of element
	using elem_type = TElem;

	///	type of reference element
	using ref_elem_type = typename reference_element_traits<TElem>::reference_element_type;

	///	used traits
	using traits = fv1_traits<ref_elem_type, TWorldDim>;

	public:
	///	dimension of reference element
		static constexpr int dim = ref_elem_type::dim;

	///	dimension of world
		static constexpr int worldDim = TWorldDim;

	///	Hanging node flag: this Geometry does not support hanging nodes
		static constexpr bool usesHangingNodes = false;

	/// flag indicating if local data may change
		static constexpr bool staticLocalData = true;
		
	///	whether the scheme shifts the scvf ip's to midpoints of the edges
		static constexpr bool condensed_scvf_ips = TCondensed;

	public:
	///	order
		static constexpr int order = 1;

	///	number of SubControlVolumes
		static constexpr size_t numSCV = traits::numSCV;

	///	type of SubControlVolume
		using scv_type = typename traits::scv_type;

	///	number of SubControlVolumeFaces
		static constexpr size_t numSCVF = traits::numSCVF;

	///	type of Shape function used
		using local_shape_fct_set_type = LagrangeP1<ref_elem_type>;

	///	number of shape functions
		static constexpr size_t nsh = local_shape_fct_set_type::nsh;

	/// number of integration points
		static constexpr size_t nip = 1;

	public:
	///	Sub-Control Volume Face structure
	/**
	 * Each finite element is cut by several sub-control volume faces. The idea
	 * is that the "dual" skeleton formed by the sub control volume faces of
	 * all elements again gives rise to a regular mesh with closed
	 * (lipschitz-bounded) control volumes. The SCVF are the boundary of the
	 * control volume. In computation the flux over each SCVF must be the same
	 * in both directions over the face in order to guarantee the conservation
	 * property.
	 */
		class SCVF
		{
			public:
			///	Number of corners of scvf
				static constexpr size_t numCo =	traits::NumCornersOfSCVF;

			public:
				SCVF() = default;

			/// index of SubControlVolume on one side of the scvf
				inline size_t from() const {return From;}

			/// index of SubControlVolume on one side of the scvf
				inline size_t to() const {return To;}

			/// normal on scvf (points direction "from"->"to"). Norm is equal to area
				inline const MathVector<worldDim>& normal() const {return Normal;}

			/// number of integration points on scvf
				inline size_t num_ip() const {return nip;}

			/// local integration point of scvf
				inline const MathVector<dim>& local_ip() const {return localIP;}

			/// global integration point of scvf
				inline const MathVector<worldDim>& global_ip() const {return globalIP;}

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim,dim>& JTInv() const {return JtInv;}

			/// Determinant of Jacobian in integration point
				inline number detJ() const {return detj;}

			/// number of shape functions
				inline size_t num_sh() const {return nsh;}

			/// value of shape function i in integration point
				inline number shape(size_t sh) const {return vShape[sh];}

			/// vector of shape functions in ip point
				inline const number* shape_vector() const {return vShape;}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return vLocalGrad[sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector() const {return vLocalGrad;}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return vGlobalGrad[sh];}

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector() const {return vGlobalGrad;}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return numCo;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return global corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			private:
			// 	let outer class access private members
				friend class FV1Geometry_gen<TElem, TWorldDim, TCondensed>;

			// This scvf separates the scv with the ids given in "from" and "to"
			// The computed normal points in direction from->to
				size_t From, To;

			//	The normal on the SCVF pointing (from -> to)
				MathVector<worldDim> Normal; // normal (incl. area)

			// ordering is:
			// 1D: edgeMidPoint
			// 2D: edgeMidPoint, CenterOfElement
			// 3D: edgeMidPoint, Side #1, CenterOfElement, Side #2
				MathVector<dim> vLocPos[numCo]; // local corners of scvf
				MathVector<worldDim> vGloPos[numCo]; // global corners of scvf
				MidID vMidID[numCo]; // dimension and id of object, that's midpoint bounds the scvf

			// scvf part
				MathVector<dim> localIP; // local integration point
				MathVector<worldDim> globalIP; // global integration point

			// shapes and derivatives
				number vShape[nsh]; // shapes at ip
				MathVector<dim> vLocalGrad[nsh]; // local grad at ip
				MathVector<worldDim> vGlobalGrad[nsh]; // global grad at ip
				MathMatrix<worldDim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

	///	sub control volume structure
		class SCV
		{
			public:
			/// Number of corners of scvf
				static constexpr size_t numCo = traits::NumCornersOfSCV;

			public:
				SCV() {};

			/// volume of scv
				inline number volume() const {return Vol;}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return numCo;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return global corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			/// return local corners
				inline const MathVector<dim>* local_corners() const
					{return &vLocPos[0];}

			/// return global corners
				inline const MathVector<worldDim>* global_corners() const
					{return &vGloPos[0];}

			/// node id that this scv is associated to
				inline size_t node_id() const {return nodeId;}

			/// number of integration points
				inline size_t num_ip() const {return nip;}

			/// local integration point of scv
				inline const MathVector<dim>& local_ip() const {return vLocPos[0];}

			/// global integration point
				inline const MathVector<worldDim>& global_ip() const {return vGloPos[0];}

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim,dim>& JTInv() const {return JtInv;}

			/// Determinant of Jacobian in integration point
				inline number detJ() const {return detj;}

			/// number of shape functions
				inline size_t num_sh() const {return nsh;}

			/// value of shape function i in integration point
				inline number shape(size_t sh) const {return vShape[sh];}

			/// vector of shape functions in ip point
				inline const number* shape_vector() const {return vShape;}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return vLocalGrad[sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector() const {return vLocalGrad;}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return vGlobalGrad[sh];}

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector() const {return vGlobalGrad;}

			private:
			// 	let outer class access private members
				friend class FV1Geometry_gen<TElem, TWorldDim, TCondensed>;

			//  node id of associated node
				size_t nodeId;

			//	volume of scv
				number Vol;

			//	local and global positions of this element
				MathVector<dim> vLocPos[numCo]; // local position of node
				MathVector<worldDim> vGloPos[numCo]; // global position of node
				MidID midId[numCo]; // dimension and id of object, that's midpoint bounds the scv

			// shapes and derivatives
				number vShape[nsh]; // shapes at ip
				MathVector<dim> vLocalGrad[nsh]; // local grad at ip
				MathVector<worldDim> vGlobalGrad[nsh]; // global grad at ip
				MathMatrix<worldDim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

	///	boundary face
		class BF
		{
			public:
			/// Number of corners of bf
				static constexpr size_t numCo =	traits::NumCornersOfBF;

			public:
				BF() = default;

			/// index of SubControlVolume of the bf
				inline size_t node_id() const {return nodeId;}

			/// number of integration points on bf
				inline size_t num_ip() const {return nip;}

			/// local integration point of bf
				inline const MathVector<dim>& local_ip() const {return localIP;}

			/// global integration point of bf
				inline const MathVector<worldDim>& global_ip() const {return globalIP;}

			/// outer normal on bf. Norm is equal to area
				inline const MathVector<worldDim>& normal() const {return Normal;} // includes area

			/// volume of bf
				inline number volume() const {return Vol;}

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim, dim>& JTInv() const {return JtInv;}

			/// Determinant of Jacobian in integration point
				inline number detJ() const {return detj;}

			/// number of shape functions
				inline size_t num_sh() const {return nsh;}

			/// value of shape function i in integration point
				inline number shape(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return vShape[sh];}

			/// vector of local gradients in ip point
				inline const number* shape_vector() const {return vShape;}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return vLocalGrad[sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector() const {return vLocalGrad;}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return vGlobalGrad[sh];}

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector() const {return vGlobalGrad;}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return numCo;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return global corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			private:
			/// let outer class access private members
				friend class FV1Geometry_gen<TElem, TWorldDim, TCondensed>;

			// 	id of scv this bf belongs to
				size_t nodeId;

			// ordering is:
			// 1D: edgeMidPoint
			// 2D: edgeMidPoint, CenterOfElement
			// 3D: edgeMidPoint, Side one, CenterOfElement, Side two
				MathVector<dim> vLocPos[numCo]; // local corners of bf
				MathVector<worldDim> vGloPos[numCo]; // global corners of bf
				MidID vMidID[numCo]; // dimension and id of object, that's midpoint bounds the bf

			// 	scvf part
				MathVector<dim> localIP; // local integration point
				MathVector<worldDim> globalIP; // global integration point
				MathVector<worldDim> Normal; // normal (incl. area)
				number Vol; // volume of bf

			// 	shapes and derivatives
				number vShape[nsh]; // shapes at ip
				MathVector<dim> vLocalGrad[nsh]; // local grad at ip
				MathVector<worldDim> vGlobalGrad[nsh]; // global grad at ip
				MathMatrix<worldDim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

	public:
	/// construct object and initialize local values and sizes
		FV1Geometry_gen();

	///	update local data
		void update_local_data();

	/// update data for given element
		void update(GridObject* elem, const MathVector<worldDim>* vCornerCoords,
		            const ISubsetHandler* ish = nullptr);

	/// update boundary data for given element
		void update_boundary_faces(GridObject* elem,
		                           const MathVector<worldDim>* vCornerCoords,
		                           const ISubsetHandler* ish = nullptr);

	///	get the element
		TElem* elem() const {return m_pElem;}
		
	/// get vector of the global coordinates of corners for current element
		const MathVector<worldDim>* corners() const {return m_vvGloMid[0];}

	/// number of SubControlVolumeFaces
		size_t num_scvf() const {return numSCVF;};

	/// const access to SubControlVolumeFace number i
		const SCVF& scvf(size_t i) const
			{UG_ASSERT(i < num_scvf(), "Invalid Index."); return m_vSCVF[i];}

	/// number of SubControlVolumes
		// do not use this method to obtain the number of shape functions,
		// since this is NOT the same for pyramids; use num_sh() instead.
		size_t num_scv() const {return numSCV;}

	/// const access to SubControlVolume number i
		const SCV& scv(size_t i) const
			{UG_ASSERT(i < num_scv(), "Invalid Index."); return m_vSCV[i];}

	/// number of shape functions
		size_t num_sh() const {return nsh;};

	///	returns reference object id
		ReferenceObjectID roid() const {return ref_elem_type::REFERENCE_OBJECT_ID;}


	public:
	/// returns number of all scvf ips
		size_t num_scvf_ips() const {return numSCVF;}

	/// returns all ips of scvf as they appear in scv loop
		const MathVector<dim>* scvf_local_ips() const {return m_vLocSCVF_IP;}

	/// returns all ips of scvf as they appear in scv loop
		const MathVector<worldDim>* scvf_global_ips() const {return m_vGlobSCVF_IP;}

	/// returns number of all scv ips
		size_t num_scv_ips() const {return numSCV;}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<dim>* scv_local_ips() const {
			if(ref_elem_type::REFERENCE_OBJECT_ID == ROID_PYRAMID || ref_elem_type::REFERENCE_OBJECT_ID == ROID_OCTAHEDRON)
				return &(m_vLocSCV_IP[0]);
			else
				return &(m_vvLocMid[0][0]);
		}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<worldDim>* scv_global_ips() const {
			if(ref_elem_type::REFERENCE_OBJECT_ID == ROID_PYRAMID || ref_elem_type::REFERENCE_OBJECT_ID == ROID_OCTAHEDRON)
				return &(m_vGlobSCV_IP[0]);
			else
				return &(m_vvGloMid[0][0]);
		}

	/// return local coords for node ID
		const MathVector<dim>& local_node_position(size_t nodeID) const
		{
			UG_ASSERT(nodeID < (size_t) ref_elem_type::numCorners, "Invalid node id.");
			return m_vvLocMid[0][nodeID];
		}

	/// return global coords for node ID
		const MathVector<worldDim>& global_node_position(size_t nodeID) const
		{
			UG_ASSERT(nodeID < (size_t) ref_elem_type::numCorners, "Invalid node id.");
			return m_vvGloMid[0][nodeID];
		}

	///	returns the local coordinates of the center of mass of the element
		const MathVector<dim>* coe_local() const {return &(m_vvLocMid[dim][0]);}
		
	///	returns the global coordinates of the center of mass of the element
		const MathVector<worldDim>* coe_global() const {return &(m_vvGloMid[dim][0]);}

	protected:
	//	global and local ips on SCVF
		MathVector<worldDim> m_vGlobSCVF_IP[numSCVF];
		MathVector<dim> m_vLocSCVF_IP[numSCVF];

	//	global and local ips on SCV (only needed for Pyramid and Octahedron)
		MathVector<worldDim> m_vGlobSCV_IP[numSCV];
		MathVector<dim> m_vLocSCV_IP[numSCV];

	public:
	/// add subset that is interpreted as boundary subset.
		inline void add_boundary_subset(int subsetIndex) {m_mapVectorBF[subsetIndex];}

	/// removes subset that is interpreted as boundary subset.
		inline void remove_boundary_subset(int subsetIndex) {m_mapVectorBF.erase(subsetIndex);}

	/// reset all boundary subsets
		inline void clear_boundary_subsets() {m_mapVectorBF.clear();}

	/// number of registered boundary subsets
		inline size_t num_boundary_subsets() {return m_mapVectorBF.size();}

	/// number of all boundary faces
		inline size_t num_bf() const
		{
			typename std::map<int, std::vector<BF> >::const_iterator it;
			size_t num = 0;
			for ( it=m_mapVectorBF.begin() ; it != m_mapVectorBF.end(); it++ )
				num += (*it).second.size();
			return num;
		}

	/// number of boundary faces on subset 'subsetIndex'
		inline size_t num_bf(int si) const
		{
			typename std::map<int, std::vector<BF> >::const_iterator it;
			it = m_mapVectorBF.find(si);
			if(it == m_mapVectorBF.end()) return 0;
			else return (*it).second.size();
		}

	/// returns the boundary face i for subsetIndex
		inline const BF& bf(int si, size_t i) const
		{
			UG_ASSERT(i < num_bf(si), "Invalid index.");
			typename std::map<int, std::vector<BF> >::const_iterator it;
			it = m_mapVectorBF.find(si);
			if(it == m_mapVectorBF.end()) UG_THROW("FVGeom: No bnd face for subset"<<si);
			return (*it).second[i];
		}

	/// returns reference to vector of boundary faces for subsetIndex
		inline const std::vector<BF>& bf(int si) const
		{
			typename std::map<int, std::vector<BF> >::const_iterator it;
			it = m_mapVectorBF.find(si);
			if(it == m_mapVectorBF.end()) return m_vEmptyVectorBF;
			return (*it).second;
		}

		void reset_curr_elem() {m_pElem = nullptr;}

	protected:
		std::map<int, std::vector<BF> > m_mapVectorBF;
		std::vector<BF> m_vEmptyVectorBF;

	private:
	///	pointer to current element
		TElem* m_pElem;

	///	max number of geom objects in all dimensions
	// 	(most objects in 1 dim, i.e. number of edges, but +1 for 1D)
		static constexpr int maxMid = numSCVF + 1;

	///	local and global geom object midpoints for each dimension
		MathVector<dim> m_vvLocMid[dim+1][maxMid];
		MathVector<worldDim> m_vvGloMid[dim+1][maxMid];

	///	SubControlVolumeFaces
		SCVF m_vSCVF[numSCVF];

	///	SubControlVolumes
		SCV m_vSCV[numSCV];

	///	Reference Mapping
		ReferenceMapping<ref_elem_type, worldDim> m_mapping;

	///	Reference Element
		const ref_elem_type& m_rRefElem;

	///	Shape function set
		const local_shape_fct_set_type& m_rTrialSpace;
};

////////////////////////////////////////////////////////////////////////////////
// Dim-dependent Finite Volume Geometry
////////////////////////////////////////////////////////////////////////////////

/// Geometry and shape functions for 1st order Vertex-Centered Finite Volume
/**
 * \tparam	TDim		reference element dim
 * \tparam	TWorldDim	(physical) world dimension
 */
template <int TDim, int TWorldDim = TDim>
class DimFV1Geometry : public FVGeometryBase
{
	public:
	///	used traits
		using traits = fv1_dim_traits<TDim, TWorldDim>;

	public:
	///	dimension of reference element
		static constexpr int dim = TDim;

	///	dimension of world
		static constexpr int worldDim = TWorldDim;

	///	Hanging node flag: this Geometry does not support hanging nodes
		static constexpr bool usesHangingNodes = false;

	/// flag indicating if local data may change
		static constexpr bool staticLocalData = false;

	public:
	///	order
		static constexpr int order = 1;

	///	number of SubControlVolumes
		static constexpr size_t maxNumSCV = traits::maxNumSCV;

	///	type of SubControlVolume
		using scv_type = typename traits::scv_type;

	///	max number of SubControlVolumeFaces
		static constexpr size_t maxNumSCVF = traits::maxNumSCVF;

	/// max number of shape functions
		static constexpr size_t maxNSH = traits::maxNSH;

	/// number of integration points
		static constexpr size_t nip = 1;

	public:
	///	Sub-Control Volume Face structure
	/**
	 * Each finite element is cut by several sub-control volume faces. The idea
	 * is that the "dual" skeleton formed by the sub control volume faces of
	 * all elements again gives rise to a regular mesh with closed
	 * (lipschitz-bounded) control volumes. The SCVF are the boundary of the
	 * control volume. In computation the flux over each SCVF must be the same
	 * in both directions over the face in order to guarantee the conservation
	 * property.
	 */
		class SCVF
		{
			public:
			///	Number of corners of scvf
				static constexpr size_t numCo = traits::NumCornersOfSCVF;

			public:
				SCVF() {}

			/// index of SubControlVolume on one side of the scvf
				inline size_t from() const {return From;}

			/// index of SubControlVolume on one side of the scvf
				inline size_t to() const {return To;}

			/// number of integration points on scvf
				inline size_t num_ip() const {return nip;}

			/// local integration point of scvf
				inline const MathVector<dim>& local_ip() const {return localIP;}

			/// global integration point of scvf
				inline const MathVector<worldDim>& global_ip() const {return globalIP;}

			/// normal on scvf (points direction "from"->"to"). Norm is equal to area
				inline const MathVector<worldDim>& normal() const {return Normal;}

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim,dim>& JTInv() const {return JtInv;}

			/// Determinant of Jacobian in integration point
				inline number detJ() const {return detj;}

			/// number of shape functions
				inline size_t num_sh() const {return numSH;}

			/// value of shape function i in integration point
				inline number shape(size_t sh) const {return vShape[sh];}

			/// vector of shape functions in ip point
				inline const number* shape_vector() const {return vShape;}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return vLocalGrad[sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector() const {return vLocalGrad;}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return vGlobalGrad[sh];}

			/// vector of gloabl gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector() const {return vGlobalGrad;}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return numCo;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			private:
			// 	let outer class access private members
				friend class DimFV1Geometry<dim, worldDim>;

			// This scvf separates the scv with the ids given in "from" and "to"
			// The computed normal points in direction from->to
				size_t From, To;

			//	The normal on the SCVF pointing (from -> to)
				MathVector<worldDim> Normal; // normal (incl. area)

			// ordering is:
			// 1D: edgeMidPoint
			// 2D: edgeMidPoint, CenterOfElement
			// 3D: edgeMidPoint, Side one, CenterOfElement, Side two
				MathVector<dim> vLocPos[numCo]; // local corners of scvf
				MathVector<worldDim> vGloPos[numCo]; // global corners of scvf
				MidID vMidID[numCo]; // dimension and id of object, that's midpoint bounds the scvf

			// scvf part
				MathVector<dim> localIP; // local integration point
				MathVector<worldDim> globalIP; // global integration point

			// shapes and derivatives
				size_t numSH;
				number vShape[maxNSH]; // shapes at ip
				MathVector<dim> vLocalGrad[maxNSH]; // local grad at ip
				MathVector<worldDim> vGlobalGrad[maxNSH]; // global grad at ip
				MathMatrix<worldDim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

	///	sub control volume structure
		class SCV
		{
			public:
			/// Number of corners of scv
				static constexpr size_t numCo = traits::NumCornersOfSCV;

			public:
				SCV() = default;

			/// volume of scv
				inline number volume() const {return Vol;}

			/// number of corners, that bound the scv
				inline size_t num_corners() const {return numCo;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return global corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			/// node id that this scv is associated to
				inline size_t node_id() const {return nodeId;}

			/// number of integration points
				inline size_t num_ip() const {return nip;}

			/// local integration point of scv
				inline const MathVector<dim>& local_ip() const {return vLocPos[0];}

			/// global integration point
				inline const MathVector<worldDim>& global_ip() const {return vGloPos[0];}

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim,dim>& JTInv() const {return JtInv;}

			/// Determinant of Jacobian in integration point
				inline number detJ() const {return detj;}

			/// number of shape functions
				inline size_t num_sh() const {return numSH;}

			/// value of shape function i in integration point
				inline number shape(size_t sh) const {return vShape[sh];}

			/// vector of shape functions in ip point
				inline const number* shape_vector() const {return vShape;}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return vLocalGrad[sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector() const {return vLocalGrad;}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return vGlobalGrad[sh];}

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector() const {return vGlobalGrad;}

			private:
			// 	let outer class access private members
				friend class DimFV1Geometry<dim, worldDim>;

			//  node id of associated node
				size_t nodeId;

			//	volume of scv
				number Vol;

			//	local and global positions of this element
				MathVector<dim> vLocPos[numCo]; // local position of node
				MathVector<worldDim> vGloPos[numCo]; // global position of node
				MidID vMidID[numCo]; // dimension and id of object, that's midpoint bounds the scv

			// shapes and derivatives
				size_t numSH;
				number vShape[maxNSH]; // shapes at ip
				MathVector<dim> vLocalGrad[maxNSH]; // local grad at ip
				MathVector<worldDim> vGlobalGrad[maxNSH]; // global grad at ip
				MathMatrix<worldDim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

	///	boundary face
		class BF
		{
			public:
			/// Number of corners of bf
				static constexpr size_t numCo = traits::NumCornersOfSCVF;

			public:
				BF() = default;

			/// index of SubControlVolume of the bf
				inline size_t node_id() const {return nodeId;}

			/// number of integration points on bf
				inline size_t num_ip() const {return nip;}

			/// local integration point of bf
				inline const MathVector<dim>& local_ip() const {return localIP;}

			/// global integration point of bf
				inline const MathVector<worldDim>& global_ip() const {return globalIP;}

			/// outer normal on bf. Norm is equal to area
				inline const MathVector<worldDim>& normal() const {return Normal;} // includes area

			/// volume of bf
				inline number volume() const {return Vol;}

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim, dim>& JTInv() const {return JtInv;}

			/// Determinant of Jacobian in integration point
				inline number detJ() const {return detj;}

			/// number of shape functions
				inline size_t num_sh() const {return numSH;}

			/// value of shape function i in integration point
				inline number shape(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return vShape[sh];}

			/// vector of local gradients in ip point
				inline const number* shape_vector() const {return vShape;}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return vLocalGrad[sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector() const {return vLocalGrad;}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return vGlobalGrad[sh];}

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector() const {return vGlobalGrad;}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return numCo;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return global corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			private:
			/// let outer class access private members
				friend class DimFV1Geometry<dim, worldDim>;

			// 	id of scv this bf belongs to
				size_t nodeId;

			// ordering is:
			// 1D: edgeMidPoint
			// 2D: edgeMidPoint, CenterOfElement
			// 3D: edgeMidPoint, Side one, CenterOfElement, Side two
				MathVector<dim> vLocPos[numCo]; // local corners of bf
				MathVector<worldDim> vGloPos[numCo]; // global corners of bf
				MidID vMidID[numCo]; // dimension and id of object, that's midpoint bounds the bf

			// 	scvf part
				MathVector<dim> localIP; // local integration point
				MathVector<worldDim> globalIP; // global integration point
				MathVector<worldDim> Normal; // normal (incl. area)
				number Vol; // volume of bf

			// 	shapes and derivatives
				size_t numSH;
				number vShape[maxNSH]; // shapes at ip
				MathVector<dim> vLocalGrad[maxNSH]; // local grad at ip
				MathVector<worldDim> vGlobalGrad[maxNSH]; // global grad at ip
				MathMatrix<worldDim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

	public:
	/// construct object and initialize local values and sizes
		DimFV1Geometry() : m_pElem(nullptr), m_roid(ROID_UNKNOWN) {};

	/// update data for given element
		void update(GridObject* elem, const MathVector<worldDim>* vCornerCoords,
		            const ISubsetHandler* ish = nullptr);

	/// update boundary data for given element
		void update_boundary_faces(GridObject* elem,
		                           const MathVector<worldDim>* vCornerCoords,
		                           const ISubsetHandler* ish = nullptr);

	///	get the element
		GridObject* elem() const {return m_pElem;}
		
	/// get vector of corners for current element
		const MathVector<worldDim>* corners() const {return m_vvGloMid[0];}

	/// number of SubControlVolumeFaces
		size_t num_scvf() const {return m_numSCVF;};

	/// const access to SubControlVolumeFace number i
		const SCVF& scvf(size_t i) const
			{UG_ASSERT(i < num_scvf(), "Invalid Index."); return m_vSCVF[i];}

	/// number of SubControlVolumes
		// do not use this method to obtain the number of shape functions,
		// since this is NOT the same for pyramids; use num_sh() instead.
		size_t num_scv() const {return m_numSCV;}

	/// const access to SubControlVolume number i
		const SCV& scv(size_t i) const
			{UG_ASSERT(i < num_scv(), "Invalid Index."); return m_vSCV[i];}

	/// number of shape functions
		size_t num_sh() const {return m_nsh;};

	public:
	/// returns number of all scvf ips
		size_t num_scvf_ips() const {return m_numSCVF;}

	/// returns all ips of scvf as they appear in scv loop
		const MathVector<dim>* scvf_local_ips() const {return m_vLocSCVF_IP;}

	/// returns all ips of scvf as they appear in scv loop
		const MathVector<worldDim>* scvf_global_ips() const {return m_vGlobSCVF_IP;}

	/// returns number of all scv ips
		size_t num_scv_ips() const {return m_numSCV;}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<dim>* scv_local_ips() const {
			if(m_roid == ROID_PYRAMID || m_roid == ROID_OCTAHEDRON)
				return &(m_vLocSCV_IP[0]);
			else
				return &(m_vvLocMid[0][0]);
		}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<worldDim>* scv_global_ips() const {
			if(m_roid == ROID_PYRAMID || m_roid == ROID_OCTAHEDRON)
				return &(m_vGlobSCV_IP[0]);
			else
				return &(m_vvGloMid[0][0]);
		}

	/// return local coords for node ID
		const MathVector<dim>& local_node_position(size_t nodeID) const
		{
			UG_ASSERT(nodeID < (size_t) maxMid, "Invalid node id.");
			return m_vvLocMid[0][nodeID];
		}

	/// return global coords for node ID
		const MathVector<worldDim>& global_node_position(size_t nodeID) const
		{
			UG_ASSERT(nodeID < (size_t) maxMid, "Invalid node id.");
			return m_vvGloMid[0][nodeID];
		}

	///	returns the local coordinates of the center of mass of the element
		const MathVector<dim>* coe_local() const {return &(m_vvLocMid[dim][0]);}
		
	///	returns the global coordinates of the center of mass of the element
		const MathVector<worldDim>* coe_global() const {return &(m_vvGloMid[dim][0]);}

	///	returns reference object id
		ReferenceObjectID roid() const {return m_roid;}

	///	update local data
		void update_local(ReferenceObjectID roid);

	protected:
	//	global and local ips on SCVF
		MathVector<worldDim> m_vGlobSCVF_IP[maxNumSCVF];
		MathVector<dim> m_vLocSCVF_IP[maxNumSCVF];

	//	global and local ips on SCV (only needed for Pyramid and Octahedron)
		MathVector<worldDim> m_vGlobSCV_IP[maxNumSCV];
		MathVector<dim> m_vLocSCV_IP[maxNumSCV];

	public:
	/// add subset that is interpreted as boundary subset.
		inline void add_boundary_subset(int subsetIndex) {m_mapVectorBF[subsetIndex];}

	/// removes subset that is interpreted as boundary subset.
		inline void remove_boundary_subset(int subsetIndex) {m_mapVectorBF.erase(subsetIndex);}

	/// reset all boundary subsets
		inline void clear_boundary_subsets() {m_mapVectorBF.clear();}

	/// number of registered boundary subsets
		inline size_t num_boundary_subsets() {return m_mapVectorBF.size();}

	/// number of all boundary faces
		inline size_t num_bf() const
		{
			typename std::map<int, std::vector<BF> >::const_iterator it;
			size_t num = 0;
			for ( it=m_mapVectorBF.begin() ; it != m_mapVectorBF.end(); it++ )
				num += (*it).second.size();
			return num;
		}

	/// number of boundary faces on subset 'subsetIndex'
		inline size_t num_bf(int si) const
		{
			typename std::map<int, std::vector<BF> >::const_iterator it;
			it = m_mapVectorBF.find(si);
			if(it == m_mapVectorBF.end()) return 0;
			else return (*it).second.size();
		}

	/// returns the boundary face i for subsetIndex
		inline const BF& bf(int si, size_t i) const
		{
			typename std::map<int, std::vector<BF> >::const_iterator it;
			it = m_mapVectorBF.find(si);
			if(it == m_mapVectorBF.end()) UG_THROW("DimFV1Geom: No BndSubset "<<si);
			return (*it).second[i];
		}

	/// returns reference to vector of boundary faces for subsetIndex
		inline const std::vector<BF>& bf(int si) const
		{
			typename std::map<int, std::vector<BF> >::const_iterator it;
			it = m_mapVectorBF.find(si);
			if(it == m_mapVectorBF.end()) return m_vEmptyVectorBF;
			return (*it).second;
		}

	protected:
		std::map<int, std::vector<BF> > m_mapVectorBF;
		std::vector<BF> m_vEmptyVectorBF;

	private:
	///	pointer to current element
		GridObject* m_pElem;

	///	current reference object id
		ReferenceObjectID m_roid;

	///	current number of scv
		size_t m_numSCV;

	///	current number of scvf
		size_t m_numSCVF;

	/// current number of shape functions
		size_t m_nsh;

	///	max number of geometric objects in a dimension
	// 	(most objects in 1 dim, i.e. number of edges, but +1 for 1D)
		static constexpr int maxMid = maxNumSCVF + 1;

	///	local and global geom object midpoints for each dimension
		MathVector<dim> m_vvLocMid[dim+1][maxMid];
		MathVector<worldDim> m_vvGloMid[dim+1][maxMid];

	///	SubControlVolumeFaces
		SCVF m_vSCVF[maxNumSCVF];

	///	SubControlVolumes
		SCV m_vSCV[maxNumSCV];
};


////////////////////////////////////////////////////////////////////////////////
// FV1 Manifold Geometry
////////////////////////////////////////////////////////////////////////////////

template <typename TElem, int TWorldDim>
class FV1ManifoldGeometry
{
	public:
	// 	type of element
	using elem_type = TElem;

	// 	type of reference element
	using ref_elem_type = typename reference_element_traits<TElem>::reference_element_type;

	public:
	// 	order
		static constexpr int order = 1;

	// 	number of BoundaryFaces
		static constexpr size_t m_numBF = ref_elem_type::numCorners;
		
	// 	dimension of reference element
		static constexpr int dim = ref_elem_type::dim;
		
	// 	dimension of world
		static constexpr int worldDim = TWorldDim;

	// 	type of BoundaryFaces
		using bf_type = typename fv1_traits<ref_elem_type, dim>::scv_type ;

	// 	Hanging node flag: this Geometry does not support hanging nodes
		static constexpr bool usesHangingNodes = false;

	/// flag indicating if local data may change
		static constexpr bool staticLocalData = true;

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
				friend class FV1ManifoldGeometry<TElem, TWorldDim>;

			// 	number of integration points
				static constexpr size_t m_numIP = 1;
				
			// 	max number of corners of bf
				static constexpr size_t numCorners = fv1_traits<ref_elem_type, dim>::NumCornersOfSCV;

			public:
				BF() = default;

			/// node id that this bf is associated to
				inline size_t node_id() const {return nodeId;}

			/// number of integration points
				inline size_t num_ip() const {return m_numIP;}

			/// local integration point of bf
				inline const MathVector<dim>& local_ip() const
				{return vLocPos[0];}	// <-- always the vertex

			/// global integration point
				inline const MathVector<worldDim>& global_ip() const
				{return vGloPos[0];}	// <-- here too

			/// volume of bf
				inline number volume() const {return vol;}

			/// number of corners, that bound the bf
				inline size_t num_corners() const {return numCorners;}

			/// return local position of corner number i
				inline const MathVector<dim>& local_corner(size_t i) const
					{UG_ASSERT(i < num_corners(), "Invalid index."); return vLocPos[i];}

			/// return global position of corner number i
				inline const MathVector<worldDim>& global_corner(size_t i) const
					{UG_ASSERT(i < num_corners(), "Invalid index."); return vGloPos[i];}
				
			/// number of shape functions
				inline size_t num_sh() const {return vShape.size();}

			/// value of shape function i in integration point
				inline number shape(size_t i, size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid index"); return vShape[i];}
				
			private:
				size_t nodeId;										// id of associated node
				
				// CORNERS: ordering is:
				// 1D: edgeMidPoint, CenterOfElement
				// 2D: edgeMidPoint, Side one, CenterOfElement, Side two
				MathVector<dim> vLocPos[numCorners];			// local position of node
				MathVector<worldDim> vGloPos[numCorners];	// global position of node
				MidID midId[numCorners];			// dimension and id of object, whose midpoint bounds the scv
				
				// IPs & shapes
				std::vector<number> vShape; // shapes at ip
				
				number vol;
		};
	
	protected:
		void copy_local_corners(BF& bf)
		{
			for (size_t i = 0; i < bf.num_corners(); ++i)
			{
				const size_t dim = bf.midId[i].dim;
				const size_t id = bf.midId[i].id;
				bf.vLocPos[i] = m_locMid[dim][id];
			}
		}
		
		void copy_global_corners(BF& bf)
		{
			for (size_t i = 0; i < bf.num_corners(); ++i)
			{
				const size_t dim = bf.midId[i].dim;
				const size_t id = bf.midId[i].id;
				bf.vGloPos[i] = m_gloMid[dim][id];
			}
		}

		std::vector<MathVector<dim> > m_vLocBFIP;
		std::vector<MathVector<worldDim> > m_vGlobBFIP;
		
		
	public:
	/// constructor
		FV1ManifoldGeometry();
		
	///	update data for given element
		void update(GridObject* elem, const MathVector<worldDim>* vCornerCoords,
		            const ISubsetHandler* ish = nullptr);
			
	///	get the element
		GridObject* elem() const {return m_pElem;}
		
	/// get vector of corners for current element
		const MathVector<worldDim>* corners() const {return m_gloMid[0];}

	/// number of BoundaryFaces
		inline size_t num_bf() const {return m_numBF;}

	/// const access to Boundary Face number i
		inline const BF& bf(size_t i) const
			{UG_ASSERT(i < num_bf(), "Invalid Index."); return m_vBF[i];}
	
	/// returns all ips of scvf as they appear in scv loop
		const MathVector<worldDim>* bf_global_ips() const {return &m_vGlobBFIP[0];}

	/// returns number of all BF ips
		size_t num_bf_global_ips() const {return m_vGlobBFIP.size();}

	/// returns all ips of BF as they appear in scv loop
		const MathVector<dim>* bf_local_ips() const {return &m_vLocBFIP[0];}

	/// returns number of all BF ips
		size_t num_bf_local_ips() const {return m_vLocBFIP.size();}

	private:
	// 	pointer to current element
		GridObject* m_pElem;
		
	// 	local and global geom object midpoints for each dimension
		MathVector<dim> m_locMid[dim+1][m_numBF];
		MathVector<worldDim> m_gloMid[dim+1][m_numBF];
		
	// 	BndFaces
		BF m_vBF[m_numBF];
		
	// 	Reference Mapping
		ReferenceMapping<ref_elem_type, worldDim> m_rMapping;

	// 	Reference Element
		const ref_elem_type& m_rRefElem;
};

}

#endif