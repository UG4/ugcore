/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_VOLUME_HIGHER_ORDER_GEOMETRY__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_VOLUME_HIGHER_ORDER_GEOMETRY__

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
// #include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/quadrature/gauss/gauss_quad.h"
#include "fv_util.h"
#include "fv_geom_base.h"

namespace ug {

////////////////////////////////////////////////////////////////////////////////
// FV Geometry for Reference Element Type (all orders)
////////////////////////////////////////////////////////////////////////////////

/// Geometry and shape functions for any order Vertex-Centered Finite Volume
/**
 * \tparam	TOrder			order
 * \tparam	TElem			Element type
 * \tparam	TWorldDim		(physical) world dimension
 * \tparam	TQuadOrder		order for quadrature

 */
template <	int TOrder, typename TElem, int TWorldDim,
			int TQuadOrder = TOrder + 1>
class FVGeometry : public FVGeometryBase
{
	private:
	///	small abbreviation for order
		static constexpr int p = TOrder;

	public:
	///	type of element
		using elem_type = TElem;

	///	type of reference element
		using ref_elem_type = typename reference_element_traits<TElem>::reference_element_type;

	///	dimension of reference element
		static constexpr int dim = ref_elem_type::dim;

	///	dimension of world
		static constexpr int worldDim = TWorldDim;

	public:
	///	order
		static constexpr int order = TOrder;

	///	number of subelements
		static constexpr size_t numSubElem = Pow<p, dim>::value;

	///	traits used
		using traits = fv1_traits<ref_elem_type, worldDim>;

	///	type of SubControlVolumeFace
		using scvf_type = typename traits::scvf_type;

	///	number of SCVF per SubElement
		static constexpr size_t numSCVFPerSubElem = ref_elem_type::numEdges;

	///	number of SubControlVolumeFaces
		static constexpr size_t numSCVF = numSubElem * numSCVFPerSubElem;

	///	quadrature order
		static constexpr int quadOrderSCVF = TQuadOrder;

	///	type of quadrature rule
		using scvf_quad_rule_type = GaussQuadrature<scvf_type, quadOrderSCVF>;

	///	number of scvf ip
		static constexpr size_t numSCVFIP = scvf_quad_rule_type::nip * numSCVF;

	///	type of SubControlVolume
		using scv_type = typename traits::scv_type;

	///	number of SCV per SubElement
		static constexpr size_t numSCVPerSubElem = ref_elem_type::numCorners;

	///	number of SubControlVolumes
		static constexpr size_t numSCV = numSubElem * numSCVPerSubElem;

	///	quadrature order
		static constexpr int quadOrderSCV = TQuadOrder;

	///	type of quadrature rule
		using scv_quad_rule_type = GaussQuadrature<scv_type, quadOrderSCV>;

	///	number of scv ip
		static constexpr size_t numSCVIP = scv_quad_rule_type::nip * numSCV;

	/// type of Shape function used
		using local_shape_fct_set_type = LagrangeLSFS<ref_elem_type, p>;

	///	number of shape functions
		static constexpr size_t nsh = local_shape_fct_set_type::nsh;

	///	Hanging node flag: this Geometry does not support hanging nodes
		static constexpr bool usesHangingNodes = false;

	/// flag indicating if local data may change
		static constexpr bool staticLocalData = true;

	public:
	///	Sub-Control Volume Face structure
		class SCVF
		{
			public:
			///	Number of integration points
				static constexpr size_t nip = scvf_quad_rule_type::nip;

			///	Number of corners of scvf
				static constexpr size_t numCo = traits::NumCornersOfSCVF;

			public:
				SCVF() = default;

			/// index of SubControlVolume on one side of the scvf
				[[nodiscard]] inline size_t from() const {return From;}

			/// index of SubControlVolume on one side of the scvf
				[[nodiscard]] inline size_t to() const {return To;}

			/// number of integration points on scvf
				[[nodiscard]] inline size_t num_ip() const {return nip;}

			///	integration weight
				[[nodiscard]] inline number weight(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vWeight[ip]*vDetJMap[ip];}

			/// local integration point of scvf
				inline const MathVector<dim>& local_ip(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vLocalIP[ip];}

			/// global integration point of scvf
				inline const MathVector<worldDim>& global_ip(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vGlobalIP[ip];}

			/// normal on scvf (points direction "from"->"to"), normalized
				inline const MathVector<worldDim>& normal() const {return Normal;}

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim,dim>& JTInv(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vJtInv[ip];}

			/// Determinant of Jacobian in integration point
				[[nodiscard]] inline number detJ(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vDetJ[ip];}

			/// number of shape functions
				[[nodiscard]] inline size_t num_sh() const {return nsh;}

			/// value of shape function i in integration point
				[[nodiscard]] inline number shape(size_t ip, size_t sh) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vvShape[ip][sh];}

			/// vector of shape functions in ip point
				inline const number* shape_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vvShape[ip];}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t ip, size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index");
					 UG_ASSERT(ip<num_ip(), "Wrong index"); return vvLocalGrad[ip][sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vvLocalGrad[ip];}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t ip, size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index");
					 UG_ASSERT(ip<num_ip(), "Wrong index"); return vvGlobalGrad[ip][sh];}

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vvGlobalGrad[ip];}

			/// number of corners, that bound the scvf
				[[nodiscard]] inline size_t num_corners() const {return numCo;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return global corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			private:
			// 	let outer class access private members
				friend class FVGeometry<TOrder, TElem, TWorldDim, TQuadOrder>;

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
				MathVector<dim> vLocalIP[nip]; // local integration point
				MathVector<worldDim> vGlobalIP[nip]; // global integration point
				const number* vWeight; // weight at ip

			// shapes and derivatives
				number vvShape[nip][nsh]; // shapes at ip
				MathVector<dim> vvLocalGrad[nip][nsh]; // local grad at ip
				MathVector<worldDim> vvGlobalGrad[nip][nsh]; // global grad at ip
				MathMatrix<worldDim,dim> vJtInv[nip]; // Jacobian transposed at ip
				number vDetJ[nip]; // Jacobian det at ip
				number vDetJMap[nip]; // Jacobian det at ip
		};

	///	sub control volume structure
		class SCV
		{
			public:
			///	Number of integration points
				static constexpr size_t nip = scv_quad_rule_type::nip;

			/// Number of corners of scvf
				static constexpr size_t numCo =	traits::NumCornersOfSCV;

			public:
				SCV() = default;

			/// number of corners, that bound the scvf
				[[nodiscard]] inline size_t num_corners() const {return numCo;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			/// node id that this scv is associated to
				[[nodiscard]] inline size_t node_id() const {return nodeId;}

			/// number of integration points
				[[nodiscard]] inline size_t num_ip() const {return nip;}

			///	weigth of integration point
				[[nodiscard]] inline number weight(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vWeight[ip]*vDetJMap[ip];}

			/// local integration point of scv
				inline const MathVector<dim>& local_ip(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vLocalIP[ip];}

			/// global integration point
				inline const MathVector<worldDim>& global_ip(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vGlobalIP[ip];}

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim,dim>& JTInv(size_t ip) const
					{UG_ASSERT(ip<num_ip(),"Wring index."); return vJtInv[ip];}

			/// Determinant of Jacobian in integration point
				[[nodiscard]] inline number detJ(size_t ip) const
					{UG_ASSERT(ip<num_ip(),"Wring index."); return vDetJ[ip];}

			/// number of shape functions
				[[nodiscard]] inline size_t num_sh() const {return nsh;}

			/// value of shape function i in integration point
				[[nodiscard]] inline number shape(size_t ip, size_t sh) const
					{UG_ASSERT(ip<num_ip(),"Wring index."); return vvShape[ip][sh];}

			/// vector of shape functions in ip point
				[[nodiscard]] inline const number* shape_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(),"Wring index."); return vvShape[ip];}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t ip, size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index");
					 UG_ASSERT(ip<num_ip(),"Wring index."); return vvLocalGrad[ip][sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(),"Wring index."); return vvLocalGrad[ip];}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t ip, size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index");
					 UG_ASSERT(ip<num_ip(),"Wring index."); return vvGlobalGrad[ip][sh];}

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(),"Wring index."); return vvGlobalGrad[ip];}

			private:
			// 	let outer class access private members
				friend class FVGeometry<TOrder, TElem, TWorldDim, TQuadOrder>;

			//  node id of associated node
				size_t nodeId;

			//  scv part
				MathVector<dim> vLocalIP[nip]; // local integration point
				MathVector<worldDim> vGlobalIP[nip]; // global intergration point
				const number* vWeight; // weight at ip

			//	local and global positions of this element
				MathVector<dim> vLocPos[numCo]; // local position of node
				MathVector<worldDim> vGloPos[numCo]; // global position of node
				MidID vMidID[numCo]; // dimension and id of object, that's midpoint bounds the scv

			// shapes and derivatives
				number vvShape[nip][nsh]; // shapes at ip
				MathVector<dim> vvLocalGrad[nip][nsh]; // local grad at ip
				MathVector<worldDim> vvGlobalGrad[nip][nsh]; // global grad at ip
				MathMatrix<worldDim,dim> vJtInv[nip]; // Jacobian transposed at ip
				number vDetJ[nip]; // Jacobian det at ip
				number vDetJMap[nip]; // Jacobian det at ip for scv integral map
		};

	///	boundary face
		class BF
		{
			public:
			/// number of integration points
				static constexpr size_t nip = scvf_quad_rule_type::nip;

			/// Number of corners of bf
				static constexpr size_t numCo = traits::NumCornersOfSCVF;

			public:
				BF() = default;

			/// index of SubControlVolume of the bf
				[[nodiscard]] inline size_t node_id() const {return nodeId;}

			/// outer normal on bf. Norm is equal to area
				inline const MathVector<worldDim>& normal() const {return Normal;} // includes area

			/// volume of bf
				[[nodiscard]] inline number volume() const {return Vol;}

			/// number of integration points on scvf
				[[nodiscard]] inline size_t num_ip() const {return nip;}

			///	integration weight
				[[nodiscard]] inline number weight(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vWeight[ip];}

			/// local integration point of scvf
				inline const MathVector<dim>& local_ip(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vLocalIP[ip];}

			/// global integration point of scvf
				inline const MathVector<worldDim>& global_ip(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vGlobalIP[ip];}

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim,dim>& JTInv(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vJtInv[ip];}

			/// Determinant of Jacobian in integration point
				[[nodiscard]] inline number detJ(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vDetJ[ip];}

			/// number of shape functions
				[[nodiscard]] inline size_t num_sh() const {return nsh;}

			/// value of shape function i in integration point
				[[nodiscard]] inline number shape(size_t ip, size_t sh) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vvShape[ip][sh];}

			/// vector of shape functions in ip point
				[[nodiscard]] inline const number* shape_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vvShape[ip];}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t ip, size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index");
					 UG_ASSERT(ip<num_ip(), "Wrong index"); return vvLocalGrad[ip][sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vvLocalGrad[ip];}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t ip, size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index");
					 UG_ASSERT(ip<num_ip(), "Wrong index"); return vvGlobalGrad[ip][sh];}

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vvGlobalGrad[ip];}

			/// number of corners, that bound the scvf
				[[nodiscard]] inline size_t num_corners() const {return numCo;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			private:
			/// let outer class access private members
				friend class FVGeometry<TOrder, TElem, TWorldDim, TQuadOrder>;

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
				MathVector<dim> vLocalIP[nip]; // local integration point
				MathVector<worldDim> vGlobalIP[nip]; // global integration point
				const number* vWeight; // weight at ip
				MathVector<worldDim> Normal; // normal (incl. area)
				number Vol; // volume of bf

			// 	shapes and derivatives
				number vvShape[nip][nsh]; // shapes at ip
				MathVector<dim> vvLocalGrad[nip][nsh]; // local grad at ip
				MathVector<worldDim> vvGlobalGrad[nip][nsh]; // global grad at ip
				MathMatrix<worldDim,dim> vJtInv[nip]; // Jacobian transposed at ip
				number vDetJ[nip]; // Jacobian det at ip
		};


	public:
	/// construct object and initialize local values and sizes
		FVGeometry();

	///	update local data
		void update_local_data();

	/// update Geometry for roid
		void update_local(ReferenceObjectID roid,
		                  const LFEID& lfeID = LFEID(LFEID::LAGRANGE, worldDim, 1),
		                  size_t orderQuad = TQuadOrder);

	/// update data for given element
		void update(GridObject* elem, const MathVector<worldDim>* vCornerCoords,
		            const ISubsetHandler* ish = nullptr);

	/// update boundary data for given element
		void update_boundary_faces(GridObject* elem,
		                           const MathVector<worldDim>* vCornerCoords,
		                           const ISubsetHandler* ish = nullptr);

	/// number of SubControlVolumeFaces
		[[nodiscard]] inline size_t num_scvf() const {return numSCVF;};

	/// const access to SubControlVolumeFace number i
		inline const SCVF& scvf(size_t i) const
			{UG_ASSERT(i < num_scvf(), "Invalid Index."); return m_vSCVF[i];}

	/// number of SubControlVolumes
		[[nodiscard]] inline size_t num_scv() const {return numSCV;}

	/// const access to SubControlVolume number i
		inline const SCV& scv(size_t i) const
			{UG_ASSERT(i < num_scv(), "Invalid Index."); return m_vSCV[i];}

	/// number of shape functions
		[[nodiscard]] inline size_t num_sh() const {return nsh;}

	public:
	/// returns number of all scvf ips
		[[nodiscard]] size_t num_scvf_ips() const {return numSCVFIP;}

	/// returns all ips of scvf as they appear in scv loop
		const MathVector<dim>* scvf_local_ips() const {return m_vLocSCVF_IP;}

	/// returns all ips of scvf as they appear in scv loop
		const MathVector<worldDim>* scvf_global_ips() const {return m_vGlobSCVF_IP;}

	/// returns number of all scv ips
		[[nodiscard]] size_t num_scv_ips() const {return numSCVIP;}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<dim>* scv_local_ips() const {return m_vLocSCV_IP;}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<worldDim>* scv_global_ips() const {return m_vGlobSCV_IP;}


	protected:
	//	global and local ips on SCVF
		MathVector<worldDim> m_vGlobSCVF_IP[numSCVFIP];
		MathVector<dim> m_vLocSCVF_IP[numSCVFIP];

	//	global and local ips on SCVF
		MathVector<worldDim> m_vGlobSCV_IP[numSCVIP];
		MathVector<dim> m_vLocSCV_IP[numSCVIP];

	protected:
	//	maximum number of geom objects in a dimension
		static constexpr int maxMid = numSCVF + 1;

	//	subelement
		struct SubElement
		{
		//	shape number of corners of subelement
			size_t vDoFID[ref_elem_type::numCorners];

		// 	local and global geom object midpoints for each dimension
		// 	(most objects in 1 dim, i.e. number of edges, but +1 for 1D)
			MathVector<dim> vvLocMid[dim+1][maxMid];
			MathVector<worldDim> vvGloMid[dim+1][maxMid];

		//	flag is subelement has boundary sides
			bool isBndElem;

		//	-1 is no bnd side, >= 0 corresponding side of whole element
			std::vector<int> vElemBndSide;
		};

	///	subelements
		SubElement m_vSubElem[numSubElem];

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
		[[nodiscard]] inline size_t num_bf() const
		{
			typename std::map<int, std::vector<BF> >::const_iterator it;
			size_t num = 0;
			for ( it=m_mapVectorBF.begin() ; it != m_mapVectorBF.end(); it++ )
				num += (*it).second.size();
			return num;
		}

	/// number of boundary faces on subset 'subsetIndex'
		[[nodiscard]] inline size_t num_bf(int si) const
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
			if(it == m_mapVectorBF.end()) UG_THROW("FVGeom: No BndSubset"<<si);
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

	///	corners of reference element
		MathVector<dim> m_vLocCorner[ref_elem_type::numCorners];

	///	SubControlVolumeFaces
		SCVF m_vSCVF[numSCVF];

	///	SubControlVolumes
		SCV m_vSCV[numSCV];

	///	Reference Mapping
		ReferenceMapping<ref_elem_type, worldDim> m_rMapping;

	///	Reference Element
		const ref_elem_type& m_rRefElem;

	///	Shape function set
		const local_shape_fct_set_type& m_rTrialSpace;

	///	Quad Rule scvf
		const scvf_quad_rule_type& m_rSCVFQuadRule;

	///	Quad Rule scv
		const scv_quad_rule_type& m_rSCVQuadRule;
};

////////////////////////////////////////////////////////////////////////////////
// Dimension FV Geometry (all orders)   DIM FV
////////////////////////////////////////////////////////////////////////////////

/// Geometry and shape functions for any order Vertex-Centered Finite Volume
/**
 * \tparam	TDim		reference element dim
 * \tparam	TWorldDim	(physical) world dimension
 */
template <int TWorldDim, int TDim = TWorldDim>
class DimFVGeometry : public FVGeometryBase
{
	public:
	///	traits used
	using traits = fv1_dim_traits<TDim, TWorldDim>;

	///	dimension of reference element
		static constexpr int dim = TDim;

	///	dimension of world
		static constexpr int worldDim = TWorldDim;

	public:
	///	type of SubControlVolumeFace
	using scvf_type = typename traits::scvf_type;

	///	type of SubControlVolume
	using scv_type = typename traits::scv_type;

	///	Hanging node flag: this Geometry does not support hanging nodes
		static constexpr bool usesHangingNodes = false;

	/// flag indicating if local data may change
		static constexpr bool staticLocalData = false;

	public:
	///	Sub-Control Volume Face structure
		class SCVF
		{
			public:
			///	Number of corners of scvf
				static constexpr size_t numCo = traits::NumCornersOfSCVF;

			public:
				SCVF() = default;

			/// index of SubControlVolume on one side of the scvf
				[[nodiscard]] inline size_t from() const {return From;}

			/// index of SubControlVolume on one side of the scvf
				[[nodiscard]] inline size_t to() const {return To;}

			/// number of integration points on scvf
				[[nodiscard]] inline size_t num_ip() const {return nip;}

			///	integration weight
				[[nodiscard]] inline number weight(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vWeight[ip]*vDetJMap[ip];}

			/// local integration point of scvf
				inline const MathVector<dim>& local_ip(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vLocalIP[ip];}

			/// global integration point of scvf
				inline const MathVector<worldDim>& global_ip(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vGlobalIP[ip];}

			/// normal on scvf (points direction "from"->"to"), normalized
				inline const MathVector<worldDim>& normal() const {return Normal;}

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim,dim>& JTInv(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vJtInv[ip];}

			/// Determinant of Jacobian in integration point
				[[nodiscard]] inline number detJ(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vDetJ[ip];}

			/// number of shape functions
				[[nodiscard]] inline size_t num_sh() const {return nsh;}

			/// value of shape function i in integration point
				[[nodiscard]] inline number shape(size_t ip, size_t sh) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vvShape[ip][sh];}

			/// vector of shape functions in ip point
				[[nodiscard]] inline const number* shape_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return &vvShape[ip][0];}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t ip, size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index");
					 UG_ASSERT(ip<num_ip(), "Wrong index"); return vvLocalGrad[ip][sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return &vvLocalGrad[ip][0];}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t ip, size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index");
					 UG_ASSERT(ip<num_ip(), "Wrong index"); return vvGlobalGrad[ip][sh];}

			/// vector of gloabl gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return &vvGlobalGrad[ip][0];}

			/// number of corners, that bound the scvf
				[[nodiscard]] inline size_t num_corners() const {return numCo;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			private:
			// 	let outer class access private members
				friend class DimFVGeometry<worldDim, dim>;

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
				size_t nip;
				std::vector<MathVector<dim> > vLocalIP; // local integration point (size: nip)
				std::vector<MathVector<worldDim> > vGlobalIP; // global integration point (size: nip)
				const number* vWeight; // weight at ip

			// shapes and derivatives
				size_t nsh;
				std::vector<std::vector<number> > vvShape; // shapes at ip (size: nip x nsh)
				std::vector<std::vector<MathVector<dim> > > vvLocalGrad; // local grad at ip (size: nip x nsh)
				std::vector<std::vector<MathVector<worldDim> > > vvGlobalGrad; // global grad at ip (size: nip x nsh)
				std::vector<MathMatrix<worldDim,dim> > vJtInv; // Jacobian transposed at ip (size: nip)
				std::vector<number> vDetJ; // Jacobian det at ip (size: nip)
				std::vector<number> vDetJMap; // Jacobian det at ip (size: nip)
		};

	///	sub control volume structure
		class SCV
		{
			public:
			/// Number of corners of scvf
				static constexpr size_t numCo =	traits::NumCornersOfSCV;

			public:
				SCV() = default;

			/// number of corners, that bound the scvf
				[[nodiscard]] inline size_t num_corners() const {return numCo;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			/// node id that this scv is associated to
				[[nodiscard]] inline size_t node_id() const {return nodeId;}

			/// number of integration points
				[[nodiscard]] inline size_t num_ip() const {return nip;}

			///	weigth of integration point
				[[nodiscard]] inline number weight(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vWeight[ip]*vDetJMap[ip];}

			/// local integration point of scv
				inline const MathVector<dim>& local_ip(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vLocalIP[ip];}

			/// global integration point
				inline const MathVector<worldDim>& global_ip(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vGlobalIP[ip];}

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim,dim>& JTInv(size_t ip) const
					{UG_ASSERT(ip<num_ip(),"Wring index."); return vJtInv[ip];}

			/// Determinant of Jacobian in integration point
				[[nodiscard]] inline number detJ(size_t ip) const
					{UG_ASSERT(ip<num_ip(),"Wring index."); return vDetJ[ip];}

			/// number of shape functions
				[[nodiscard]] inline size_t num_sh() const {return nsh;}

			/// value of shape function i in integration point
				[[nodiscard]] inline number shape(size_t ip, size_t sh) const
					{UG_ASSERT(ip<num_ip(),"Wring index."); return vvShape[ip][sh];}

			/// vector of shape functions in ip point
				[[nodiscard]] inline const number* shape_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(),"Wring index."); return &vvShape[ip][0];}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t ip, size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index");
					 UG_ASSERT(ip<num_ip(),"Wring index."); return vvLocalGrad[ip][sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(),"Wring index."); return &vvLocalGrad[ip][0];}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t ip, size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index");
					 UG_ASSERT(ip<num_ip(),"Wring index."); return vvGlobalGrad[ip][sh];}

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(),"Wring index."); return &vvGlobalGrad[ip][0];}

			private:
			// 	let outer class access private members
				friend class DimFVGeometry<worldDim, dim>;

			//  node id of associated node
				size_t nodeId;

			//	local and global positions of this element
				MathVector<dim> vLocPos[numCo]; // local position of node
				MathVector<worldDim> vGloPos[numCo]; // global position of node
				MidID vMidID[numCo]; // dimension and id of object, that's midpoint bounds the scv

			//  scv part
				size_t nip;
				std::vector<MathVector<dim> > vLocalIP; // local integration point (size: nip)
				std::vector<MathVector<worldDim> > vGlobalIP; // global intergration point (size: nip)
				const number* vWeight; // weight at ip

			// shapes and derivatives
				size_t nsh;
				std::vector<std::vector<number> > vvShape; // shapes at ip (size: nip x nsh)
				std::vector<std::vector<MathVector<dim> > > vvLocalGrad; // local grad at ip (size: nip x nsh)
				std::vector<std::vector<MathVector<worldDim> > > vvGlobalGrad; // global grad at ip (size: nip x nsh)
				std::vector<MathMatrix<worldDim,dim> > vJtInv; // Jacobian transposed at ip (size: nip)
				std::vector<number> vDetJ; // Jacobian det at ip (size: nip)
				std::vector<number> vDetJMap; // Jacobian det at ip (size: nip)
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
				[[nodiscard]] inline size_t node_id() const {return nodeId;}

			/// outer normal on bf. Norm is equal to area
				inline const MathVector<worldDim>& normal() const {return Normal;} // includes area

			/// volume of bf
				[[nodiscard]] inline number volume() const {return Vol;}

			/// number of integration points on scvf
				[[nodiscard]] inline size_t num_ip() const {return nip;}

			///	integration weight
				[[nodiscard]] inline number weight(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vWeight[ip];}

			/// local integration point of scvf
				inline const MathVector<dim>& local_ip(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vLocalIP[ip];}

			/// global integration point of scvf
				inline const MathVector<worldDim>& global_ip(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vGlobalIP[ip];}

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim,dim>& JTInv(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vJtInv[ip];}

			/// Determinant of Jacobian in integration point
				[[nodiscard]] inline number detJ(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vDetJ[ip];}

			/// number of shape functions
				[[nodiscard]] inline size_t num_sh() const {return nsh;}

			/// value of shape function i in integration point
				[[nodiscard]] inline number shape(size_t ip, size_t sh) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return vvShape[ip][sh];}

			/// vector of shape functions in ip point
				[[nodiscard]] inline const number* shape_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return &vvShape[ip][0];}

			/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t ip, size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index");
					 UG_ASSERT(ip<num_ip(), "Wrong index"); return vvLocalGrad[ip][sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return &vvLocalGrad[ip][0];}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t ip, size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index");
					 UG_ASSERT(ip<num_ip(), "Wrong index"); return vvGlobalGrad[ip][sh];}

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector(size_t ip) const
					{UG_ASSERT(ip<num_ip(), "Wrong index"); return &vvGlobalGrad[ip][0];}

			/// number of corners, that bound the scvf
				[[nodiscard]] inline size_t num_corners() const {return numCo;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			private:
			/// let outer class access private members
				friend class DimFVGeometry<worldDim, dim>;

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
				size_t nip;
				std::vector<MathVector<dim> > vLocalIP; // local integration point (size: nip)
				std::vector<MathVector<worldDim> > vGlobalIP; // global integration point (size: nip)
				const number* vWeight; // weight at ip
				MathVector<worldDim> Normal; // normal (incl. area)
				number Vol; // volume of bf

			// 	shapes and derivatives
				size_t nsh;
				std::vector<std::vector<number> > vvShape; // shapes at ip (size: nip x nsh)
				std::vector<std::vector<MathVector<dim> > > vvLocalGrad; // local grad at ip (size: nip x nsh)
				std::vector<std::vector<MathVector<worldDim> > > vvGlobalGrad; // global grad at ip (size: nip x nsh)
				std::vector<MathMatrix<worldDim,dim> > vJtInv; // Jacobian transposed at ip (size: nip)
				std::vector<number> vDetJ; // Jacobian det at ip (size: nip)
		};


	public:
	/// construct object and initialize local values and sizes
		DimFVGeometry();

	///	update local data
		void update_local(ReferenceObjectID roid, const LFEID& lfeID, size_t orderQuad);
		void update_local(ReferenceObjectID roid, const LFEID& lfeID)
		{
			update_local(roid, lfeID, lfeID.order() + 1);
		}

	/// update data for given element
		void update(GridObject* pElem, const MathVector<worldDim>* vCornerCoords,
		            const ISubsetHandler* ish = nullptr)
		{
			update(pElem, vCornerCoords, m_lfeID, m_quadOrderSCV, ish);
		}

	/// update data for given element
		void update(GridObject* pElem, const MathVector<worldDim>* vCornerCoords,
		            const LFEID& lfeID, size_t orderQuad,
		            const ISubsetHandler* ish = nullptr);

	/// update boundary data for given element
		void update_boundary_faces(GridObject* pElem,
		                           const MathVector<worldDim>* vCornerCoords,
		                           const ISubsetHandler* ish = nullptr);

	/// number of SubControlVolumeFaces
		[[nodiscard]] inline size_t num_scvf() const {return m_numSCVF;};

	/// const access to SubControlVolumeFace number i
		inline const SCVF& scvf(size_t i) const
			{UG_ASSERT(i < num_scvf(), "Invalid Index."); return m_vSCVF[i];}

	/// number of SubControlVolumes
		[[nodiscard]] inline size_t num_scv() const {return m_numSCV;}

	/// const access to SubControlVolume number i
		inline const SCV& scv(size_t i) const
			{UG_ASSERT(i < num_scv(), "Invalid Index."); return m_vSCV[i];}

	/// number of shape functions
		[[nodiscard]] inline size_t num_sh() const {return m_nsh;};

	public:
	/// returns number of all scvf ips
		[[nodiscard]] size_t num_scvf_ips() const {return m_numSCVFIP;}

	/// returns all ips of scvf as they appear in scv loop
		const MathVector<dim>* scvf_local_ips() const {return &(m_vLocSCVF_IP[0]);}

	/// returns all ips of scvf as they appear in scv loop
		const MathVector<worldDim>* scvf_global_ips() const {return &(m_vGlobSCVF_IP[0]);}

	/// returns number of all scv ips
		[[nodiscard]] size_t num_scv_ips() const {return m_numSCVIP;}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<dim>* scv_local_ips() const {return &(m_vLocSCV_IP[0]);}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<worldDim>* scv_global_ips() const {return &(m_vGlobSCV_IP[0]);}


	protected:
	//	global and local ips on SCVF (size: numSCVFIP)
		std::vector<MathVector<worldDim> > m_vGlobSCVF_IP;
		std::vector<MathVector<dim> > m_vLocSCVF_IP;

	//	global and local ips on SCVF
		std::vector<MathVector<worldDim> > m_vGlobSCV_IP;
		std::vector<MathVector<dim> > m_vLocSCV_IP;

	protected:
	//	miximum number of geom objects in a dimension
		static constexpr int maxMid = traits::maxNumSCVF +1;

	//	subelement
		struct SubElement
		{
		//	shape number of corners of subelement
			size_t vDoFID[traits::maxNumSCV];

		// 	local and global geom object midpoints for each dimension
		// 	(most objects in 1 dim, i.e. number of edges, but +1 for 1D)
			MathVector<dim> vvLocMid[dim+1][maxMid];
			MathVector<worldDim> vvGloMid[dim+1][maxMid];

		//	flag is subelement has boundary sides
			bool isBndElem;

		//	-1 is no bnd side, >= 0 corresponding side of whole element
			std::vector<int> vElemBndSide;
		};

	///	subelements (size: numSubElem)
		std::vector<SubElement> m_vSubElem;

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
		[[nodiscard]] inline size_t num_bf() const
		{
			typename std::map<int, std::vector<BF> >::const_iterator it;
			size_t num = 0;
			for ( it=m_mapVectorBF.begin() ; it != m_mapVectorBF.end(); ++it )
				num += (*it).second.size();
			return num;
		}

	/// number of boundary faces on subset 'subsetIndex'
		[[nodiscard]] inline size_t num_bf(int si) const
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
			if(it == m_mapVectorBF.end()) UG_THROW("DimFVGeom: No BndSubset: "<<si);
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
		GridObject* m_pElem;

	///	current reference object id
		ReferenceObjectID m_roid;

	///	current order
		int m_orderShape;

	///	current trial space
		LFEID m_lfeID;

	///	current number of subelements
		size_t m_numSubElem;

	///	number of shape functions
		size_t m_nsh;

	///	current number of scvf
		size_t m_numSCVF;

	///	current number of SCVF per SubElement
		size_t m_numSCVFPerSubElem;

	///	SubControlVolumeFaces (size: numSCVF)
		std::vector<SCVF> m_vSCVF;

	///	quadrature order
		int m_quadOrderSCVF;

	///	number of scvf ip
		size_t m_numSCVFIP;


	///	current number of scv
		size_t m_numSCV;

	///	number of SCV per SubElement
		size_t m_numSCVPerSubElem;

	///	SubControlVolumes (size: numSCV)
		std::vector<SCV> m_vSCV;

	///	current quadrature order scv
		int m_quadOrderSCV;

	///	number of scv ip
		size_t m_numSCVIP;
};

}

#endif