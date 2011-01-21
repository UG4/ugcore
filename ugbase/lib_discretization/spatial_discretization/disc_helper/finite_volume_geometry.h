/*
 * finite_volume_geometry.h
 *
 *  Created on: 04.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_GEOMETRY__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_GEOMETRY__

// extern libraries
#include <cmath>
#include <map>
#include <vector>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"

// library intern includes
#include "lib_discretization/common/common.h"
#include "../../reference_element/reference_element.h"
#include "../../local_shape_function_set/local_shape_function_set_provider.h"
#include "./finite_volume_util.h"

namespace ug{


template <	typename TElem,
			int TWorldDim>
class FV1Geometry {
	public:
		// type of element
		typedef TElem elem_type;

		// type of reference element
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	public:
		// order
		static const int order = 1;

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

		// Hanging node flag: this Geometry does not support hanging nodes
		static const bool usesHangingNodes = false;

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
				// type of element
				typedef TElem elem_type;

				// type of reference element
				typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

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
				inline const MathMatrix<world_dim,dim>& JTInv(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return JtInv;}

				/// Determinante of Jacobian in integration point
				inline number detJ(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return detj;}

				/// number of shape functions
				inline size_t num_sh() const {return vShape.size();}

				/// value of shape function i in integration point
				inline number shape(size_t i, size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return vShape[i];}

				/// vector of local gradients in ip point
				inline const std::vector<number>& shape_vector(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return vShape;}

				/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t i, size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return localGrad[i];}

				/// vector of local gradients in ip point
				inline const std::vector<MathVector<dim> >& local_grad_vector(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return localGrad;}

				/// value of global gradient of shape function i in integration point
				inline const MathVector<world_dim>& global_grad(size_t i, size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return globalGrad[i];}

				/// vector of gloabl gradients in ip point
				inline const std::vector<MathVector<world_dim> >& global_grad_vector(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return globalGrad;}

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

		class BF
		{
			public:
				// type of element
				typedef TElem elem_type;

				// type of reference element
				typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

			private:
				// let outer class access private members
				friend class FV1Geometry<TElem, TWorldDim>;

				// number of integration points
				static const size_t m_numIP = 1;

				// Number of corners of bf
				static const size_t m_numCorners = finite_volume_traits<ref_elem_type, TWorldDim>::NumCornersOfSCVF;

			public:
				BF() {}

				/// index of SubControlVolume of the bf
				inline size_t node_id() const {return m_nodeId;}

				/// number of integration points on bf
				inline size_t num_ip() const {return m_numIP;}

				/// local integration point of bf
				inline const MathVector<dim>& local_ip(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return localIP;}

				/// global integration point of bf
				inline const MathVector<world_dim>& global_ip(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return globalIP;}

				/// outer normal on bf. Norm is equal to area
				inline const MathVector<world_dim>& normal() const {return Normal;} // includes area

				/// volume of bf
				inline number volume() const {return m_volume;}

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

				/// vector of local gradients in ip point
				inline const std::vector<number>& shape_vector(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return vShape;}

				/// value of local gradient of shape function i in integration point
				inline const MathVector<dim>& local_grad(size_t i, size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return localGrad[i];}

				/// vector of local gradients in ip point
				inline const std::vector<MathVector<dim> >& local_grad_vector(size_t ip) const
						{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return localGrad;}

				/// value of global gradient of shape function i in integration point
				inline const MathVector<world_dim>& global_grad(size_t i, size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return globalGrad[i];}

				/// vector of global gradients in ip point
				inline const std::vector<MathVector<world_dim> >& global_grad_vector(size_t ip) const
						{UG_ASSERT(ip < num_ip(), "Invalid ip index"); return globalGrad;}

				/// number of corners, that bound the scvf
				inline size_t num_corners() const {return m_numCorners;}

				/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t i) const
					{UG_ASSERT(i < num_corners(), "Invalid corner index."); return m_vLocPos[i];}

				/// return glbal corner number i
				inline const MathVector<world_dim>& global_corner(size_t i) const
					{UG_ASSERT(i < num_corners(), "Invalid corner index."); return m_vGloPos[i];}

			private:
				// id of scv this bf belongs to
				size_t m_nodeId;

				// ordering is:
				// 1D: edgeMidPoint
				// 2D: edgeMidPoint, CenterOfElement
				// 3D: edgeMidPoint, Side one, CenterOfElement, Side two
				MathVector<dim> m_vLocPos[m_numCorners]; // local corners of bf
				MathVector<world_dim> m_vGloPos[m_numCorners]; // global corners of bf
				MidID midId[m_numCorners]; // dimension and id of object, that's midpoint bounds the bf

				// scvf part
				MathVector<dim> localIP; // local integration point
				MathVector<world_dim> globalIP; // global intergration point
				MathVector<world_dim> Normal; // normal (incl. area)
				number m_volume; // volume of bf

				// shapes and derivatives
				std::vector<number> vShape; // shapes at ip
				std::vector<MathVector<dim> > localGrad; // local grad at ip
				std::vector<MathVector<world_dim> > globalGrad; // global grad at ip
				MathMatrix<world_dim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

	public:
		/// construct object and initialize local values and sizes
		FV1Geometry();

		/// update data for given element
		bool update(TElem* elem, const ISubsetHandler& ish, const MathVector<world_dim>* vCornerCoords);

		/// get vector of corners for current element
		const MathVector<world_dim>* corners() const {return m_gloMid[0];}

		/// number of SubControlVolumeFaces
		inline size_t num_scvf() const {return m_numSCVF;};

		/// const access to SubControlVolumeFace number i
		inline const SCVF& scvf(size_t i) const
			{UG_ASSERT(i < num_scvf(), "Invalid Index."); return m_vSCVF[i];}

		/// number of SubControlVolumes
		inline size_t num_scv() const {return m_numSCV;}

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
		inline size_t num_bf(int subsetIndex) const
		{
			typename std::map<int, std::vector<BF> >::const_iterator it;
			it = m_mapVectorBF.find(subsetIndex);
			return (*it).second.size();
		}

		/// returns the boundary face i for subsetIndex
		inline const BF& bf(int subsetIndex, size_t i) const
		{
			UG_ASSERT(i < num_bf(subsetIndex), "Invalid index.");
			typename std::map<int, std::vector<BF> >::const_iterator it;
			it = m_mapVectorBF.find(subsetIndex);
			return (*it).second[i];
		}

		/// returns reference to vector of boundary faces for subsetIndex
		inline const std::vector<BF>& bf(int subsetIndex) const
		{
			typename std::map<int, std::vector<BF> >::const_iterator it;
			it = m_mapVectorBF.find(subsetIndex);
			return (*it).second;
		}

	protected:
		std::map<int, std::vector<BF> > m_mapVectorBF;

	protected:
		void copy_local_corners(SCVF& scvf)
		{
			for(size_t i = 0; i < scvf.num_corners(); ++i)
			{
				const size_t dim = scvf.midId[i].dim;
				const size_t id = scvf.midId[i].id;
				scvf.m_vLocPos[i] = m_locMid[dim][id];
			}
		}

		void copy_global_corners(SCVF& scvf)
		{
			for(size_t i = 0; i < scvf.num_corners(); ++i)
			{
				const size_t dim = scvf.midId[i].dim;
				const size_t id = scvf.midId[i].id;
				scvf.m_vGloPos[i] = m_gloMid[dim][id];
			}
		}

		void copy_local_corners(SCV& scv)
		{
			for(size_t i = 0; i < scv.num_corners(); ++i)
			{
				const size_t dim = scv.midId[i].dim;
				const size_t id = scv.midId[i].id;
				scv.m_vLocPos[i] = m_locMid[dim][id];
			}
		}

		void copy_global_corners(SCV& scv)
		{
			for(size_t i = 0; i < scv.num_corners(); ++i)
			{
				const size_t dim = scv.midId[i].dim;
				const size_t id = scv.midId[i].id;
				scv.m_vGloPos[i] = m_gloMid[dim][id];
			}
		}

		void copy_local_corners(BF& bf)
		{
			for(size_t i = 0; i < bf.num_corners(); ++i)
			{
				const size_t dim = bf.midId[i].dim;
				const size_t id = bf.midId[i].id;
				bf.m_vLocPos[i] = m_locMid[dim][id];
			}
		}

		void copy_global_corners(BF& bf)
		{
			for(size_t i = 0; i < bf.num_corners(); ++i)
			{
				const size_t dim = bf.midId[i].dim;
				const size_t id = bf.midId[i].id;
				bf.m_vGloPos[i] = m_gloMid[dim][id];
			}
		}

	private:
		// pointer to current element
		TElem* m_pElem;

		// local and global geom object midpoints for each dimension
		// (most objects in 1 dim, i.e. number of edges)
		MathVector<dim> m_locMid[dim+1][m_numSCVF];
		MathVector<world_dim> m_gloMid[dim+1][m_numSCVF];

		// SubControlVolumeFaces
		SCVF m_vSCVF[m_numSCVF];

		// SubControlVolumes
		SCV m_vSCV[m_numSCV];

		// Reference Mapping
		ReferenceMapping<ref_elem_type, world_dim> m_rMapping;

		// Reference Element
		ref_elem_type m_rRefElem;
};

}

#include "./finite_volume_geometry_impl.h"

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_GEOMETRY__ */
