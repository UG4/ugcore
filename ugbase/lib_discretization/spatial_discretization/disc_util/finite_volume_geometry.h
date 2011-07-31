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
#include "lib_discretization/reference_element/reference_element.h"
#include "lib_discretization/reference_element/reference_element_traits.h"
#include "lib_discretization/local_finite_element/local_shape_function_set.h"
#include "lib_discretization/local_finite_element/lagrange/lagrangep1.h"
#include "finite_volume_util.h"

namespace ug{


///	a singleton class that returns a new id for each type
class UniqueFVGeomIDProvider{
	public:
		static UniqueFVGeomIDProvider& inst(){
			static UniqueFVGeomIDProvider instance;
			return instance;
		}

		size_t new_id()	{return ++m_id;}

	private:
		UniqueFVGeomIDProvider() : m_id(0)	{}
		size_t m_id;
};

///	This method associates a unique unsigned integer value with each type.
template <class TType>
size_t GetUniqueFVGeomID()
{
	static size_t typeID = UniqueFVGeomIDProvider::inst().new_id();
	return typeID;
}


/// base class for all FVGeometries
class FVGeometryBase {};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// FV1 Geometry for Reference Element Type
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// helper class to store dimnesion and id of a midpoint of a sub-element
struct MidID
{
		MidID() : dim(0), id(0) {};
		MidID(size_t dim_, size_t id_) : dim(dim_), id(id_) {};
		size_t dim;
		size_t id;
};

/// Geometry and shape functions for 1st order Vertex-Centered Finite Volume
/**
 * \tparam	TElem		Element type
 * \tparam	TWorldDim	(physical) world dimension
 */
template <	typename TElem, int TWorldDim>
class FV1Geometry : public FVGeometryBase
{
	public:
	// 	type of element
		typedef TElem elem_type;

	// 	type of reference element
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	public:
	// 	order
		static const int order = 1;

	// 	number of SubControlVolumes
		static const size_t numSCV = ref_elem_type::num_corners;

	// 	type of SubControlVolume
		typedef typename fv1_traits<ref_elem_type, TWorldDim>::scv_type scv_type;

	// 	number of SubControlVolumeFaces
		static const size_t numSCVF = ref_elem_type::num_edges;

	//	type of Shape function used
		typedef LagrangeP1<ref_elem_type, 1> local_shape_fct_set_type;

	//	number of shape functions
		static const size_t nsh = local_shape_fct_set_type::nsh;

	public:
	// 	dimension of reference element
		static const int dim = ref_elem_type::dim;

	// 	dimension of world
		static const int worldDim = TWorldDim;

	// 	Hanging node flag: this Geometry does not support hanging nodes
		static const bool usesHangingNodes = false;

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
			///	Number of integration points
				static const size_t numIP = 1;

			///	Number of corners of scvf
				static const size_t numCorners =
						fv1_traits<ref_elem_type, TWorldDim>::NumCornersOfSCVF;

			private:
			// 	let outer class access private members
				friend class FV1Geometry<TElem, TWorldDim>;

			public:
				SCVF() {}

			/// index of SubControlVolume on one side of the scvf
				inline size_t from() const {return m_from;}

			/// index of SubControlVolume on one side of the scvf
				inline size_t to() const {return m_to;}

			/// number of integration points on scvf
				inline size_t num_ip() const {return numIP;}

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
				inline size_t num_sh() const {return nsh;}

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

			/// vector of gloabl gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector() const {return globalGrad;}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return numCorners;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			private:
			// This scvf separates the scv with the ids given in "from" and "to"
			// The computed normal points in direction from->to
				size_t m_from, m_to;

			//	The normal on the SCVF pointing (from -> to)
				MathVector<worldDim> Normal; // normal (incl. area)

			// ordering is:
			// 1D: edgeMidPoint
			// 2D: edgeMidPoint, CenterOfElement
			// 3D: edgeMidPoint, Side one, CenterOfElement, Side two
				MathVector<dim> vLocPos[numCorners]; // local corners of scvf
				MathVector<worldDim> vGloPos[numCorners]; // global corners of scvf
				MidID midId[numCorners]; // dimension and id of object, that's midpoint bounds the scvf

			// scvf part
				MathVector<dim> localIP; // local integration point
				MathVector<worldDim> globalIP; // global intergration point

			// shapes and derivatives
				number vShape[nsh]; // shapes at ip
				MathVector<dim> localGrad[nsh]; // local grad at ip
				MathVector<worldDim> globalGrad[nsh]; // global grad at ip
				MathMatrix<worldDim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

	///	sub control volume structure
		class SCV
		{
			public:
			/// Number of integration points
				static const size_t numIP = 1;

			/// Number of corners of scvf
				static const size_t maxNumCorners =
						fv1_traits<ref_elem_type, TWorldDim>::MaxNumCornersOfSCV;

			private:
			// 	let outer class access private members
				friend class FV1Geometry<TElem, TWorldDim>;

			public:
				SCV() : m_numCorners(maxNumCorners) {};

			/// node id that this scv is associated to
				inline size_t node_id() const {return nodeId;}

			/// number of integration points
				inline size_t num_ip() const {return numIP;}

			/// local integration point of scv
				inline const MathVector<dim>& local_ip() const {return vLocPos[0];}

			/// global integration point
				inline const MathVector<worldDim>& global_ip() const {return vGloPos[0];}

			/// volume of scv
				inline number volume() const {return vol;}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return m_numCorners;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			private:
			//  node id of associated node
				size_t nodeId;

			//	volume of scv
				number vol;

			//	number of corners of this element
				size_t m_numCorners;

			//	local and global positions of this element
				MathVector<dim> vLocPos[maxNumCorners]; // local position of node
				MathVector<worldDim> vGloPos[maxNumCorners]; // global position of node
				MidID midId[maxNumCorners]; // dimension and id of object, that's midpoint bounds the scv
		};

	///	boundary face
		class BF
		{
			public:
			/// number of integration points
				static const size_t m_numIP = 1;

			/// Number of corners of bf
				static const size_t m_numCorners =
						fv1_traits<ref_elem_type, TWorldDim>::NumCornersOfSCVF;

			private:
			/// let outer class access private members
				friend class FV1Geometry<TElem, TWorldDim>;

			public:
				BF() {}

			/// index of SubControlVolume of the bf
				inline size_t node_id() const {return m_nodeId;}

			/// number of integration points on bf
				inline size_t num_ip() const {return m_numIP;}

			/// local integration point of bf
				inline const MathVector<dim>& local_ip() const {return localIP;}

			/// global integration point of bf
				inline const MathVector<worldDim>& global_ip() const {return globalIP;}

			/// outer normal on bf. Norm is equal to area
				inline const MathVector<worldDim>& normal() const {return Normal;} // includes area

			/// volume of bf
				inline number volume() const {return m_volume;}

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
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return localGrad[sh];}

			/// vector of local gradients in ip point
				inline const MathVector<dim>* local_grad_vector() const {return localGrad;}

			/// value of global gradient of shape function i in integration point
				inline const MathVector<worldDim>& global_grad(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return globalGrad[sh];}

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector() const {return globalGrad;}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return m_numCorners;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			private:
			// 	id of scv this bf belongs to
				size_t m_nodeId;

			// ordering is:
			// 1D: edgeMidPoint
			// 2D: edgeMidPoint, CenterOfElement
			// 3D: edgeMidPoint, Side one, CenterOfElement, Side two
				MathVector<dim> vLocPos[m_numCorners]; // local corners of bf
				MathVector<worldDim> vGloPos[m_numCorners]; // global corners of bf
				MidID midId[m_numCorners]; // dimension and id of object, that's midpoint bounds the bf

			// 	scvf part
				MathVector<dim> localIP; // local integration point
				MathVector<worldDim> globalIP; // global intergration point
				MathVector<worldDim> Normal; // normal (incl. area)
				number m_volume; // volume of bf

			// 	shapes and derivatives
				number vShape[nsh]; // shapes at ip
				MathVector<dim> localGrad[nsh]; // local grad at ip
				MathVector<worldDim> globalGrad[nsh]; // global grad at ip
				MathMatrix<worldDim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

	public:
	/// construct object and initialize local values and sizes
		FV1Geometry();

	///	update local data
		bool update_local_data();

	/// update data for given element
		bool update(TElem* elem, const MathVector<worldDim>* vCornerCoords,
		            const ISubsetHandler* ish = NULL);

	/// get vector of corners for current element
		const MathVector<worldDim>* corners() const {return m_gloMid[0];}

	/// number of SubControlVolumeFaces
		inline size_t num_scvf() const {return numSCVF;};

	/// const access to SubControlVolumeFace number i
		inline const SCVF& scvf(size_t i) const
			{UG_ASSERT(i < num_scvf(), "Invalid Index."); return m_vSCVF[i];}

	/// number of SubControlVolumes
		inline size_t num_scv() const {return numSCV;}

	/// const access to SubControlVolume number i
		inline const SCV& scv(size_t i) const
			{UG_ASSERT(i < num_scv(), "Invalid Index."); return m_vSCV[i];}

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
		const MathVector<dim>* scv_local_ips() const {return &(m_locMid[0][0]);}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<worldDim>* scv_global_ips() const {return &(m_gloMid[0][0]);}


	protected:
	//	global and local ips on SCVF
		MathVector<worldDim> m_vGlobSCVF_IP[numSCVF];
		MathVector<dim> m_vLocSCVF_IP[numSCVF];

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

	private:
	///	pointer to current element
		TElem* m_pElem;

	// 	local and global geom object midpoints for each dimension
	// 	(most objects in 1 dim, i.e. number of edges, but +1 for 1D)
		MathVector<dim> m_locMid[dim+1][numSCVF + 1];
		MathVector<worldDim> m_gloMid[dim+1][numSCVF +1];

	// 	SubControlVolumeFaces
		SCVF m_vSCVF[numSCVF];

	// 	SubControlVolumes
		SCV m_vSCV[numSCV];

	// 	Reference Mapping
		ReferenceMapping<ref_elem_type, worldDim> m_rMapping;

	// 	Reference Element
		const ref_elem_type& m_rRefElem;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Dim-dependent Finite Volume Geometry
////////////////////////////////////////////////////////////////////////////////
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
	// 	dimension of reference element
		static const int dim = TDim;

	// 	dimension of world
		static const int worldDim = TWorldDim;

	// 	Hanging node flag: this Geometry does not support hanging nodes
		static const bool usesHangingNodes = false;

	public:
	// 	order
		static const int order = 1;

	// 	number of SubControlVolumes
		static const size_t maxNumSCV = fv1_dim_traits<dim, worldDim>::maxNumSCV;

	// 	type of SubControlVolume
		typedef typename fv1_dim_traits<dim, worldDim>::scv_type scv_type;

	// 	number of SubControlVolumeFaces
		static const size_t maxNumSCVF = fv1_dim_traits<dim, worldDim>::maxNumSCVF;

	//	number of shape functions
		static const size_t maxNSH = fv1_dim_traits<dim, worldDim>::maxNSH;

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
			///	Number of integration points
				static const size_t numIP = 1;

			///	Number of corners of scvf
				static const size_t numCorners = fv1_dim_traits<dim, worldDim>::NumCornersOfSCVF;

			private:
			// 	let outer class access private members
				friend class DimFV1Geometry<dim, worldDim>;

			public:
				SCVF() {}

			/// index of SubControlVolume on one side of the scvf
				inline size_t from() const {return m_from;}

			/// index of SubControlVolume on one side of the scvf
				inline size_t to() const {return m_to;}

			/// number of integration points on scvf
				inline size_t num_ip() const {return numIP;}

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
				inline size_t num_sh() const {return m_numSH;}

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

			/// vector of gloabl gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector() const {return globalGrad;}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return numCorners;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			private:
			// This scvf separates the scv with the ids given in "from" and "to"
			// The computed normal points in direction from->to
				size_t m_from, m_to;

			//	The normal on the SCVF pointing (from -> to)
				MathVector<worldDim> Normal; // normal (incl. area)

			// ordering is:
			// 1D: edgeMidPoint
			// 2D: edgeMidPoint, CenterOfElement
			// 3D: edgeMidPoint, Side one, CenterOfElement, Side two
				MathVector<dim> vLocPos[numCorners]; // local corners of scvf
				MathVector<worldDim> vGloPos[numCorners]; // global corners of scvf
				MidID midId[numCorners]; // dimension and id of object, that's midpoint bounds the scvf

			// scvf part
				MathVector<dim> localIP; // local integration point
				MathVector<worldDim> globalIP; // global intergration point

			///	current number of shape functions
				size_t m_numSH;

			// shapes and derivatives
				number vShape[maxNSH]; // shapes at ip
				MathVector<dim> localGrad[maxNSH]; // local grad at ip
				MathVector<worldDim> globalGrad[maxNSH]; // global grad at ip
				MathMatrix<worldDim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

	///	sub control volume structure
		class SCV
		{
			public:
			/// Number of integration points
				static const size_t numIP = 1;

			/// Number of corners of scvf
				static const size_t maxNumCorners = fv1_dim_traits<dim, worldDim>::MaxNumCornersOfSCV;

			private:
			// 	let outer class access private members
				friend class DimFV1Geometry<dim, worldDim>;

			public:
				SCV() : m_numCorners(maxNumCorners) {};

			/// node id that this scv is associated to
				inline size_t node_id() const {return nodeId;}

			/// number of integration points
				inline size_t num_ip() const {return numIP;}

			/// local integration point of scv
				inline const MathVector<dim>& local_ip() const {return vLocPos[0];}

			/// global integration point
				inline const MathVector<worldDim>& global_ip() const {return vGloPos[0];}

			/// volume of scv
				inline number volume() const {return vol;}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return m_numCorners;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim,dim>& JTInv() const {return JtInv;}

			/// Determinant of Jacobian in integration point
				inline number detJ() const {return detj;}

			/// number of shape functions
				inline size_t num_sh() const {return m_numSH;}

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

			/// vector of gloabl gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector() const {return globalGrad;}

			private:
			//  node id of associated node
				size_t nodeId;

			//	volume of scv
				number vol;

			//	number of corners of this element
				size_t m_numCorners;

			///	current number of shape functions
				size_t m_numSH;

			//	local and global positions of this element
				MathVector<dim> vLocPos[maxNumCorners]; // local position of node
				MathVector<worldDim> vGloPos[maxNumCorners]; // global position of node
				MidID midId[maxNumCorners]; // dimension and id of object, that's midpoint bounds the scv

			// shapes and derivatives
				number vShape[maxNSH]; // shapes at ip
				MathVector<dim> localGrad[maxNSH]; // local grad at ip
				MathVector<worldDim> globalGrad[maxNSH]; // global grad at ip
				MathMatrix<worldDim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

	///	boundary face
		class BF
		{
			public:
			/// number of integration points
				static const size_t m_numIP = 1;

			/// Number of corners of bf
				static const size_t m_numCorners = fv1_dim_traits<dim, worldDim>::NumCornersOfSCVF;

			private:
			/// let outer class access private members
				friend class DimFV1Geometry<dim, worldDim>;

			public:
				BF() {}

			/// index of SubControlVolume of the bf
				inline size_t node_id() const {return m_nodeId;}

			/// number of integration points on bf
				inline size_t num_ip() const {return m_numIP;}

			/// local integration point of bf
				inline const MathVector<dim>& local_ip() const {return localIP;}

			/// global integration point of bf
				inline const MathVector<worldDim>& global_ip() const {return globalIP;}

			/// outer normal on bf. Norm is equal to area
				inline const MathVector<worldDim>& normal() const {return Normal;} // includes area

			/// volume of bf
				inline number volume() const {return m_volume;}

			/// Transposed Inverse of Jacobian in integration point
				inline const MathMatrix<worldDim, dim>& JTInv() const {return JtInv;}

			/// Determinant of Jacobian in integration point
				inline number detJ() const {return detj;}

			/// number of shape functions
				inline size_t num_sh() const {return m_numSH;}

			/// value of shape function i in integration point
				inline number shape(size_t sh) const
					{UG_ASSERT(sh < num_sh(), "Invalid index"); return vShape[sh];}

			/// vector of local gradients in ip point
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

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return m_numCorners;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			private:
			// 	id of scv this bf belongs to
				size_t m_nodeId;

			// ordering is:
			// 1D: edgeMidPoint
			// 2D: edgeMidPoint, CenterOfElement
			// 3D: edgeMidPoint, Side one, CenterOfElement, Side two
				MathVector<dim> vLocPos[m_numCorners]; // local corners of bf
				MathVector<worldDim> vGloPos[m_numCorners]; // global corners of bf
				MidID midId[m_numCorners]; // dimension and id of object, that's midpoint bounds the bf

			// 	scvf part
				MathVector<dim> localIP; // local integration point
				MathVector<worldDim> globalIP; // global intergration point
				MathVector<worldDim> Normal; // normal (incl. area)
				number m_volume; // volume of bf

			///	current number of shape functions
				size_t m_numSH;

			// 	shapes and derivatives
				number vShape[maxNSH]; // shapes at ip
				MathVector<dim> localGrad[maxNSH]; // local grad at ip
				MathVector<worldDim> globalGrad[maxNSH]; // global grad at ip
				MathMatrix<worldDim,dim> JtInv; // Jacobian transposed at ip
				number detj; // Jacobian det at ip
		};

	public:
	/// construct object and initialize local values and sizes
		DimFV1Geometry() : m_pElem(NULL), m_roid(ROID_INVALID) {};

	///	update local data
		bool update_local_data();

	/// update data for given element
		bool update(GeometricObject* elem, const MathVector<worldDim>* vCornerCoords,
		            const ISubsetHandler* ish = NULL);

	/// get vector of corners for current element
		const MathVector<worldDim>* corners() const {return m_gloMid[0];}

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
	/// returns number of all scvf ips
		size_t num_scvf_ips() const {return m_numSCVF;}

	/// returns all ips of scvf as they appear in scv loop
		const MathVector<dim>* scvf_local_ips() const {return m_vLocSCVF_IP;}

	/// returns all ips of scvf as they appear in scv loop
		const MathVector<worldDim>* scvf_global_ips() const {return m_vGlobSCVF_IP;}

	/// returns number of all scv ips
		size_t num_scv_ips() const {return m_numSCV;}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<dim>* scv_local_ips() const {return &(m_locMid[0][0]);}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<worldDim>* scv_global_ips() const {return &(m_gloMid[0][0]);}


	protected:
	//	global and local ips on SCVF
		MathVector<worldDim> m_vGlobSCVF_IP[maxNumSCVF];
		MathVector<dim> m_vLocSCVF_IP[maxNumSCVF];

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

	private:
	///	pointer to current element
		GeometricObject* m_pElem;

	///	current reference object id
		ReferenceObjectID m_roid;

	///	current number of scvf
		size_t m_numSCVF;

	///	current number of scv
		size_t m_numSCV;

	// 	local and global geom object midpoints for each dimension
	// 	(most objects in 1 dim, i.e. number of edges, but +1 for 1D)
		MathVector<dim> m_locMid[dim+1][maxNumSCVF + 1];
		MathVector<worldDim> m_gloMid[dim+1][maxNumSCVF +1];

	// 	SubControlVolumeFaces
		SCVF m_vSCVF[maxNumSCVF];

	// 	SubControlVolumes
		SCV m_vSCV[maxNumSCV];
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// FV1 Monifold Boundary
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <	typename TElem,
			int TWorldDim>
class FV1ManifoldBoundary
{
	public:
		// type of element
		typedef TElem elem_type;

		// type of reference element
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	public:
		// order
		static const int order = 1;

		// number of BoundaryFaces
		static const size_t m_numBF = ref_elem_type::num_corners;

		// type of BoundaryFaces
		typedef typename fv1_traits<ref_elem_type, TWorldDim>::scv_type bf_type;

	public:
		// dimension of reference element
		static const int dim = ref_elem_type::dim;

		// dimension of world
		static const int worldDim = TWorldDim;

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
		class BF
		{
			private:
				// let outer class access private members
				friend class FV1ManifoldBoundary<TElem, TWorldDim>;

				// number of integration points
				static const size_t m_numIP = 1;
				
				// max number of corners of bf
				static const size_t m_maxNumCorners = fv1_traits<ref_elem_type, TWorldDim>::MaxNumCornersOfSCV;

			public:
				BF() : m_numCorners(m_maxNumCorners) {};

				/// node id that this bf is associated to
				inline size_t node_id() const {return nodeId;}

				/// number of integration points
				inline size_t num_ip() const {return m_numIP;}

				/// local integration point of bf
				inline const MathVector<dim>& local_ip(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid index"); return vLocPos[0];}

				/// global integration point
				inline const MathVector<worldDim>& global_ip(size_t ip) const
					{UG_ASSERT(ip < num_ip(), "Invalid index"); return vGloPos[0];}

				/// volume of bf
				inline number volume() const {return vol;}

				/// number of corners, that bound the bf
				inline size_t num_corners() const {return m_numCorners;}

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
				size_t m_numCorners;
				MathVector<dim> vLocPos[m_maxNumCorners];			// local position of node
				MathVector<worldDim> vGloPos[m_maxNumCorners];	// global position of node
				MidID midId[m_maxNumCorners];			// dimension and id of object, whose midpoint bounds the scv
				
				//IPs & shapes
				MathVector<dim> localIP; // local integration point
				MathVector<worldDim> globalIP; // global intergration point
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
		FV1ManifoldBoundary();
		
	///	update data for given element
		bool update(TElem* elem, const MathVector<worldDim>* vCornerCoords,
		            const ISubsetHandler* ish = NULL);
			
	/// get vector of corners for current element
		const MathVector<worldDim>* corners() const {return m_gloMid[0];}

	/// number of BoundaryFaces
		inline size_t num_bf() const {return m_numBF;}

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

	
	private:
		// pointer to current element
		TElem* m_pElem;
		
		// local and global geom object midpoints for each dimension
		MathVector<dim> m_locMid[dim+1][m_numBF];
		MathVector<worldDim> m_gloMid[dim+1][m_numBF];
		
		// SubControlVolumes
		BF m_vBF[m_numBF];
		
		// Reference Mapping
		ReferenceMapping<ref_elem_type, worldDim> m_rMapping;

		// Reference Element
		const ref_elem_type& m_rRefElem;

		
};

}

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_GEOMETRY__ */
