/*
 * fvcr_geom.h
 *
 *  Created on: 21.06.2012
 *      Author: Christian Wehner
 *
 * Node centered finite volume geometry for Crouzeix-Raviart-Elements
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__CR_FINITE_VOLUME_GEOMETRY__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__CR_FINITE_VOLUME_GEOMETRY__

// extern libraries
#include <cmath>
#include <map>
#include <vector>

// other ug4 modules
#include "common/common.h"

// library intern includes
#include "lib_disc/reference_element/reference_element.h"
#include "lib_disc/reference_element/reference_mapping.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/reference_element/reference_element_traits.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "lib_disc/local_finite_element/crouzeix-raviart/crouzeix_raviart.h"
#include "lib_disc/quadrature/gauss/gauss_quad.h"
#include "lib_disc/common/geometry_util.h"
#include "lib_disc/domain_util_impl.h"

#include "fv_geom_base.h"
#include "fv_util.h"

namespace ug{

/* crfv traits */

template <int TWorldDim,int nrfaceco> struct crfv_traits
{
    typedef void scv_type;
    typedef void face_type;
	static const size_t maxNumSCVF;
	static const size_t maxNumSCV;
	static const size_t maxNSH;
	static const size_t maxNumCo;
};

template <> struct crfv_traits<1, 1>
{
    typedef ReferenceEdge scv_type;
    typedef ReferenceVertex face_type;
	static const size_t maxNumSCVF = 1;
	static const size_t maxNumSCV = 2;
	static const size_t maxNSH = maxNumSCV;
	static const size_t maxNumCo = 2;
};

template <> struct crfv_traits<1, 2>
{
	typedef ReferenceEdge scv_type;
    typedef ReferenceVertex face_type;
	static const size_t maxNumSCVF = 1;
	static const size_t maxNumSCV = 2;
	static const size_t maxNSH = maxNumSCV;
	static const size_t maxNumCo = 2;
};

template <> struct crfv_traits<2, 2>
{
    typedef ReferenceTriangle scv_type;
    typedef ReferenceEdge face_type;
	static const size_t maxNumSCVF = 4;
	static const size_t maxNumSCV = 4;
	static const size_t maxNSH = maxNumSCV;
	static const size_t maxNumCo = 4;
};

template <> struct crfv_traits<2, 3>
{
    typedef ReferenceTriangle scv_type;
    typedef ReferenceEdge face_type;
	static const size_t maxNumSCVF = 4;
	static const size_t maxNumSCV = 4;
	static const size_t maxNSH = maxNumSCV;
	static const size_t maxNumCo = 4;
};

template <> struct crfv_traits<3, 3>
{
	typedef ReferenceTetrahedron scv_type;
    typedef ReferenceTriangle face_type;
	static const size_t maxNumSCVF = 10;
	static const size_t maxNumSCV = 6;
	static const size_t maxNSH = maxNumSCV;
	static const size_t maxNumCo = 8;
};

template <> struct crfv_traits<3, 4>
{
    typedef ReferencePyramid scv_type;
    typedef ReferenceQuadrilateral face_type;
	static const size_t maxNumSCVF = 10;
	static const size_t maxNumSCV = 6;
	static const size_t maxNSH = maxNumSCV;
	static const size_t maxNumCo = 8;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Dim-dependent Finite Volume Geometry
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

///
/// Geometry and shape functions for 1st order Vertex-Centered Finite Volume with Crouzeix-Raviart elements
//
//  basics:
//  Unknowns are centered in element vertices(1d)/edges(2d)/faces(3d)
//  scv i corners are given by face i corners + element barycenter
//  scvf i corners are given by edge i + element barycenter
//
//  important differences to standard linear Lagrange FV class:
//
//  new scv class parameter "Normal" is normal on associated face  (this is computed for
//  Navier-Stokes-Crouzeix-Raviart-staggered grid continuity equation discretization)
//  midpoint computations and related parameters eliminated (not necessary for scv computation)
//  parameter nodeID refers to node located in edge/face for the unknown, not to corner node
//  ip of scv is not corner but associated edge/face midpoint
//  location of unknowns (edge/face barycenter) is stored in m_vLocUnkCoords, m_vGlobUnkCoords
//
/**
 * \tparam	TDim		reference element dim
 * \tparam	TWorldDim	(physical) world dimension
 */
template <int TDim, int TWorldDim = TDim>
class DimCRFVGeometry : public FVGeometryBase
{
	public:
	///	used traits
		typedef crfv_traits<TDim, TDim> traits;

	public:
	///	dimension of reference element
		static const int dim = TDim;

	///	dimension of world
		static const int worldDim = TWorldDim;

	///	Hanging node flag: this Geometry does not support hanging nodes
		static const bool usesHangingNodes = false;

	/// flag indicating if local data may change
		static const bool staticLocalData = false;

	public:
	///	order
		static const int order = 1;

	///	number of SubControlVolumes
		static const size_t maxNumSCV = traits::maxNumSCV;

	///	traits
		typedef typename crfv_traits<TWorldDim,TWorldDim>::scv_type scv_type0;
	    typedef typename crfv_traits<TWorldDim,TWorldDim>::face_type face_type0;
	    typedef typename crfv_traits<TWorldDim,TWorldDim+1>::scv_type scv_type1;
	    typedef typename crfv_traits<TWorldDim,TWorldDim+1>::face_type face_type1;

	///	max number of SubControlVolumeFaces
		static const size_t maxNumSCVF = traits::maxNumSCVF;

	/// max number of shape functions
		static const size_t maxNSH = traits::maxNSH;

		static const size_t maxNumCo = traits::maxNumCo;

	/// number of integration points
		static const size_t nip = 1;

	// local/global element barycenter
		MathVector<dim> localBary;
		MathVector<worldDim> globalBary;

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
				static const size_t numCo = dim;

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
				friend class DimCRFVGeometry<dim, worldDim>;

			// This scvf separates the scv with the ids given in "from" and "to"
			// The computed normal points in direction from->to
				size_t From, To;

			//	The normal on the SCVF pointing (from -> to)
				MathVector<worldDim> Normal; // normal (incl. area)

				MathVector<dim> vLocPos[numCo]; // local corners of scvf
				MathVector<worldDim> vGloPos[numCo]; // global corners of scvf

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
				static const size_t maxNumCo = 5;

			public:
				SCV() : numCorners(maxNumCo) {};

			/// volume of scv
				inline number volume() const {return Vol;}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return numCorners;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return global corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			/// node id that this scv is associated to
				inline size_t node_id() const {return nodeID;}
				
			/// number of integration points
				inline size_t num_ip() const {return nip;}

			/// local integration point of scv
				inline const MathVector<dim>& local_ip() const {return vLocIP;}

			/// global integration point
				inline const MathVector<worldDim>& global_ip() const {return vGlobIP;}
				
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

			/// vector of global gradients in ip point
				inline const MathVector<worldDim>* global_grad_vector() const {return vGlobalGrad;}

			private:
			// 	let outer class access private members
				friend class DimCRFVGeometry<dim, worldDim>;
				
			//	The normal on the associated face pointing outward
				MathVector<worldDim> Normal; // normal (incl. area)

			//  node id of associated node
				size_t nodeID;

			//	volume of scv
				number Vol;

			//	number of corners of this element
				int numCorners;

			//	local and global positions of this element
				MathVector<dim> vLocPos[maxNumCo]; // local position of node
				MathVector<worldDim> vGloPos[maxNumCo]; // global position of node
				MathVector<dim> vLocIP;
				MathVector<dim> vGlobIP;

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
			/// max number of corners of bf
				static const size_t maxNumCo=4;

			public:
				BF() {}

			/// index of SubControlVolume of the bf
				inline size_t node_id() const {return nodeID;}

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
				friend class DimCRFVGeometry<dim, worldDim>;

			// 	id of scv this bf belongs to
				size_t nodeID;

				size_t numCo;

				MathVector<dim> vLocPos[maxNumCo]; // local corners of bf
				MathVector<worldDim> vGloPos[maxNumCo]; // global corners of bf

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
		DimCRFVGeometry() : m_pElem(NULL), m_roid(ROID_UNKNOWN) {};

	///	update local data
		void update_local_data();

	/// update data for given element
		void update(GeometricObject* elem, const MathVector<worldDim>* vCornerCoords,
		            const ISubsetHandler* ish = NULL);

	/// update boundary data for given element
		void update_boundary_faces(GeometricObject* elem,
		                           const MathVector<worldDim>* vCornerCoords,
		                           const ISubsetHandler* ish = NULL);

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

	/// number of shape functions
		inline size_t num_sh() const {return m_nsh;};

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
		const MathVector<dim>* scv_local_ips() const {return m_vLocUnkCoords;}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<worldDim>* scv_global_ips() const {return m_vGlobUnkCoords;}
		
	/// returns local barycenter
		const MathVector<dim> local_bary() const {return localBary;}
		
    /// returns global barycenter
		const MathVector<worldDim> global_bary() const {return globalBary;}

	protected:
	//	global and local ips on SCVF
		MathVector<worldDim> m_vGlobSCVF_IP[maxNumSCVF];
		MathVector<dim> m_vLocSCVF_IP[maxNumSCVF];
	// coord of location for unknowns in faces (edge/face barycenter)
		MathVector<worldDim> m_vGlobUnkCoords[maxNumSCV];
		MathVector<dim> m_vLocUnkCoords[maxNumSCV];

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
			if(it == m_mapVectorBF.end()) return 0;
			else return (*it).second.size();
		}

	/// returns the boundary face i for subsetIndex
		inline const BF& bf(int subsetIndex, size_t i) const
		{
			UG_ASSERT(i < num_bf(subsetIndex), "Invalid index.");
			typename std::map<int, std::vector<BF> >::const_iterator it;
			it = m_mapVectorBF.find(subsetIndex);
			UG_ASSERT(it != m_mapVectorBF.end(), "Bnd Subset Index not requested.")
			return (*it).second[i];
		}

	/// returns reference to vector of boundary faces for subsetIndex
		inline const std::vector<BF>& bf(int subsetIndex) const
		{
			typename std::map<int, std::vector<BF> >::const_iterator it;
			it = m_mapVectorBF.find(subsetIndex);
			UG_ASSERT(it != m_mapVectorBF.end(), "Bnd Subset Index not requested.")
			return (*it).second;
		}

	protected:
		std::map<int, std::vector<BF> > m_mapVectorBF;

	private:
	///	pointer to current element
		GeometricObject* m_pElem;

	///	current reference object id
		ReferenceObjectID m_roid;

	///	current number of scv
		size_t m_numSCV;

	///	current number of scvf
		size_t m_numSCVF;

	/// current number of shape functions
		size_t m_nsh;

		MathVector<dim> m_ipCoord[maxNumSCVF];
		MathVector<dim> faceBaryCoord[maxNumSCV];
		MathVector<dim> cornerCoord[maxNumCo];

	///	SubControlVolumeFaces
		SCVF m_vSCVF[maxNumSCVF];

	///	SubControlVolumes
		SCV m_vSCV[maxNumSCV];
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Dimension-indipendent Finite Volume Geometry
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <	typename TElem, int TWorldDim>
class CRFVGeometry : public FVGeometryBase
{
	public:
	///	type of element
		typedef TElem elem_type;

	///	type of reference element
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	///	dimension of reference element
		static const int dim = ref_elem_type::dim;

	///	dimension of world
		static const int worldDim = TWorldDim;

	///	Hanging node flag: this Geometry does not support hanging nodes
		static const bool usesHangingNodes = false;

	/// flag indicating if local data may change
		static const bool staticLocalData = true;

	///	type of Shape function used
		typedef CrouzeixRaviartLSFS<ref_elem_type> local_shape_fct_set_type;

	///	number of shape functions
		static const size_t nsh = local_shape_fct_set_type::nsh;

	///	number of SubControlVolumes
		static const size_t numSCV = nsh;

	///	number of SubControlVolumeFaces
		static const size_t numSCVF = ref_elem_type::numEdges;

	public:
	///	order
		static const int order = 1;

	///	traits
		typedef typename crfv_traits<TWorldDim,TWorldDim>::scv_type scv_type0;
	    typedef typename crfv_traits<TWorldDim,TWorldDim>::face_type face_type0;
	    typedef typename crfv_traits<TWorldDim,TWorldDim+1>::scv_type scv_type1;
	    typedef typename crfv_traits<TWorldDim,TWorldDim+1>::face_type face_type1;

	/// number of integration points
		static const size_t nip = 1;

	// local/global element barycenter
		MathVector<dim> localBary;
		MathVector<worldDim> globalBary;

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
				static const size_t numCo = dim;

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

			/// return glbal corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			private:
			// 	let outer class access private members
				friend class CRFVGeometry<TElem, TWorldDim>;

			// This scvf separates the scv with the ids given in "from" and "to"
			// The computed normal points in direction from->to
				size_t From, To;

			//	The normal on the SCVF pointing (from -> to)
				MathVector<worldDim> Normal; // normal (incl. area)

				MathVector<dim> vLocPos[numCo]; // local corners of scvf
				MathVector<worldDim> vGloPos[numCo]; // global corners of scvf

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
			/// Number of corners of scv
				static const size_t maxNumCo = 5;

			public:
				SCV() : numCorners(maxNumCo) {};

			/// volume of scv
				inline number volume() const {return Vol;}

			/// number of corners, that bound the scvf
				inline size_t num_corners() const {return numCorners;}

			/// return local corner number i
				inline const MathVector<dim>& local_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vLocPos[co];}

			/// return global corner number i
				inline const MathVector<worldDim>& global_corner(size_t co) const
					{UG_ASSERT(co < num_corners(), "Invalid index."); return vGloPos[co];}

			/// node id that this scv is associated to
				inline size_t node_id() const {return nodeID;}

			/// number of integration points
				inline size_t num_ip() const {return nip;}

			/// local integration point of scv
				inline const MathVector<dim>& local_ip() const {return vLocIP;}

			/// global integration point
				inline const MathVector<worldDim>& global_ip() const {return vGlobIP;}

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
				friend class CRFVGeometry<TElem, TWorldDim>;

			//	The normal on the associated face pointing outward
				MathVector<worldDim> Normal; // normal (incl. area)

			//  node id of associated node
				size_t nodeID;

			//	volume of scv
				number Vol;

			//	number of corners of this element
				int numCorners;

			//	local and global positions of this element
				MathVector<dim> vLocPos[maxNumCo]; // local position of node
				MathVector<worldDim> vGloPos[maxNumCo]; // global position of node

				MathVector<dim> vLocIP; // local position of node
				MathVector<worldDim> vGlobIP; // global position of node

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
			/// max number of corners of bf
				static const size_t maxNumCo=4;

			public:
				BF() {}

			/// index of SubControlVolume of the bf
				inline size_t node_id() const {return nodeID;}

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
				friend class CRFVGeometry<TElem, TWorldDim>;

			// 	id of scv this bf belongs to
				size_t nodeID;

				size_t numCo;

				MathVector<dim> vLocPos[maxNumCo]; // local corners of bf
				MathVector<worldDim> vGloPos[maxNumCo]; // global corners of bf

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
		CRFVGeometry();

	///	update local data
		void update_local_data();

	/// update data for given element
		void update(GeometricObject* elem, const MathVector<worldDim>* vCornerCoords,
		            const ISubsetHandler* ish = NULL);

	/// update boundary data for given element
		void update_boundary_faces(GeometricObject* elem,
		                           const MathVector<worldDim>* vCornerCoords,
		                           const ISubsetHandler* ish = NULL);
								   
		const MathVector<worldDim>* corners() const {return m_vCo;}							

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

	/// number of shape functions
		inline size_t num_sh() const {return nsh;};

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
		const MathVector<dim>* scv_local_ips() const {return m_vLocUnkCoords;}

	/// returns all ips of scv as they appear in scv loop
		const MathVector<worldDim>* scv_global_ips() const {return m_vGlobUnkCoords;}

	/// returns local barycenter
		const MathVector<dim> local_bary() const {return localBary;}

	/// returns global barycenter
		const MathVector<worldDim> global_bary() const {return globalBary;}

	protected:
	//	global and local ips on SCVF
		MathVector<worldDim> m_vGlobSCVF_IP[numSCVF];
		MathVector<dim> m_vLocSCVF_IP[numSCVF];
	 // coord of location for unknowns in faces (edge/face barycenter)
		MathVector<worldDim> m_vGlobUnkCoords[numSCV];
		MathVector<dim> m_vLocUnkCoords[numSCV];
		
		static const size_t numMaxCo = 8;
	 // corner coordinates
		MathVector<worldDim> m_vCo[numMaxCo];

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
			if(it == m_mapVectorBF.end()) return 0;
			else return (*it).second.size();
		}

	/// returns the boundary face i for subsetIndex
		inline const BF& bf(int subsetIndex, size_t i) const
		{
			UG_ASSERT(i < num_bf(subsetIndex), "Invalid index.");
			typename std::map<int, std::vector<BF> >::const_iterator it;
			it = m_mapVectorBF.find(subsetIndex);
			UG_ASSERT(it != m_mapVectorBF.end(), "Bnd Subset Index not requested.")
			return (*it).second[i];
		}

	/// returns reference to vector of boundary faces for subsetIndex
		inline const std::vector<BF>& bf(int subsetIndex) const
		{
			typename std::map<int, std::vector<BF> >::const_iterator it;
			it = m_mapVectorBF.find(subsetIndex);
			UG_ASSERT(it != m_mapVectorBF.end(), "Bnd Subset Index not requested.")
			return (*it).second;
		}

	protected:
		std::map<int, std::vector<BF> > m_mapVectorBF;

	private:
		MathVector<dim> m_ipCoord[numSCVF];
		MathVector<dim> faceBaryCoord[numSCV];

		///	SubControlVolumeFaces
		SCVF m_vSCVF[numSCVF];

		///	SubControlVolumes
		SCV m_vSCV[numSCV];

		///	pointer to current element
		TElem* m_pElem;

		///	Reference Mapping
		ReferenceMapping<ref_elem_type, worldDim> m_mapping;

		///	Reference Element
		const ref_elem_type& m_rRefElem;

		///	Shape function set
		const local_shape_fct_set_type& m_rTrialSpace;
};


}

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__CR_FINITE_VOLUME_GEOMETRY__ */
