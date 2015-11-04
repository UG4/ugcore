/*
 * hfvcr_geom.h
 *
 *  Created on: 17.09.2013
 *      Author: Christian Wehner
 *
 * Node centered finite volume geometry for Crouzeix-Raviart-Elements
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__HANGING_CR_FINITE_VOLUME_GEOMETRY__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__HANGING_CR_FINITE_VOLUME_GEOMETRY__

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
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/local_finite_element/crouzeix-raviart/crouzeix_raviart.h"
#include "lib_disc/quadrature/gauss/gauss_quad.h"
#include "lib_disc/common/geometry_util.h"
#include "lib_disc/domain_util_impl.h"

#include "fv_geom_base.h"
#include "fv_util.h"

namespace ug{

/* hcrfv traits */

template <int TWorldDim,int nrfaceco> struct hcrfv_traits
{
    typedef void scv_type;
    typedef void face_type;
	static const size_t maxNumSCVF;
	static const size_t maxNumSCV;
	static const size_t maxNSH;
	static const size_t maxNumCo;
};

template <> struct hcrfv_traits<2, 2>
{
    typedef ReferenceTriangle scv_type;
    typedef ReferenceEdge face_type;
	static const size_t maxNumSCVF = 8;
	static const size_t maxNumSCV = 8;
	static const size_t maxNSH = maxNumSCV;
	static const size_t maxNumCo = 4;
};

template <> struct hcrfv_traits<2, 3>
{
    typedef ReferenceTriangle scv_type;
    typedef ReferenceEdge face_type;
	static const size_t maxNumSCVF = 8;
	static const size_t maxNumSCV = 8;
	static const size_t maxNSH = maxNumSCV;
	static const size_t maxNumCo = 4;
};

template <> struct hcrfv_traits<3, 3>
{
	typedef ReferenceTetrahedron scv_type;
    typedef ReferenceTriangle face_type;
	static const size_t maxNumSCVF = 40;
	static const size_t maxNumSCV = 24;
	static const size_t maxNSH = maxNumSCV;
	static const size_t maxNumCo = 8;
};

template <> struct hcrfv_traits<3, 4>
{
    typedef ReferencePyramid scv_type;
    typedef ReferenceQuadrilateral face_type;
	static const size_t maxNumSCVF = 40;
	static const size_t maxNumSCV = 24;
	static const size_t maxNSH = maxNumSCV;
	static const size_t maxNumCo = 8;
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Dimension-indipendent Finite Volume Geometry
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <	typename TElem, int TWorldDim>
class HCRFVGeometry : public FVGeometryBase
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

	///	Hanging node flag: this geometry supports hanging nodes
		static const bool usesHangingNodes = true;

	/// flag indicating if local data may change
		static const bool staticLocalData = true;

	///	type of Shape function used
		typedef CrouzeixRaviartLSFS<ref_elem_type> local_shape_fct_set_type;

	///	number of shape functions
		static const size_t nsh = local_shape_fct_set_type::nsh;

	///	number of SubControlVolumes
		static const size_t numNaturalSCV = nsh;

	///	number of SubControlVolumeFaces
		static const size_t numNaturalSCVF = ref_elem_type::numEdges;
		
	///	used traits
		typedef hcrfv_traits<dim,worldDim> traits;

		static const size_t maxNumSCV = traits::maxNumSCV;
		
		static const size_t maxNumSCVF = traits::maxNumSCVF;

	public:
	///	order
		static const int order = 1;

	///	traits
		typedef typename hcrfv_traits<TWorldDim,TWorldDim>::scv_type scv_type0;
	    typedef typename hcrfv_traits<TWorldDim,TWorldDim>::face_type face_type0;
	    typedef typename hcrfv_traits<TWorldDim,TWorldDim+1>::scv_type scv_type1;
	    typedef typename hcrfv_traits<TWorldDim,TWorldDim+1>::face_type face_type1;

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
				friend class HCRFVGeometry<TElem, TWorldDim>;

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
				friend class HCRFVGeometry<TElem, TWorldDim>;

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

		class CONSTRAINED_DOF
		{
			public:
				static const size_t maxNumConstrainingDofs = 4;
				inline size_t constraining_dofs_index(size_t i) const{
					return cDofInd[i];
				}
				inline number constraining_dofs_weight(size_t i) const{
					return cDofWeights[i];
				}
				inline size_t index() const{
					return i;
				}
				inline size_t num_constraining_dofs() const{
					return numConstrainingDofs;
				}
			private:
				// 	let outer class access private members
				friend class HCRFVGeometry<TElem, TWorldDim>;

				//  constraining dofs indices
				size_t cDofInd[maxNumConstrainingDofs];
				//  weights
				number cDofWeights[maxNumConstrainingDofs];
				//  local index of dof in element
				size_t i;
				//  nr of constraining dofs
				size_t numConstrainingDofs;
		};
	
		class HandledEdge
		{
		public:
			size_t index;
			size_t associatedSCV[2];
			size_t scvfIndex;
			// indicates if the handled side is from or to
			bool from;
		};

	public:
	/// construct object and initialize local values and sizes
		HCRFVGeometry();

	///	update local data
		void update_local_data();

	/// update data for given element
		void update(GridObject* elem, const MathVector<worldDim>* vCornerCoords,
		            const ISubsetHandler* ish = NULL);

	///	debug output
		void print();
								   
		const MathVector<worldDim>* corners() const {return m_vCo;}							

	/// number of SubControlVolumeFaces
		inline size_t num_scvf() const {return numSCVF;};

	/// const access to SubControlVolumeFace number i
		inline const SCVF& scvf(size_t i) const
			{UG_ASSERT(i < maxNumSCVF, "Invalid Index."); return m_vSCVF[i];}

	/// number of SubControlVolumes
		inline size_t num_scv() const {return numSCV;}

	/// const access to SubControlVolume number i
		inline const SCV& scv(size_t i) const
			{UG_ASSERT(i < maxNumSCV, "Invalid Index."); return m_vSCV[i];}

	/// number of constrained dofs
		inline size_t num_constrained_dofs() const {return numConstrainedDofs;}

	/// const access to constrained dof i
		inline const CONSTRAINED_DOF& constrained_dof(size_t i) const
			{UG_ASSERT(i < numConstrainedDofs, "Invalid Index."); return m_vCD[i];}

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
		MathVector<worldDim> m_vGlobSCVF_IP[maxNumSCVF];
		MathVector<dim> m_vLocSCVF_IP[maxNumSCVF];
	 // coord of location for unknowns in faces (edge/face barycenter)
		MathVector<worldDim> m_vGlobUnkCoords[maxNumSCV];
		MathVector<dim> m_vLocUnkCoords[maxNumSCV];
		
		static const size_t numMaxCo = 8;
	 // corner coordinates
		MathVector<worldDim> m_vCo[numMaxCo];

	private:
		MathVector<dim> m_ipCoord[maxNumSCVF];

		///	SubControlVolumeFaces
		SCVF m_vSCVF[maxNumSCVF];

		///	SubControlVolumes
		SCV m_vSCV[maxNumSCV];

		/// constrained Dofs
		CONSTRAINED_DOF m_vCD[maxNumSCV];
	
		std::vector<HandledEdge> handledEdges;

		///	pointer to current element
		TElem* m_pElem;

		///	Reference Mapping
		ReferenceMapping<ref_elem_type, worldDim> m_mapping;

		///	Reference Element
		const ref_elem_type& m_rRefElem;

		///	Shape function set
		const local_shape_fct_set_type& m_rTrialSpace;
		
		size_t numSCV;
		size_t numSCVF;
		size_t numConstrainedDofs;
		// numDofs number of all dofs including constraining and constrained dofs
		size_t numDofs;

		bool localUpdateNecessary;
	
		static const size_t deleted = 117;

};


}

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__HANGING_CR_FINITE_VOLUME_GEOMETRY__ */
