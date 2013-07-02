/*
 * fe_geom.h
 *
 *  Created on: 08.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__

#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/reference_element/reference_mapping.h"
#include "common/util/provider.h"

#include <cmath>

namespace ug{

template <	typename TElem,	int TWorldDim,
			typename TTrialSpace, typename TQuadratureRule>
class FEGeometry
{
	public:
	///	type of reference element
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	/// reference element dimension
		static const int dim = ref_elem_type::dim;

	/// world dimension
		static const int worldDim = TWorldDim;

	///	type of trial space
		typedef TTrialSpace trial_space_type;

	///	type of quadrature rule
		typedef TQuadratureRule quad_rule_type;

	///	number of shape functions
		static const size_t nsh = trial_space_type::nsh;

	///	number of integration points
		static const size_t nip = quad_rule_type::nip;

	/// flag indicating if local data may change
		static const bool staticLocalData = true;

	public:
	///	Constructor
		FEGeometry();

	/// number of integration points
		size_t num_ip() const {return nip;}

	/// number of shape functions
		size_t num_sh() const {return nsh;}

	/// weight for integration point
		number weight(size_t ip) const {return m_vDetJ[ip] * m_rQuadRule.weight(ip);}

	/// local integration point
		const MathVector<dim>& local_ip(size_t ip) const {return m_rQuadRule.point(ip);}

	/// global integration point
		const MathVector<worldDim>& global_ip(size_t ip) const
		{
			UG_ASSERT(ip < nip, "Wrong ip.");
			return m_vIPGlobal[ip];
		}

	/// local integration point
		const MathVector<dim>* local_ips() const {return m_rQuadRule.points();}

	/// global integration point
		const MathVector<worldDim>* global_ips() const{return &m_vIPGlobal[0];}

		/// shape function at ip
		number shape(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < nip, "Wrong index"); UG_ASSERT(sh < nsh, "Wrong index");
			return m_vvShape[ip][sh];
		}

	/// local gradient at ip
		const MathVector<dim>& local_grad(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < nip, "Wrong index"); UG_ASSERT(sh < nsh, "Wrong index");
			return m_vvGradLocal[ip][sh];
		}

	/// global gradient at ip
		const MathVector<worldDim>& global_grad(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < nip, "Wrong index"); UG_ASSERT(sh < nsh, "Wrong index");
			return m_vvGradGlobal[ip][sh];
		}

	public:
	/// update Geometry for roid
		void update_local(ReferenceObjectID roid, const LFEID& lfeID, size_t orderQuad);
		void update_local(ReferenceObjectID roid, const LFEID& lfeID){
			update_local(roid, lfeID, 2*lfeID.order() + 1);
		}

	/// update Geometry for corners
		void update(GeometricObject* pElem, const MathVector<worldDim>* vCorner)
		{
			update(pElem, vCorner, LFEID(), m_rQuadRule.order());
		}

	/// update Geometry for corners
		void update(GeometricObject* pElem, const MathVector<worldDim>* vCorner,
		            const LFEID& lfeID, size_t orderQuad);
		void update(GeometricObject* pElem, const MathVector<worldDim>* vCorner,
		            const LFEID& lfeID){
			update(pElem, vCorner, lfeID, 2*lfeID.order() + 1);
		}

	protected:
	///	current element
		TElem* m_pElem;

	///	Quadrature rule
		const quad_rule_type& m_rQuadRule;

	///	Quadrature rule
		const trial_space_type& m_rTrialSpace;

	///	reference mapping
		ReferenceMapping<ref_elem_type, worldDim> m_mapping;

	///	global integration points
		MathVector<worldDim> m_vIPGlobal[nip];

	///	shape functions evaluated at ip
		number m_vvShape[nip][nsh];

	///	global gradient evaluated at ip
		MathVector<dim> m_vvGradLocal[nip][nsh];

	///	local gradient evaluated at ip
		MathVector<worldDim> m_vvGradGlobal[nip][nsh];

	///	jacobian of transformation at ip
		MathMatrix<worldDim,dim> m_vJTInv[nip];

	///	determinate of transformation at ip
		number m_vDetJ[nip];
};



template <int TWorldDim, int TRefDim = TWorldDim>
class DimFEGeometry
{
	public:
	/// reference element dimension
		static const int dim = TRefDim;

	/// world dimension
		static const int worldDim = TWorldDim;

	/// flag indicating if local data may change
		static const bool staticLocalData = false;

	public:
	///	default Constructor
		DimFEGeometry();

	///	Constructor
		DimFEGeometry(size_t order, LFEID lfeid);

	///	Constructor
		DimFEGeometry(ReferenceObjectID roid, size_t order, LFEID lfeid);

	/// number of integration points
		size_t num_ip() const {return m_nip;}

	/// number of shape functions
		size_t num_sh() const {return m_nsh;}

	/// weight for integration point
		number weight(size_t ip) const
		{
			UG_ASSERT(ip < m_nip, "Wrong index");
			return m_vDetJ[ip] * m_vQuadWeight[ip];
		}

	/// local integration point
		const MathVector<dim>& local_ip(size_t ip) const
		{
			UG_ASSERT(ip < m_nip, "Wrong index");
			return m_vIPLocal[ip];
		}

	/// global integration point
		const MathVector<worldDim>& global_ip(size_t ip) const
		{
			UG_ASSERT(ip < m_vIPGlobal.size(), "Wrong ip.");
			return m_vIPGlobal[ip];
		}

	/// local integration point
		const MathVector<dim>* local_ips() const {return m_vIPLocal;}

	/// global integration point
		const MathVector<worldDim>* global_ips() const{return &m_vIPGlobal[0];}

	/// shape function at ip
		number shape(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < m_vvShape.size(), "Wrong index");
			UG_ASSERT(sh < m_vvShape[ip].size(), "Wrong index");
			return m_vvShape[ip][sh];
		}

	/// local gradient at ip
		const MathVector<dim>& local_grad(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < m_vvGradLocal.size(), "Wrong index");
			UG_ASSERT(sh < m_vvGradLocal[ip].size(), "Wrong index");
			return m_vvGradLocal[ip][sh];
		}

	/// global gradient at ip
		const MathVector<worldDim>& global_grad(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < m_vvGradGlobal.size(), "Wrong index");
			UG_ASSERT(sh < m_vvGradGlobal[ip].size(), "Wrong index");
			return m_vvGradGlobal[ip][sh];
		}

	/// update Geometry for roid
		void update_local(ReferenceObjectID roid, const LFEID& lfeID, size_t orderQuad);
		void update_local(ReferenceObjectID roid, const LFEID& lfeID){
			update_local(roid, lfeID, 2*lfeID.order() + 1);
		}

	/// update Geometry for corners
		void update(GeometricObject* pElem, const MathVector<worldDim>* vCorner)
		{
			update(pElem, vCorner, m_lfeID, m_quadOrder);
		}

	/// update Geometry for corners
		void update(GeometricObject* pElem, const MathVector<worldDim>* vCorner,
					const LFEID& lfeID, size_t orderQuad);
		void update(GeometricObject* pElem, const MathVector<worldDim>* vCorner,
					const LFEID& lfeID){
			update(pElem, vCorner, lfeID, 2*lfeID.order() + 1);
		}

	public:
		///	boundary face
		class BF
		{
			public:
				BF() {}

				/// outer normal on bf. Norm is equal to area
				inline const MathVector<worldDim>& normal() const {return Normal;} // includes area

				/// volume of bf
				inline number volume() const {return Vol;}

				/// number of integration points on scvf
				inline size_t num_ip() const {return nip;}

				///	integration weight
				inline number weight(size_t ip) const
				{UG_ASSERT(ip<num_ip(), "Wrong index"); return vDetJ[ip] * vWeight[ip];}

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
				inline number detJ(size_t ip) const
				{UG_ASSERT(ip<num_ip(), "Wrong index"); return vDetJ[ip];}

				/// number of shape functions
				inline size_t num_sh() const {return nsh;}

				/// value of shape function i in integration point
				inline number shape(size_t ip, size_t sh) const
				{UG_ASSERT(ip<num_ip(), "Wrong index"); return vvShape[ip][sh];}

				/// vector of shape functions in ip point
				inline const number* shape_vector(size_t ip) const
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

			protected:
				/// let outer class access private members
				friend class DimFEGeometry<worldDim, dim>;

				size_t nip;
				std::vector<MathVector<dim> > vLocalIP; // local integration point (size: nip)
				std::vector<MathVector<worldDim> > vGlobalIP; // global integration point (size: nip)
				const number* vWeight; // weight at ip
				MathVector<worldDim> Normal; // normal (incl. area)
				number Vol; // volume of bf
				std::vector<MathMatrix<worldDim,dim> > vJtInv; // Jacobian transposed at ip (size: nip)
				std::vector<number> vDetJ; // Jacobian det at ip (size: nip)

				size_t nsh;
				std::vector<std::vector<number> > vvShape; // shapes at ip (size: nip x nsh)
				std::vector<std::vector<MathVector<dim> > > vvLocalGrad; // local grad at ip (size: nip x nsh)
				std::vector<std::vector<MathVector<worldDim> > > vvGlobalGrad; // global grad at ip (size: nip x nsh)
		};

		/// update boundary data for given element
		void update_boundary_faces(GeometricObject* pElem,
								   const MathVector<worldDim>* vCornerCoords,
								   size_t orderQuad,
								   const ISubsetHandler* ish);

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

	protected:
		std::map<int, std::vector<BF> > m_mapVectorBF;
		std::vector<BF> m_vEmptyVectorBF;

	protected:
	///	current reference object id the local values are prepared for
		ReferenceObjectID m_roid;

	///	current element
		GeometricObject* m_pElem;

	///	current integration order
		int m_quadOrder;

	///	current local finite element id
		LFEID m_lfeID;

	///	number of integration point
		size_t m_nip;

	///	local quadrature points
		const MathVector<dim>* m_vIPLocal;

	///	local quadrature weights
		const number* m_vQuadWeight;

	///	global integration points (size = nip)
		std::vector<MathVector<worldDim> > m_vIPGlobal;

	///	jacobian of transformation at ip (size = nip)
		std::vector<MathMatrix<worldDim,dim> > m_vJTInv;

	///	determinant of transformation at ip (size = nip)
		std::vector<number> m_vDetJ;

	///	number of shape functions
		size_t m_nsh;

	///	shape functions evaluated at ip (size = nip x nsh)
		std::vector<std::vector<number> > m_vvShape;

	///	global gradient evaluated at ip (size = nip x nsh)
		std::vector<std::vector<MathVector<dim> > > m_vvGradLocal;

	///	local gradient evaluated at ip (size = nip x nsh)
		std::vector<std::vector<MathVector<worldDim> > > m_vvGradGlobal;
};

} // end namespace ug

#include "fe_geom_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__ */
