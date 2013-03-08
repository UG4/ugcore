/*
 * fe_geom.h
 *
 *  Created on: 08.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__

#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
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
		void update_local(ReferenceObjectID roid, LFEID lfeID, size_t orderQuad);

	/// update Geometry for corners
		void update(TElem* pElem, const MathVector<worldDim>* vCorner)
		{
			update(pElem, vCorner, m_rTrialSpace.type(), m_rQuadRule.order());
		}

	/// update Geometry for corners
		void update(TElem* pElem, const MathVector<worldDim>* vCorner,
		            LFEID lfeID, size_t orderQuad);

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
		void update_local(ReferenceObjectID roid, LFEID lfeID, size_t orderQuad);

	/// update Geometry for corners
		void update(GeometricObject* pElem, const MathVector<worldDim>* vCorner)
		{
			update(pElem, vCorner, m_lfeID, m_quadOrder);
		}

	/// update Geometry for corners
		void update(GeometricObject* pElem, const MathVector<worldDim>* vCorner,
					LFEID lfeID, size_t orderQuad);

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
