/*
 * finite_element_geometry.h
 *
 *  Created on: 08.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__

#include "lib_discretization/quadrature/quadrature.h"
#include "lib_discretization/local_shape_function_set/local_shape_function_set_provider.h"
#include <cmath>

namespace ug{

template <	typename TElem,
			int TWorldDim,
			template <class TRefElem, int p> class TTrialSpace,
			int TOrderTrialSpace,
			template <class TRefElem, int p> class TQuadRule,
			int TOrderQuadRule>
class FEGeometry
{
	public:
	///	type of reference element
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	/// reference element dimension
		static const int dim = ref_elem_type::dim;

	/// world dimension
		static const int world_dim = TWorldDim;

	///	type of trial space
		typedef TTrialSpace<ref_elem_type, TOrderTrialSpace> trial_space_type;

	///	type of quadrature rule
		typedef TQuadRule<ref_elem_type, TOrderQuadRule> quad_rule_type;

	///	number of shape functions
		static const size_t nsh = trial_space_type::nsh;

	///	number of integration points
		static const size_t nip = quad_rule_type::nip;

	public:
	///	Constructor
		FEGeometry()
			: m_rQuadRule(QuadRuleProvider::get<quad_rule_type>()),
			  m_rTrialSpace(LSFSProvider::get<trial_space_type>())
		{
		//	evaluate local shapes and gradients
			for(size_t ip = 0; ip < nip; ++ip)
				for(size_t sh = 0; sh < nsh; ++sh)
				{
					m_vvShape[ip][sh] = m_rTrialSpace.shape(sh, m_rQuadRule.point(ip));
					m_vvGradLocal[ip][sh] = m_rTrialSpace.grad(sh, m_rQuadRule.point(ip));
				}
		}

	/// number of integration points
		size_t num_ip() const {return nip;}

	/// number of shape functions
		size_t num_sh() const {return nsh;}

	/// weight for integration point
		number weight(size_t ip) const {return fabs(m_detJ[ip]) * m_rQuadRule.weight(ip);}

	/// local integration point
		const MathVector<dim>& ip_local(size_t ip) const {return m_rQuadRule.point(ip);}

	/// global integration point
		const MathVector<world_dim>& ip_global(size_t ip) const
		{
			UG_ASSERT(ip < m_vIPGlobal.size(), "Wrong ip.");
			return m_vIPGlobal[ip];
		}

	/// local integration point
		const MathVector<dim>* local_ips() const {return m_rQuadRule.points();}

	/// global integration point
		const MathVector<world_dim>* global_ips() const{return &m_vIPGlobal[0];}

		/// shape function at ip
		number shape(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < nip, "Wrong index"); UG_ASSERT(sh < nsh, "Wrong index");
			return m_vvShape[ip][sh];
		}

	/// local gradient at ip
		const MathVector<dim>& grad_local(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < nip, "Wrong index"); UG_ASSERT(sh < nsh, "Wrong index");
			return m_vvGradLocal[ip][sh];
		}

	/// global gradient at ip
		const MathVector<world_dim>& grad_global(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < nip, "Wrong index"); UG_ASSERT(sh < nsh, "Wrong index");
			return m_vvGradGlobal[ip][sh];
		}

	/// update Geometry for corners
		bool update(const MathVector<world_dim>* vCorner)
		{
		//	update the mapping for the new corners
			m_mapping.update(vCorner);

		//	evaluate global data
			for(size_t ip = 0; ip < nip; ++ip)
			{
			// 	compute transformation inverse and determinate at ip
				if(!m_mapping.jacobian_transposed_inverse(m_vIPLocal[ip], m_JTInv[ip]))
				{
					UG_LOG("FE1Geometry:update(..): "
							"Cannot compute jacobian transposed inverse.\n");
					return false;
				}

			//	compute determinant
				if(!m_mapping.jacobian_det(m_vIPLocal[ip], m_detJ[ip]))
				{
					UG_LOG("FE1Geometry:update(..): "
							"Cannot compute jacobian determinante.\n");
					return false;
				}

			//	compute global integration points
				if(!m_mapping.local_to_global(m_vIPLocal[ip], m_vIPGlobal[ip]))
				{
					UG_LOG("FE1Geometry:update(..):"
							" Cannot compute global ip pos.\n");
					return false;
				}

			// 	compute global gradients
				for(size_t sh = 0; sh < nsh; ++sh)
				{
					MatVecMult(m_vvGradGlobal[ip][sh], m_JTInv[ip], m_vvGradLocal[ip][sh]);
				}
			}

		//	we're done
			return true;
		}

	protected:
	///	Quadrature rule
		const quad_rule_type& m_rQuadRule;

	///	Quadrature rule
		const trial_space_type& m_rTrialSpace;

	///	reference mapping
		ReferenceMapping<ref_elem_type, world_dim> m_mapping;

	///	local integration points
		MathVector<dim> m_vIPLocal[nip];

	///	global integration points
		MathVector<world_dim> m_vIPGlobal[nip];

	///	shape functions evaluated at ip
		number m_vvShape[nip][nsh];

	///	global gradient evaluated at ip
		MathVector<dim> m_vvGradLocal[nip][nsh];

	///	local gradient evaluated at ip
		MathVector<world_dim> m_vvGradGlobal[nip][nsh];

	///	jacobian of transformation at ip
		MathMatrix<world_dim,dim> m_JTInv[nip];

	///	determinate of transformation at ip
		number m_detJ[nip];
};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__ */
