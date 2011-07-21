/*
 * finite_element_geometry.h
 *
 *  Created on: 08.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__

#include "lib_discretization/quadrature/quadrature.h"
#include "lib_discretization/local_finite_element/local_shape_function_set.h"
#include "common/util/provider.h"

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
		static const int worldDim = TWorldDim;

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
			: m_rQuadRule(Provider::get<quad_rule_type>()),
			  m_rTrialSpace(Provider::get<trial_space_type>())
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
		number weight(size_t ip) const {return fabs(m_vDetJ[ip]) * m_rQuadRule.weight(ip);}

	/// local integration point
		const MathVector<dim>& ip_local(size_t ip) const {return m_rQuadRule.point(ip);}

	/// global integration point
		const MathVector<worldDim>& ip_global(size_t ip) const
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
		const MathVector<dim>& grad_local(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < nip, "Wrong index"); UG_ASSERT(sh < nsh, "Wrong index");
			return m_vvGradLocal[ip][sh];
		}

	/// global gradient at ip
		const MathVector<worldDim>& grad_global(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < nip, "Wrong index"); UG_ASSERT(sh < nsh, "Wrong index");
			return m_vvGradGlobal[ip][sh];
		}

	/// update Geometry for corners
		bool update(const MathVector<worldDim>* vCorner)
		{
		//	update the mapping for the new corners
			m_mapping.update(vCorner);

		//	compute global integration points
			for(size_t ip = 0; ip < nip; ++ip)
				m_mapping.local_to_global(m_vIPGlobal[ip], ip_local(ip));

		//	evaluate global data
		//	if reference mapping is linear,
			if(ReferenceMapping<ref_elem_type, worldDim>::isLinear)
			{
			// 	compute transformation inverse and determinate at first ip
				m_mapping.jacobian_transposed_inverse(m_vJTInv[0], ip_local(0));
				m_vDetJ[0] = m_mapping.jacobian_det(ip_local(0));

			//	copy values
				for(size_t ip = 1; ip < nip; ++ip)
				{
					m_vJTInv[ip] = m_vJTInv[0];
					m_vDetJ[ip] = m_vDetJ[0];
				}
			}
		//	else compute jacobian for each point
			{
				for(size_t ip = 0; ip < nip; ++ip)
				{
				// 	compute transformation inverse and determinate at ip
					m_mapping.jacobian_transposed_inverse(m_vJTInv[ip], ip_local(ip));

				//	compute determinant
					m_vDetJ[ip] = m_mapping.jacobian_det(ip_local(ip));
				}
			}

		// 	compute global gradients
			for(size_t ip = 0; ip < nip; ++ip)
				for(size_t sh = 0; sh < nsh; ++sh)
					MatVecMult(m_vvGradGlobal[ip][sh], m_vJTInv[ip], m_vvGradLocal[ip][sh]);

		//	we're done
			return true;
		}

	protected:
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

	public:
	///	Constructor
		DimFEGeometry() :
			m_roid(ROID_INVALID), m_order(0),
			m_vIPLocal(NULL), m_vQuadWeight(NULL)
		{}

	/// number of integration points
		size_t num_ip() const {return m_nip;}

	/// number of shape functions
		size_t num_sh() const {return m_nsh;}

	/// weight for integration point
		number weight(size_t ip) const
		{
			UG_ASSERT(ip < m_nip, "Wrong index");
			return fabs(m_vDetJ[ip]) * m_vQuadWeight[ip];
		}

	/// local integration point
		const MathVector<dim>& ip_local(size_t ip) const
		{
			UG_ASSERT(ip < m_nip, "Wrong index");
			return m_vIPLocal[ip];
		}

	/// global integration point
		const MathVector<worldDim>& ip_global(size_t ip) const
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
		const MathVector<dim>& grad_local(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < m_vvGradLocal.size(), "Wrong index");
			UG_ASSERT(sh < m_vvGradLocal[ip].size(), "Wrong index");
			return m_vvGradLocal[ip][sh];
		}

	/// global gradient at ip
		const MathVector<worldDim>& grad_global(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < m_vvGradGlobal.size(), "Wrong index");
			UG_ASSERT(sh < m_vvGradGlobal[ip].size(), "Wrong index");
			return m_vvGradGlobal[ip][sh];
		}

	/// update Geometry for corners
		bool update(GeometricObject* pElem, const MathVector<worldDim>* vCorner)
		{
		//	get reference element type
			ReferenceObjectID roid = pElem->reference_object_id();

			////////////////////////
			//	local values
			////////////////////////

		//	if already prepared for this roid, skip update of local values
			if(m_roid != roid)
			{
			//	remember current roid
				m_roid = roid;

			//	request for quadrature rule
				const QuadratureRule<dim>& quadRule
						= QuadratureRuleProvider<dim>::get_rule(roid, m_order);

			//	copy quad informations
				m_nip = quadRule.size();
				m_vIPLocal = quadRule.points();
				m_vQuadWeight = quadRule.weights();

			//	resize
				m_vIPGlobal.resize(m_nip);
				m_vJTInv.resize(m_nip);
				m_vDetJ.resize(m_nip);

				m_vvGradGlobal.resize(m_nip);
				m_vvGradLocal.resize(m_nip);
				m_vvShape.resize(m_nip);

			//	request for trial space
				const LocalShapeFunctionSet<dim>& lsfs
				 	 = LocalShapeFunctionSetProvider::get<dim>(roid, m_lfeID);

			//	get all shapes by on call
				for(size_t ip = 0; ip < m_nip; ++ip)
				{
					lsfs.shapes(m_vvShape[ip], m_vIPLocal[ip]);
					lsfs.grads(m_vvGradLocal[ip], m_vIPLocal[ip]);
				}
			}

			////////////////////////
			//	global values
			////////////////////////

/*		//	update the mapping for the new corners
			m_mapping.update(vCorner);

		//	compute global integration points
			for(size_t ip = 0; ip < nip; ++ip)
				m_mapping.local_to_global(m_vIPLocal[ip], m_vIPGlobal[ip]);

		//	evaluate global data
		//	if reference mapping is linear,
			if(ReferenceMapping<ref_elem_type, worldDim>::isLinear)
			{
			// 	compute transformation inverse and determinate at first ip
				m_mapping.jacobian_transposed_inverse(m_vIPLocal[0], m_JTInv[0]);
				m_detJ[0] = m_mapping.jacobian_det(m_vIPLocal[0]);

			//	copy values
				for(size_t ip = 1; ip < nip; ++ip)
				{
					m_JTInv[ip] = m_JTInv[0];
					m_detJ[ip] = m_detJ[0];
				}
			}
		//	else compute jacobian for each point
			{
				for(size_t ip = 0; ip < nip; ++ip)
				{
				// 	compute transformation inverse and determinate at ip
					m_mapping.jacobian_transposed_inverse(m_vIPLocal[ip], m_JTInv[ip]);

				//	compute determinant
					m_detJ[ip] = m_mapping.jacobian_det(m_vIPLocal[ip]);
				}
			}

		// 	compute global gradients
			for(size_t ip = 0; ip < nip; ++ip)
				for(size_t sh = 0; sh < nsh; ++sh)
					MatVecMult(m_vvGradGlobal[ip][sh], m_JTInv[ip], m_vvGradLocal[ip][sh]);
*/
		//	we're done
			return true;
		}

	protected:
	///	current reference object id the local values are prepared for
		ReferenceObjectID m_roid;

	///	current integration order
		int m_order;

	///	current local finite element id
		LFEID m_lfeID;

	///	number of intergration point
		size_t m_nip;

	///	local quadrature weights
		const number* m_vQuadWeight;

	///	local quadrature points
		const MathVector<dim>* m_vIPLocal;

	///	global integration points (size = nip)
		std::vector<MathVector<worldDim> > m_vIPGlobal;

	///	jacobian of transformation at ip (size = nip)
		std::vector<MathMatrix<worldDim,dim> > m_vJTInv;

	///	determinate of transformation at ip (size = nip)
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

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__ */
