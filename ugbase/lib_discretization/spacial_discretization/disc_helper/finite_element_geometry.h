/*
 * finite_element_geometry.h
 *
 *  Created on: 08.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__

namespace ug{

template <	typename TElem,
			int TWorldDim = reference_element_traits<TElem>::reference_element_type::dim>
class FE1Geometry
{
	public:
		/// reference dim
		static const int dim = reference_element_traits<TElem>::reference_element_type::dim;

		/// world dim
		static const int world_dim = TWorldDim;

	protected:
		// reference element
		typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

	public:
		FE1Geometry()
			: m_rQuadRule(QuadratureRuleFactory<ref_elem_type>::get_quadrature_rule(1))
			{
				const LocalShapeFunctionSet<ref_elem_type>& TrialSpace =
						LocalShapeFunctionSetFactory::inst().get_local_shape_function_set<ref_elem_type>(LSFS_LAGRANGEP1);

				m_ip_local = m_rQuadRule.point(0);

				for(size_t sh = 0; sh < (size_t) ref_elem_type::num_corners; ++sh)
				{
					if(!TrialSpace.evaluate(sh, m_ip_local, m_shape[sh]))
						{UG_ASSERT(0, "Cannot evaluate shape.");}
					if(!TrialSpace.evaluate_grad(sh, m_ip_local, m_grad_local[sh]))
						{UG_ASSERT(0, "Cannot evaluate gradient.");}
				}
			}

		// number of integration points
		size_t num_ip() const {return 1;}

		// number of shape functions
		size_t num_sh() const {return ref_elem_type::num_corners;}

		// weight for integration point
		number weight(size_t ip) const
		{
			UG_ASSERT(ip == 0, "Wrong ip.");
			return m_detJ;
		}

		/// local integration point
		const MathVector<dim>& ip_local(size_t ip) const
		{
			UG_ASSERT(ip == 0, "Wrong ip.");
			return m_ip_local;
		}

		/// global integration point
		const MathVector<world_dim>& ip_global(size_t ip) const
		{
			UG_ASSERT(ip == 0, "Wrong ip.");
			return m_ip_global;
		}

		/// global gradient at ip
		const MathVector<world_dim>& grad_global(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip == 0, "Wrong ip.");
			UG_ASSERT(sh < num_sh(), "Wrong sh.");
			return m_grad_global[sh];
		}

		/// shape function at ip
		number shape(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip == 0, "Wrong ip.");
			UG_ASSERT(sh < num_sh(), "Wrong sh.");
			return m_shape[sh];
		}

		/// update Geometry for corners
		bool update(MathVector<world_dim> corners[])
		{
			m_mapping.update(corners);

			// compute transformation inverse and determinate at ip
			if(!m_mapping.jacobian_transposed_inverse(m_ip_local, m_JTInv))
				{UG_LOG("FE1Geometry:update(..): Cannot compute jacobian transposed inverse.\n"); return false;}
			if(!m_mapping.jacobian_det(m_ip_local, m_detJ))
				{UG_LOG("FE1Geometry:update(..): Cannot compute jacobian determinante.\n"); return false;}
			if(!m_mapping.local_to_global(m_ip_local, m_ip_global))
			{UG_LOG("FE1Geometry:update(..): Cannot compute global ip pos.\n"); return false;}

			// compute Shape Gradient
			for(size_t sh = 0; sh < (size_t) ref_elem_type::num_corners; ++sh)
			{
				MatVecMult(m_grad_global[sh], m_JTInv, m_grad_local[sh]);
			}

			return true;
		}

	protected:
		const QuadratureRule<ref_elem_type>& m_rQuadRule;

		ReferenceMapping<ref_elem_type, world_dim> m_mapping;

		MathVector<dim> m_ip_local;
		MathVector<world_dim> m_ip_global;

		number m_shape[ref_elem_type::num_corners];
		MathVector<dim> m_grad_local[ref_elem_type::num_corners];
		MathVector<world_dim> m_grad_global[ref_elem_type::num_corners];

		MathMatrix<world_dim,dim> m_JTInv;
		number m_detJ;
};

// Singeton
template <	typename TElem,
			int TWorldDim = reference_element_traits<TElem>::reference_element_type::dim>
class FE1Geom
{
	private:
		FE1Geom(){};
		FE1Geom(const FE1Geom&){};
		FE1Geom& operator=(const FE1Geom&);

	public:
		static FE1Geometry<TElem, TWorldDim>& instance()
		{
			static FE1Geometry<TElem, TWorldDim> inst;
			return inst;
		}
};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__ */
