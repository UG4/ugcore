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
class FEGeometry
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
		FEGeometry(int order)
			: m_rQuadRule(QuadratureRuleFactory<ref_elem_type>::get_quadrature_rule(order))
			{
				UG_ASSERT(order == 1, " Currently only first order implemented.");
				const LocalShapeFunctionSet<ref_elem_type>& TrialSpace =
						LocalShapeFunctionSetFactory::inst().get_local_shape_function_set<ref_elem_type>(LSFS_LAGRANGEP1);

				// number of ips
				m_numIP = m_rQuadRule.num_points();
				m_numSH = TrialSpace.num_sh();

				// resize
				m_vIPLocal.resize(m_numIP);
				m_vIPGlobal.resize(m_numIP);
				m_JTInv.resize(m_numIP);
				m_detJ.resize(m_numIP);
				m_vvShape.resize(m_numIP);
				m_vvGradLocal.resize(m_numIP);
				m_vvGradGlobal.resize(m_numIP);
				for(size_t ip = 0; ip < m_numIP; ++ip)
				{
					m_vvShape[ip].resize(m_numSH);
					m_vvGradLocal[ip].resize(m_numSH);
					m_vvGradGlobal[ip].resize(m_numSH);
				}

				// write local Data
				for(size_t ip = 0; ip < m_numIP; ++ip)
				{
					m_vIPLocal[ip] = m_rQuadRule.point(ip);
					for(size_t sh = 0; sh < m_numSH; ++sh)
					{
						if(!TrialSpace.evaluate(sh, m_vIPLocal[ip], m_vvShape[ip][sh]))
							{UG_ASSERT(0, "Cannot evaluate shape.");}
						if(!TrialSpace.evaluate_grad(sh, m_vIPLocal[ip], m_vvGradLocal[ip][sh]))
							{UG_ASSERT(0, "Cannot evaluate gradient.");}
					}
				}
			}

		// number of integration points
		size_t num_ip() const {return m_numIP;}

		// number of shape functions
		size_t num_sh() const {return m_numSH;}

		// weight for integration point
		number weight(size_t ip) const
		{
			UG_ASSERT(ip < m_detJ.size(), "Wrong ip.");
			return m_detJ[ip];
		}

		/// local integration point
		const MathVector<dim>& ip_local(size_t ip) const
		{
			UG_ASSERT(ip < m_vIPLocal.size(), "Wrong ip.");
			return m_vIPLocal[ip];
		}

		/// global integration point
		const MathVector<world_dim>& ip_global(size_t ip) const
		{
			UG_ASSERT(ip < m_vIPGlobal.size(), "Wrong ip.");
			return m_vIPGlobal[ip];
		}

		/// global gradient at ip
		const MathVector<world_dim>& grad_global(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < m_vvGradGlobal.size(), "Wrong ip.");
			UG_ASSERT(sh < m_vvGradGlobal[ip].size(), "Wrong sh.");
			return m_vvGradGlobal[ip][sh];
		}

		/// shape function at ip
		number shape(size_t ip, size_t sh) const
		{
			UG_ASSERT(ip < m_vvShape.size(), "Wrong ip.");
			UG_ASSERT(sh < m_vvShape[ip].size(), "Wrong sh.");
			return m_vvShape[ip][sh];
		}

		/// update Geometry for corners
		bool update(MathVector<world_dim> corners[])
		{
			m_mapping.update(corners);

			// write local Data
			for(size_t ip = 0; ip < m_numIP; ++ip)
			{
				// compute transformation inverse and determinate at ip
				if(!m_mapping.jacobian_transposed_inverse(m_vIPLocal[ip], m_JTInv[ip]))
					{UG_LOG("FE1Geometry:update(..): Cannot compute jacobian transposed inverse.\n"); return false;}
				if(!m_mapping.jacobian_det(m_vIPLocal[ip], m_detJ[ip]))
					{UG_LOG("FE1Geometry:update(..): Cannot compute jacobian determinante.\n"); return false;}
				if(!m_mapping.local_to_global(m_vIPLocal[ip], m_vIPGlobal[ip]))
					{UG_LOG("FE1Geometry:update(..): Cannot compute global ip pos.\n"); return false;}

				// compute global gradients
				for(size_t sh = 0; sh < m_numSH; ++sh)
				{
					MatVecMult(m_vvGradGlobal[ip][sh], m_JTInv[ip], m_vvGradLocal[ip][sh]);
				}
			}
			return true;
		}

	protected:
		const QuadratureRule<ref_elem_type>& m_rQuadRule;

		ReferenceMapping<ref_elem_type, world_dim> m_mapping;

		size_t m_numIP;
		size_t m_numSH;
		std::vector<MathVector<dim> > m_vIPLocal;
		std::vector<MathVector<world_dim> > m_vIPGlobal;

		std::vector<std::vector<number> > m_vvShape;
		std::vector<std::vector<MathVector<dim> > > m_vvGradLocal;
		std::vector<std::vector<MathVector<world_dim> > > m_vvGradGlobal;

		std::vector<MathMatrix<world_dim,dim> > m_JTInv;
		std::vector<number> m_detJ;
};

// Singeton
template <	typename TElem,
			int TWorldDim = reference_element_traits<TElem>::reference_element_type::dim>
class FEGeometryProvider
{
	private:
		FEGeometryProvider(){};
		FEGeometryProvider(const FEGeometryProvider&){};
		FEGeometryProvider& operator=(const FEGeometryProvider&);

		template <int TOrder>
		static FEGeometry<TElem, TWorldDim>& geom()
		{
			static FEGeometry<TElem, TWorldDim> inst(TOrder);
			return inst;
		}

	public:
		static FEGeometry<TElem, TWorldDim>& get_geom(int order)
		{
			UG_ASSERT(order == 1, "Currently only order 1 implemented.");
			return geom<1>();
		}
};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_ELEMENT_GEOMETRY__ */
