/*
 * linear_operator.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__LINEAR_OPERATOR__INTERPOLATION_OPERATOR__
#define __H__LIBDISCRETIZATION__LINEAR_OPERATOR__INTERPOLATION_OPERATOR__

#include "common/common.h"

#include "lib_discretization/function_spaces/function_spaces.h"
#include "lib_discretization/linear_operator/linear_operator.h"
#include "lib_discretization/local_shape_function_set/local_shape_function_set_factory.h"

namespace ug{

template <typename TDomain, typename TAlgebra, typename TDoFManager>
class LagrangeInterpolationOperator : public LinearOperator<ContinuousFunctionSpace<TDomain>, ApproximationSpace<TDomain, TAlgebra, TDoFManager> > {
	public:
		// domain space
		typedef ContinuousFunctionSpace<TDomain> domain_type;

		// range space
		typedef ApproximationSpace<TDomain,TAlgebra,TDoFManager>  codomain_type;

		// domain function type
		typedef typename domain_type::function_type domain_function_type;

		// codomain function type
		typedef typename codomain_type::function_type codomain_function_type;

	public:
		LagrangeInterpolationOperator() : m_id(LSFS_INVALID) {};

		// Init Operator
		bool init(uint fct, uint level, ApproximationSpace<TDomain, TAlgebra, TDoFManager>& approxSpace)
		{
			m_approxSpace = &approxSpace;

			m_fct = fct;
			m_level = level;
			m_id = approxSpace.get_local_shape_function_set_id(fct);

			if(m_id == LSFS_INVALID)
			{
				UG_LOG("Fundamental discrete function nr " << fct << " is not contained in Approximation Space. Can not determine Local Shape Function Set.\n");
				return false;
			}

			return true;
		}

		// prepare Operator
		bool prepare()
		{
			return true;
		}

		// apply Operator, i.e. v = L(u);
		bool apply(const domain_function_type& u, codomain_function_type& v)
		{
			if(m_id == LSFS_INVALID)
			{
				UG_LOG("Operator not correctly prepared.\n");
				return false;
			}

			if(v.dim_fct(m_fct) == 2)
			{
				for(int subsetIndex = 0; subsetIndex < v.num_subsets(); ++subsetIndex)
				{
					if(v.fct_is_def_in_subset(m_fct, subsetIndex) != true) continue;

					if(interpolate_values<Triangle>(u, v, m_level, subsetIndex) != true) return false;
					if(interpolate_values<Quadrilateral>(u, v, m_level, subsetIndex) != true) return false;
				}
			}
			else
			{
				UG_LOG("Currently only 2D interpolation implemented.\n");
				return false;
			}

			return true;
		}

	protected:
		template <typename TElem>
		bool interpolate_values(const domain_function_type& u, codomain_function_type& v, uint level, int subsetIndex)
		{
			const int refDim = reference_element_traits<TElem>::dim;
			typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;
			const LocalShapeFunctionSet<reference_element_type>& trialSpace = LocalShapeFunctionSetFactory::inst().get_local_shape_function_set<reference_element_type>(m_id);
			typename geometry_traits<TElem>::iterator iterBegin, iterEnd, iter;

			iterBegin = v.template begin<TElem>(level, subsetIndex);
			iterEnd = v.template end<TElem>(level, subsetIndex);

			// get local positions of interpolation points
			const int num_sh = trialSpace.num_shape_functions();
			std::vector<MathVector<refDim> > loc_pos(num_sh);
			for(int i = 0; i < num_sh; ++i)
			{
				if(trialSpace.position_of_dof(i, loc_pos[i]) != true)
				{
					UG_LOG("Chosen Local Shape function Set does not provide senseful interpolation points. Can not use Lagrange interpolation.\n");
					return false;
				}
			}

			TDomain& domain = m_approxSpace->get_domain();
			typename TDomain::position_accessor_type aaPos = domain.get_position_accessor();
			typename TDomain::position_type corners[reference_element_traits<TElem>::num_corners];

			reference_element_type refElem;
			typename codomain_function_type::local_vector_type val(num_sh);

			// iterate over all elements
			for(iter = iterBegin; iter != iterEnd; ++iter)
			{
				TElem* elem = *iter;

				// load corners of this element
				for(int i = 0; i < reference_element_traits<TElem>::num_corners; ++i)
				{
					VertexBase* vert = elem->vertex(i);
					corners[i] = aaPos[vert];
				}

				for(int i = 0; i < num_sh; ++i)
				{
					typename TDomain::position_type glob_pos;
					refElem.template mapLocalToGlobal<TDomain::dim>(corners, loc_pos[i], glob_pos);
					if(u(glob_pos, val[i]) != true)
					{
						UG_LOG("Error while evaluating function u. Aborting interpolation.\n");
						return false;
					}
				}
				if(v.set_dof_values(elem, m_fct, val) != true)
				{
					UG_LOG("Error while writing values of degrees of freedom. Aborting interpolation.\n");
					return false;
				}
			}
			return true;
		}

	protected:
		ApproximationSpace<TDomain, TAlgebra, TDoFManager>* m_approxSpace;
		LocalShapeFunctionSetID m_id;
		uint m_fct;
		uint m_level;
};

} // namespace ug

#endif
