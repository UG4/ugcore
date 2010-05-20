/*
 * interpolation_operator.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__INTERPOLATION_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__INTERPOLATION_OPERATOR__

#include "common/common.h"

#include "lib_discretization/function_spaces/function_spaces.h"
#include "lib_discretization/operator/operator.h"
#include "lib_discretization/local_shape_function_set/local_shape_function_set_factory.h"

namespace ug{

template <typename TDiscreteFunction>
class LagrangeInterpolationOperator : public IDiscreteLinearOperator<typename ContinuousFunctionSpace<typename TDiscreteFunction::domain_type>::function_type, TDiscreteFunction> {
	protected:
		// type of physical domain
		typedef typename TDiscreteFunction::domain_type domain_type;

	public:
		// domain space
		typedef typename ContinuousFunctionSpace<domain_type>::function_type domain_function_type;

		// range space
		typedef TDiscreteFunction  codomain_function_type;

	public:
		LagrangeInterpolationOperator() : m_fct(-1) {};

		// Init Operator
		bool init(uint fct)
		{
			m_fct = fct;

			return true;
		}

		// implement interface
		bool init()
		{
			// should not be called
			return false;
		}

		// prepare Operator
		bool prepare(domain_function_type& u, codomain_function_type& v)
		{
			return true;
		}

		// apply Operator, i.e. v = L(u);
		bool apply(domain_function_type& u, codomain_function_type& v)
		{
			uint level = v.get_level();
			LocalShapeFunctionSetID id = v.get_local_shape_function_set_id(m_fct);
			if(id == LSFS_INVALID)
			{
				UG_LOG("Fundamental discrete function nr " << m_fct << " is not contained in Approximation Space. Can not determine Local Shape Function Set.\n");
				return false;
			}
			if(v.dim_fct(m_fct) == 2)
			{
				for(int subsetIndex = 0; subsetIndex < v.num_subsets(); ++subsetIndex)
				{
					if(v.fct_is_def_in_subset(m_fct, subsetIndex) != true) continue;
					if(interpolate_values<Triangle>(u, v, id, level, subsetIndex) != true) return false;
					if(interpolate_values<Quadrilateral>(u, v, id, level, subsetIndex) != true) return false;
				}
			}
			else
			{
				UG_LOG("Currently only 2D interpolation implemented.\n");
				return false;
			}

			return true;
		}

		// apply Operator, i.e. v := v - L(u);
		bool apply_sub(domain_function_type& u, codomain_function_type& v)
		{
			UG_ASSERT(0, "Not Implemented.");
			return true;
		}

	protected:
		template <typename TElem>
		bool interpolate_values(const domain_function_type& u, codomain_function_type& v, LocalShapeFunctionSetID id, uint level, int subsetIndex)
		{
			const int refDim = reference_element_traits<TElem>::dim;
			typedef typename reference_element_traits<TElem>::reference_element_type reference_element_type;



			const LocalShapeFunctionSet<reference_element_type>& trialSpace = LocalShapeFunctionSetFactory::inst().get_local_shape_function_set<reference_element_type>(id);
			typename geometry_traits<TElem>::iterator iterBegin, iterEnd, iter;

			iterBegin = v.template begin<TElem>(subsetIndex);
			iterEnd = v.template end<TElem>(subsetIndex);

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

			domain_type& domain = v.get_domain();
			typename domain_type::position_accessor_type aaPos = domain.get_position_accessor();
			typename domain_type::position_type corners[reference_element_traits<TElem>::num_corners];

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
					typename domain_type::position_type glob_pos;
					refElem.template mapLocalToGlobal<domain_type::dim>(corners, loc_pos[i], glob_pos);
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
		uint m_fct;
};

} // namespace ug

#endif /*__H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__INTERPOLATION_OPERATOR__*/
