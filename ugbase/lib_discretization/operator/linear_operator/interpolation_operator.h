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
class LagrangeInterpolationOperator : public ILinearOperator<typename ContinuousFunctionSpace<typename TDiscreteFunction::domain_type>::function_type, TDiscreteFunction> {
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
		bool init(size_t fct)
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
			if(m_fct >= v.num_fct())
			{
				UG_LOG("Function space does not contain a function with index " << m_fct << ".\n");
				return false;
			}
			if(v.dim(m_fct) == 2)
			{
				for(int si = 0; si < v.num_subsets(); ++si)
				{
					if(v.is_def_in_subset(m_fct, si) != true) continue;

					if(interpolate_values<Triangle>(u, v, si) != true) return false;
					if(interpolate_values<Quadrilateral>(u, v, si) != true) return false;
				}
			}
			else if(v.dim(m_fct) == 3)
			{
				for(int si = 0; si < v.num_subsets(); ++si)
				{
					if(v.is_def_in_subset(m_fct, si) != true) continue;

					if(interpolate_values<Triangle>(u, v, si) != true) return false;
					if(interpolate_values<Quadrilateral>(u, v, si) != true) return false;
				}
			}
			else
			{
				UG_LOG("Not implemented.\n");
				return false;
			}

			v.set_storage_type(PST_CONSISTENT);
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
		bool interpolate_values(const domain_function_type& u, codomain_function_type& v, int si)
		{
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
			const int dim = ref_elem_type::dim;

			LocalShapeFunctionSetID id = v.local_shape_function_set_id(m_fct);

			const LocalShapeFunctionSet<ref_elem_type>& trialSpace =
					LocalShapeFunctionSetFactory::inst().get_local_shape_function_set<ref_elem_type>(id);

			// get local positions of interpolation points
			const size_t num_sh = trialSpace.num_shape_functions();
			std::vector<MathVector<dim> > loc_pos(num_sh);
			for(size_t i = 0; i < num_sh; ++i)
			{
				if(trialSpace.position_of_dof(i, loc_pos[i]) != true)
				{
					UG_LOG("Chosen Local Shape function Set does not provide senseful interpolation points. Can not use Lagrange interpolation.\n");
					return false;
				}
			}

			domain_type& domain = v.get_approximation_space().get_domain();
			typename domain_type::position_accessor_type aaPos = domain.get_position_accessor();
			typename domain_type::position_type corners[ref_elem_type::num_corners];

			typename codomain_function_type::local_vector_type val(num_sh);
			//reference_element_type refElem;
			ReferenceMapping<ref_elem_type, domain_type::dim> mapping;

			// iterate over all elements
			typename geometry_traits<TElem>::iterator iterBegin, iterEnd, iter;
			iterBegin = v.template begin<TElem>(si);
			iterEnd = v.template end<TElem>(si);
			for(iter = iterBegin; iter != iterEnd; ++iter)
			{
				TElem* elem = *iter;

				// load corners of this element
				for(int i = 0; i < ref_elem_type::num_corners; ++i)
				{
					VertexBase* vert = elem->vertex(i);
					corners[i] = aaPos[vert];
				}

				mapping.update(corners);

				// get global positions
				for(size_t i = 0; i < num_sh; ++i)
				{
					typename domain_type::position_type glob_pos;
					//refElem.template mapLocalToGlobal<domain_type::dim>(corners, loc_pos[i], glob_pos);
					if(mapping.local_to_global(loc_pos[i], glob_pos) != true) return false;
					if(u(glob_pos, val[i]) != true)
					{
						UG_LOG("Error while evaluating function u. Aborting interpolation.\n");
						return false;
					}
				}

				// set values
				if(v.set_dof_values(elem, m_fct, val) != true)
				{
					UG_LOG("Error while writing values of degrees of freedom. Aborting interpolation.\n");
					return false;
				}
			}
			return true;
		}

	protected:
		size_t m_fct;
};

} // namespace ug

#endif /*__H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__INTERPOLATION_OPERATOR__*/
