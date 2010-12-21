/*
 * interpolation_operator.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__INTERPOLATION_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__INTERPOLATION_OPERATOR__

#include "common/common.h"

#include "lib_algebra/operator/operator_interface.h"
#include "lib_discretization/function_spaces/function_spaces.h"
#include "lib_discretization/operator/operator.h"
#include "lib_discretization/local_shape_function_set/local_shape_function_set_provider.h"
#include "lib_discretization/spatial_discretization/ip_data/user_data.h"
#include <boost/function.hpp>

namespace ug{

template <typename TElem, typename TGridFunction>
bool InterpolateFunctionOnElem( boost::function<void (	number& res,
														const MathVector<TGridFunction::domain_type::dim>& x,
														number& time)>
														InterpolFunction,
								TGridFunction& u, size_t fct, int si, number time)
{
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
	const int dim = ref_elem_type::dim;
	typedef typename TGridFunction::domain_type domain_type;

	LocalShapeFunctionSetID id = u.local_shape_function_set_id(fct);

	const LocalShapeFunctionSet<ref_elem_type>& trialSpace =
			LocalShapeFunctionSetProvider::get_local_shape_function_set<ref_elem_type>(id);

	// get local positions of interpolation points
	const size_t num_sh = trialSpace.num_sh();
	std::vector<MathVector<dim> > loc_pos(num_sh);
	for(size_t i = 0; i < num_sh; ++i)
	{
		if(!trialSpace.position(i, loc_pos[i]))
		{
			UG_LOG("Chosen Local Shape function Set does not provide "
					"senseful interpolation points. Can not use Lagrange interpolation.\n");
			return false;
		}
	}

	domain_type& domain = u.get_approximation_space().get_domain();
	typename domain_type::position_accessor_type aaPos = domain.get_position_accessor();
	typename domain_type::position_type corners[ref_elem_type::num_corners];

	ReferenceMapping<ref_elem_type, domain_type::dim> mapping;

	// iterate over all elements
	typename geometry_traits<TElem>::const_iterator iterBegin, iterEnd, iter;
	iterBegin = u.template begin<TElem>(si);
	iterEnd = u.template end<TElem>(si);

	number val;
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

		typename TGridFunction::vector_type& v_vec = u.get_vector();
		typename TGridFunction::multi_index_vector_type ind;
		u.get_multi_indices(elem, fct, ind);
		if(ind.size() != num_sh)
			{UG_LOG("Num dofs on element != num_sh.\n"); return false;}

		// get global positions
		for(size_t i = 0; i < num_sh; ++i)
		{
			typename domain_type::position_type glob_pos;

			if(!mapping.local_to_global(loc_pos[i], glob_pos))
				return false;

		//	evaluate functor
			InterpolFunction(val, glob_pos, time);

		//	set value
			BlockRef(v_vec[ind[i][0]], ind[i][1]) = val;
		}
	}
	return true;
}


template <typename TGridFunction>
bool InterpolateFunctionHelp(	boost::function<void (	number& res,
									const MathVector<TGridFunction::domain_type::dim>& x,
									number& time)>
									InterpolFunction,
							TGridFunction& u, const char* name, number time)
{
//	get Function Pattern
	const FunctionPattern& pattern = u.get_approximation_space().get_function_pattern();

//	get function id of name
	const size_t fct = pattern.fct_id_by_name(name);

//	check that function found
	if(fct == (size_t)-1)
	{
		UG_LOG("In InterpolateFunction: Name of function not found.\n");
		return false;
	}

//	check that function exists
	if(fct >= u.num_fct())
	{
		UG_LOG("Function space does not contain a function with index " << fct << ".\n");
		return false;
	}

//	loop for dimensions
	switch(u.dim(fct))
	{
	case 2:
		for(int si = 0; si < u.num_subsets(); ++si)
		{
			if(!u.is_def_in_subset(fct, si)) continue;

			if(!InterpolateFunctionOnElem<Triangle, TGridFunction>(InterpolFunction, u, fct, si, time)) return false;
			if(!InterpolateFunctionOnElem<Quadrilateral, TGridFunction>(InterpolFunction, u, fct, si, time)) return false;
		}
		break;
	case 3:
		for(int si = 0; si < u.num_subsets(); ++si)
		{
			if(u.is_def_in_subset(fct, si) != true) continue;

			// todo: more elements
			if(!InterpolateFunctionOnElem<Tetrahedron, TGridFunction>(InterpolFunction, u, fct, si, time)) return false;
			if(!InterpolateFunctionOnElem<Hexahedron, TGridFunction>(InterpolFunction, u, fct, si, time)) return false;
		}
		break;
	default: UG_LOG("InterpolateFunction: Dimension not implemented.\n"); return false;
	}

//	adjust parallel storage state
	#ifdef UG_PARALLEL
	u.set_storage_type(PST_CONSISTENT);
	#endif
	return true;

}

template <typename TGridFunction>
bool InterpolateFunction(	IUserData<number, TGridFunction::domain_type::dim>& InterpolFunctionProvider,
							TGridFunction& u, const char* name, number time)
{
//	extract functor
	static const int dim = TGridFunction::domain_type::dim;
	typedef typename IUserData<number, dim>::functor_type functor_type;
	functor_type InterpolFunction = InterpolFunctionProvider.get_functor();

	return InterpolateFunctionHelp(InterpolFunction, u, name, time);
}





// todo: we need a rework here !
template <typename TDiscreteFunction>
class LagrangeInterpolationOperator //: public IOperator<typename ContinuousFunctionSpace<typename TDiscreteFunction::domain_type>::function_type, TDiscreteFunction>
{
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
		bool prepare(codomain_function_type& v, domain_function_type& u)
		{
			return true;
		}

		// apply Operator, i.e. v = L(u);
		bool apply(codomain_function_type& v, const domain_function_type& u)
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

			#ifdef UG_PARALLEL
			v.set_storage_type(PST_CONSISTENT);
			#endif
			return true;
		}

		// apply Operator, i.e. v := v - L(u);
		bool apply_sub(codomain_function_type& v, const domain_function_type& u)
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
					LocalShapeFunctionSetProvider::get_local_shape_function_set<ref_elem_type>(id);

			// get local positions of interpolation points
			const size_t num_sh = trialSpace.num_sh();
			std::vector<MathVector<dim> > loc_pos(num_sh);
			for(size_t i = 0; i < num_sh; ++i)
			{
				if(trialSpace.position(i, loc_pos[i]) != true)
				{
					UG_LOG("Chosen Local Shape function Set does not provide senseful interpolation points. Can not use Lagrange interpolation.\n");
					return false;
				}
			}

			domain_type& domain = v.get_approximation_space().get_domain();
			typename domain_type::position_accessor_type aaPos = domain.get_position_accessor();
			typename domain_type::position_type corners[ref_elem_type::num_corners];

			ReferenceMapping<ref_elem_type, domain_type::dim> mapping;

			// iterate over all elements
			typename geometry_traits<TElem>::iterator iterBegin, iterEnd, iter;
			iterBegin = v.template begin<TElem>(si);
			iterEnd = v.template end<TElem>(si);

			number val;
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

				typename TDiscreteFunction::vector_type& v_vec = v.get_vector();
				typename TDiscreteFunction::multi_index_vector_type ind;
				v.get_multi_indices(elem, m_fct, ind);
				if(ind.size() != num_sh)
					{UG_LOG("Num dofs on element != num_sh.\n"); return false;}

				// get global positions
				for(size_t i = 0; i < num_sh; ++i)
				{
					typename domain_type::position_type glob_pos;
					//refElem.template mapLocalToGlobal<domain_type::dim>(corners, loc_pos[i], glob_pos);
					if(mapping.local_to_global(loc_pos[i], glob_pos) != true) return false;
					if(u(glob_pos, val) != true)
					{
						UG_LOG("Error while evaluating function u. Aborting interpolation.\n");
						return false;
					}
					BlockRef(v_vec[ind[i][0]], ind[i][1]) = val;
				}

			}
			return true;
		}

	protected:
		size_t m_fct;
};

} // namespace ug

#endif /*__H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__INTERPOLATION_OPERATOR__*/
