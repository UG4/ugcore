/*
 * interpolate.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACES__INTERPOLATE__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACES__INTERPOLATE__

#include "common/common.h"

#include "lib_discretization/common/subset_group.h"
#include "lib_discretization/domain_util.h"
#include "lib_discretization/local_shape_function_set/local_shape_function_set_provider.h"
#include "lib_discretization/spatial_discretization/ip_data/user_data.h"
#include <boost/function.hpp>

namespace ug{

/// interpolates a function on an element
template <typename TElem, typename TGridFunction>
bool InterpolateFunctionOnElem( boost::function<void (
									number& res,
									const MathVector<TGridFunction::domain_type::dim>& x,
									number& time)> InterpolFunction,
								TGridFunction& u, size_t fct, int si, number time)
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				ref_elem_type;

//	dimension of reference element
	const int dim = ref_elem_type::dim;

//	domain type and position_type
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_type position_type;

//	id of shape functions used
	LocalShapeFunctionSetID id = u.local_shape_function_set_id(fct);

//	get trial space
	const LocalShapeFunctionSet<ref_elem_type>& trialSpace =
			LocalShapeFunctionSetProvider::get_local_shape_function_set<ref_elem_type>(id);

//	number of dofs on element
	const size_t num_sh = trialSpace.num_sh();

// 	load local positions of dofs for the trial space on element
	std::vector<MathVector<dim> > loc_pos(num_sh);
	for(size_t i = 0; i < num_sh; ++i)
		if(!trialSpace.position(i, loc_pos[i]))
		{
			UG_LOG("ERROR in 'InterpolateFunctionOnElem': Cannot find meaningful"
					" local positions of dofs.\n");
			return false;
		}

//	create a reference mapping
	ReferenceMapping<ref_elem_type, domain_type::dim> mapping;

// 	iterate over all elements
	typename geometry_traits<TElem>::const_iterator iterEnd, iter;
	iterEnd = u.template end<TElem>(si);
	for(iter = u.template begin<TElem>(si); iter != iterEnd; ++iter)
	{
	//	get element
		TElem* elem = *iter;

	//	get all corner coordinates
		std::vector<position_type> vCorner;
		CollectCornerCoordinates(vCorner, *elem, u.get_domain());

	//	update the reference mapping for the corners
		mapping.update(&vCorner[0]);

//		typename TGridFunction::vector_type& v_vec = *dynamic_cast<typename TGridFunction::vector_type*>(&u);

	//	get multiindices of element
		typename TGridFunction::multi_index_vector_type ind;
		u.get_multi_indices(elem, fct, ind);

	//	check multi indices
		if(ind.size() != num_sh)
		{
			UG_LOG("ERROR in 'InterpolateFunctionOnElem': Wrong number of"
					" multi indices.\n");
			return false;
		}

	// 	loop all dofs
		for(size_t i = 0; i < num_sh; ++i)
		{
		//	global position
			position_type glob_pos;

		//  map local dof position to global position
			if(!mapping.local_to_global(loc_pos[i], glob_pos))
			{
				UG_LOG("ERROR in 'InterpolateFunctionOnElem': Cannot compute"
						" global dof position.\n");
				return false;
			}

		//	value at position
			number val;
			InterpolFunction(val, glob_pos, time);

		//	set value
			BlockRef(u[ind[i][0]], ind[i][1]) = val;
		}
	}

//	we're done
	return true;
}


template <typename TGridFunction>
bool InterpolateFunctionHelp(	boost::function<void (
									number& res,
									const MathVector<TGridFunction::domain_type::dim>& x,
									number& time)>
									InterpolFunction,
								TGridFunction& u,
								size_t fct,
								number time,
								const SubsetGroup& ssGrp)
{
//	loop subsets
	for(size_t i = 0; i < ssGrp.num_subsets(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	skip if function is not defined in subset
		if(!u.is_def_in_subset(fct, si)) continue;

	//	switch dimensions
		bool bRes = true;
		switch(u.dim(fct))
		{
		case 1:
			bRes &= InterpolateFunctionOnElem<Edge, TGridFunction>(InterpolFunction, u, fct, si, time);
			break;
		case 2:
			bRes &= InterpolateFunctionOnElem<Triangle, TGridFunction>(InterpolFunction, u, fct, si, time);
			bRes &= InterpolateFunctionOnElem<Quadrilateral, TGridFunction>(InterpolFunction, u, fct, si, time);
			break;
		case 3:
			bRes &= InterpolateFunctionOnElem<Tetrahedron, TGridFunction>(InterpolFunction, u, fct, si, time);
			bRes &= InterpolateFunctionOnElem<Hexahedron, TGridFunction>(InterpolFunction, u, fct, si, time);
			bRes &= InterpolateFunctionOnElem<Prism, TGridFunction>(InterpolFunction, u, fct, si, time);
			bRes &= InterpolateFunctionOnElem<Pyramid, TGridFunction>(InterpolFunction, u, fct, si, time);
			break;
		default: UG_LOG("InterpolateFunction: Dimension not implemented.\n"); return false;
		}

	//	check success
		if(!bRes)
		{
			UG_LOG("ERROR in InterpolateFunction: Cannot interpolate on elements.\n");
			return false;
		}
	}

//	adjust parallel storage state
#ifdef UG_PARALLEL
	u.set_storage_type(PST_CONSISTENT);
#endif

//	we're done
	return true;
}

/// interpolates a function on a subset
/**
 * This function interpolates a GridFunction. To evaluate the function on every
 * point a IUserData must be passed.
 *
 * \param[out]		u			interpolated grid function
 * \param[in]		data		data evaluator
 * \param[in]		name		symbolic name of function
 * \param[in]		time		time point
 * \param[in]		subsets		subsets, where to interpolate
 */
template <typename TGridFunction>
bool InterpolateFunction(	IUserData<number, TGridFunction::domain_type::dim>& data,
							TGridFunction& u, const char* name, number time,
							const char* subsets)
{
//	world dimension
	static const int dim = TGridFunction::domain_type::dim;

//	extract functor
	typedef typename IUserData<number, dim>::functor_type functor_type;
	functor_type InterpolFunction = data.get_functor();

//	get Function Pattern
	const typename TGridFunction::approximation_space_type& approxSpace
				= u.get_approximation_space();

//	get function id of name
	const size_t fct = approxSpace.fct_id_by_name(name);

//	check that function found
	if(fct == (size_t)-1)
	{
		UG_LOG("ERROR in InterpolateFunction: Name of function not found.\n");
		return false;
	}

//	check that function exists
	if(fct >= u.num_fct())
	{
		UG_LOG("ERROR in InterpolateFunction: Function space does not contain"
				" a function with index " << fct << ".\n");
		return false;
	}

//	create subset group
	SubsetGroup ssGrp; ssGrp.set_subset_handler(*approxSpace.get_subset_handler());

//	read subsets
	if(subsets != NULL)
		ConvertStringToSubsetGroup(ssGrp, *approxSpace.get_subset_handler(), subsets);
	else // add all if no subset specified
		ssGrp.add_all();

//	forward
	return InterpolateFunctionHelp(InterpolFunction, u, fct, time, ssGrp);
}

/// interpolates a function on the whole domain
template <typename TGridFunction>
bool InterpolateFunction(	IUserData<number, TGridFunction::domain_type::dim>& InterpolFunctionProvider,
							TGridFunction& u, const char* name, number time)
{
//	forward
	return InterpolateFunction(InterpolFunctionProvider, u, name, time, NULL);
}

} // namespace ug

#endif /*__H__LIBDISCRETIZATION__FUNCTION_SPACES__INTERPOLATE__*/
