/*
 * interpolate.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACES__INTERPOLATE__
#define __H__UG__LIB_DISC__FUNCTION_SPACES__INTERPOLATE__

#include "common/common.h"

#include "lib_disc/common/subset_group.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include <boost/function.hpp>
#include "lib_disc/spatial_disc/ip_data/const_user_data.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{

/// interpolates a function on an element
template <typename TElem, typename TGridFunction>
void InterpolateFunctionOnElem(
		IPData<number, TGridFunction::domain_type::dim>& InterpolFunction,
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

//	get iterators
	typename TGridFunction::template traits<TElem>::const_iterator iterEnd, iter;
	iterEnd = u.template end<TElem>(si);
	iter = u.template begin<TElem>(si);

//	check if something to do:
	if(iter == iterEnd) return;

//	id of shape functions used
	LFEID id = u.local_finite_element_id(fct);

//	get trial space
	const LocalShapeFunctionSet<ref_elem_type>& trialSpace =
			LocalShapeFunctionSetProvider::get<ref_elem_type>(id);

//	number of dofs on element
	const size_t nsh = trialSpace.num_sh();

// 	load local positions of dofs for the trial space on element
	std::vector<MathVector<dim> > loc_pos(nsh);
	for(size_t i = 0; i < nsh; ++i)
		if(!trialSpace.position(i, loc_pos[i]))
			UG_THROW_FATAL("InterpolateFunctionOnElem: Cannot find meaningful"
					" local positions of dofs.");

//	create a reference mapping
	ReferenceMapping<ref_elem_type, domain_type::dim> mapping;

// 	iterate over all elements
	for( ; iter != iterEnd; ++iter)
	{
	//	get element
		TElem* elem = *iter;

	//	get all corner coordinates
		std::vector<position_type> vCorner;
		CollectCornerCoordinates(vCorner, *elem, *u.domain());

	//	update the reference mapping for the corners
		mapping.update(&vCorner[0]);

	//	get multiindices of element
		std::vector<MultiIndex<2> > ind;
		u.multi_indices(elem, fct, ind);

	//	check multi indices
		if(ind.size() != nsh)
			UG_THROW_FATAL("InterpolateFunctionOnElem: Number of shapes is "
					<<nsh<<", but got "<<ind.size()<<" multi indices.");

	// 	loop all dofs
		for(size_t i = 0; i < nsh; ++i)
		{
		//	global position
			position_type glob_pos;

		//  map local dof position to global position
			mapping.local_to_global(glob_pos, loc_pos[i]);

		//	value at position
			number val;
			InterpolFunction(val, glob_pos, time, si);

		//	set value
			BlockRef(u[ind[i][0]], ind[i][1]) = val;
		}
	}
}

template <int dim, typename TGridFunction>
struct InterpolationLoopHelp
{};

template <typename TGridFunction>
struct InterpolationLoopHelp<3, TGridFunction>{
static void invoke(IPData<number, TGridFunction::domain_type::dim>& InterpolFunction,
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
		try
		{
		const int dim = ssGrp.dim(i);
		switch(dim)
		{
		case DIM_SUBSET_EMPTY_GRID: break;
		case 0: /* \TODO: do nothing may be wrong */	break;
		case 1:
			InterpolateFunctionOnElem<Edge, TGridFunction>(InterpolFunction, u, fct, si, time);
			break;
		case 2:
			InterpolateFunctionOnElem<Triangle, TGridFunction>(InterpolFunction, u, fct, si, time);
			InterpolateFunctionOnElem<Quadrilateral, TGridFunction>(InterpolFunction, u, fct, si, time);
			break;
		case 3:
			InterpolateFunctionOnElem<Tetrahedron, TGridFunction>(InterpolFunction, u, fct, si, time);
			InterpolateFunctionOnElem<Hexahedron, TGridFunction>(InterpolFunction, u, fct, si, time);
			InterpolateFunctionOnElem<Prism, TGridFunction>(InterpolFunction, u, fct, si, time);
			InterpolateFunctionOnElem<Pyramid, TGridFunction>(InterpolFunction, u, fct, si, time);
			break;
		default: UG_THROW_FATAL("InterpolateFunction: Dimension " <<dim<<
		                " not possible for world dim "<<3<<".");
		}
		}
		UG_CATCH_THROW("InterpolateFunction: Cannot interpolate on elements.");
	}

//	adjust parallel storage state
#ifdef UG_PARALLEL
	u.set_storage_type(PST_CONSISTENT);
#endif
}
};

template <typename TGridFunction>
struct InterpolationLoopHelp<2, TGridFunction>{
static void invoke(IPData<number, TGridFunction::domain_type::dim>& InterpolFunction,
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
		try
		{
		const int dim = ssGrp.dim(i);
		switch(dim)
		{
		case DIM_SUBSET_EMPTY_GRID: break;
		case 0: /* \TODO: do nothing may be wrong */	break;
		case 1:
			InterpolateFunctionOnElem<Edge, TGridFunction>(InterpolFunction, u, fct, si, time);
			break;
		case 2:
			InterpolateFunctionOnElem<Triangle, TGridFunction>(InterpolFunction, u, fct, si, time);
			InterpolateFunctionOnElem<Quadrilateral, TGridFunction>(InterpolFunction, u, fct, si, time);
			break;
		default: UG_THROW_FATAL("InterpolateFunction: Dimension " << dim <<
		                        " not possible for world dim "<<2<<".");
		}
		}
		UG_CATCH_THROW("InterpolateFunction: Cannot interpolate on elements.");
	}

//	adjust parallel storage state
#ifdef UG_PARALLEL
	u.set_storage_type(PST_CONSISTENT);
#endif
}
};

template <typename TGridFunction>
struct InterpolationLoopHelp<1, TGridFunction>{
static void invoke(	IPData<number, TGridFunction::domain_type::dim>& InterpolFunction,
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
		try
		{
		const int dim = ssGrp.dim(i);
		switch(dim)
		{
		case DIM_SUBSET_EMPTY_GRID: break;
		case 0: /* \TODO: do nothing may be wrong */	break;
		case 1:
			InterpolateFunctionOnElem<Edge, TGridFunction>(InterpolFunction, u, fct, si, time);
			break;
		default: UG_THROW_FATAL("InterpolateFunction: Dimension " <<dim<<
		                        " not possible for world dim "<<1<<".");
		}
		}
		UG_CATCH_THROW("InterpolateFunction: Cannot interpolate on elements.");
	}

//	adjust parallel storage state
#ifdef UG_PARALLEL
	u.set_storage_type(PST_CONSISTENT);
#endif
}
};


/// interpolates a function on a subset
/**
 * This function interpolates a GridFunction. To evaluate the function on every
 * point a functor must be passed.
 *
 * \param[out]		u			interpolated grid function
 * \param[in]		data		data evaluator
 * \param[in]		name		symbolic name of function
 * \param[in]		time		time point
 * \param[in]		subsets		subsets, where to interpolate
 */
template <typename TGridFunction>
void InterpolateFunction(IPData<number, TGridFunction::domain_type::dim>& InterpolFunction,
                         TGridFunction& u, const char* name, number time, const char* subsets)
{
//	get function id of name
	const size_t fct = u.fct_id_by_name(name);

//	check that function found
	if(fct == (size_t)-1)
		UG_THROW_FATAL("InterpolateFunction: Name of function not found.");

//	check that function exists
	if(fct >= u.num_fct())
		UG_THROW_FATAL("InterpolateFunction: Function space does not contain"
						" a function with index " << fct);

//	create subset group
	SubsetGroup ssGrp; ssGrp.set_subset_handler(u.domain()->subset_handler());

//	read subsets
	if(subsets != NULL)
		ConvertStringToSubsetGroup(ssGrp, u.domain()->subset_handler(), subsets);
	else // add all if no subset specified
		ssGrp.add_all();


//	forward
	InterpolationLoopHelp<TGridFunction::dim, TGridFunction>::
								invoke(InterpolFunction, u, fct, time, ssGrp);
}

/// interpolates a function on the whole domain
template <typename TGridFunction>
void InterpolateFunction(IPData<number, TGridFunction::domain_type::dim>& InterpolFunction,
                         TGridFunction& u, const char* name, number time)
{
//	forward
	InterpolateFunction(InterpolFunction, u, name, time, NULL);
}

template <typename TGridFunction>
void InterpolateFunction(number val,
                         TGridFunction& u, const char* name, number time)
{
	ConstUserNumber<TGridFunction::domain_type::dim> p(val);
	InterpolateFunction(p, u, name, time);
}

#ifdef UG_FOR_LUA
template <typename TGridFunction>
void InterpolateFunction(const char* InterpolFunction,
                         TGridFunction& u, const char* name, number time)
{
	LuaUserData<number, TGridFunction::domain_type::dim> p(InterpolFunction);
	InterpolateFunction(p, u, name, time);
}

template <typename TGridFunction>
void InterpolateFunction(const char* InterpolFunction,
                         TGridFunction& u, const char* name, number time, const char* subsets)
{
	LuaUserData<number, TGridFunction::domain_type::dim> p(InterpolFunction);
	InterpolateFunction(p, u, name, time, subsets);
}
#endif


/////////////////////////////////////////////////////
// Vertices only section
/////////////////////////////////////////////////////

/// interpolates a function on vertices
template <typename TGridFunction>
void InterpolateFunctionOnVertices(IPData<number, TGridFunction::domain_type::dim>& InterpolFunction,
                                   TGridFunction& u, size_t fct, int si, number time)
{
//	domain type and position_type
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_type position_type;

//	id of shape functions used
	LFEID id = u.local_finite_element_id(fct);

// get position accessor
	const typename domain_type::position_accessor_type& aaPos
										= u.domain().position_accessor();

// 	iterate over all elements
	typename geometry_traits<VertexBase>::const_iterator iterEnd, iter;
	iterEnd = u.template end<VertexBase>(si);
	for(iter = u.template begin<VertexBase>(si); iter != iterEnd; ++iter)
	{
	//	get element
		VertexBase* vrt = *iter;

	//	global position
		position_type glob_pos = aaPos[vrt];

	//	value at position
		number val;
		(*InterpolFunction)(val, glob_pos, time);

	//	get multiindices of element
		typename TGridFunction::multi_index_vector_type ind;
		u.multi_indices(vrt, fct, ind);

	// 	loop all dofs
		for(size_t i = 0; i < ind.size(); ++i)
		{
		//	set value
			BlockRef(u[ind[i][0]], ind[i][1]) = val;
		}
	}
}



/// interpolates a function on subsets, only for vertex dofs
/**
 * This function interpolates a GridFunction. To evaluate the function on every
 * point a functor must be passed.
 *
 * \param[out]		u			interpolated grid function
 * \param[in]		data		data evaluator
 * \param[in]		name		symbolic name of function
 * \param[in]		time		time point
 * \param[in]		subsets		subsets, where to interpolate
 */
template <typename TGridFunction>
void InterpolateFunctionOnVertices(IPData<number, TGridFunction::domain_type::dim>& InterpolFunction,
                                   TGridFunction& u, const char* name, number time, const char* subsets)
{
//	get Function Pattern
	const typename TGridFunction::approximation_space_type& approxSpace
				= u.approximation_space();

//	get function id of name
	const size_t fct = approxSpace.fct_id_by_name(name);

//	check that function found
	if(fct == (size_t)-1)
		UG_THROW_FATAL("InterpolateFunctionOnVertices: Name of function not found.");

//	check that function exists
	if(fct >= u.num_fct())
		UG_THROW_FATAL("InterpolateFunctionOnVertices: Function space does not contain"
				" a function with index " << fct);

//	create subset group
	SubsetGroup ssGrp; ssGrp.set_subset_handler(*approxSpace.subset_handler());

//	read subsets
	if(subsets != NULL)
		ConvertStringToSubsetGroup(ssGrp, *approxSpace.subset_handler(), subsets);
	else // add all if no subset specified
		ssGrp.add_all();

//	loop subsets
	for(size_t i = 0; i < ssGrp.num_subsets(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	skip if function is not defined in subset
		if(!u.is_def_in_subset(fct, si)) continue;

		InterpolateFunctionOnVertices<TGridFunction>
							(InterpolFunction, u, fct, si, time);
	}

//	adjust parallel storage state
#ifdef UG_PARALLEL
	u.set_storage_type(PST_CONSISTENT);
#endif
}


} // namespace ug

#endif /*__H__UG__LIB_DISC__FUNCTION_SPACES__INTERPOLATE__*/
