/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Andreas Vogel, Christian Wehner
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/*
 *      an adaption of interpolate.h to maximum error computation
 */

#ifndef __H__UG__LIB_DISC__MAX__ERROR__
#define __H__UG__LIB_DISC__MAX__ERROR__

#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_grid/tools/subset_group.h"

#include "lib_disc/common/groups_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/reference_element/reference_mapping.h"
#include "lib_disc/domain_util.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Compute max error on Vertices only
////////////////////////////////////////////////////////////////////////////////
	
/**
 * This function computes max error of a grid function on a vertex loop. Thus, it can only
 * be used if all degrees of freedom are located in the vertices only (e.g. P1
 * finite elements). In those cases it is faster than the element by element
 * loop.
 *
 * @param[in] globalMaxError		reference to global max error
 * @param[in] spInterpolFunction	data providing exact solution
 * @param[out] spGridFct			grid function
 * @param[in] fct					symbolic name of function component
 * @param[in] ssGrp					subsets, where to compute max error
 * @param[in] time					time point
 */
template <typename TGridFunction>
void MaxErrorOnVertices(number& globalMaxError,SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                           SmartPtr<TGridFunction> spGridFct,
                           size_t fct,
                           number time,
                           const SubsetGroup& ssGrp)
{
//	domain type and position_type
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_type position_type;

// get position accessor
	const typename domain_type::position_accessor_type& aaPos
										= spGridFct->domain()->position_accessor();

	std::vector<DoFIndex> ind;
	typename TGridFunction::template dim_traits<0>::const_iterator iterEnd, iter;

	for(size_t i = 0; i < ssGrp.size(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	skip if function is not defined in subset
		if(!spGridFct->is_def_in_subset(fct, si)) continue;

	// 	iterate over all elements
		iterEnd = spGridFct->template end<Vertex>(si);
		iter = spGridFct->template begin<Vertex>(si);
		for(; iter != iterEnd; ++iter)
		{
		//	get element
			Vertex* vrt = *iter;

		//	global position
			position_type glob_pos = aaPos[vrt];

		//	value at position
			number val;
			(*spInterpolFunction)(val, glob_pos, time, si);

		//	get multiindices of element
			spGridFct->dof_indices(vrt, fct, ind);

		// 	loop all dofs
			for(size_t i = 0; i < ind.size(); ++i)
			{
			//	compute error
				number localError = std::abs(DoFRef((*spGridFct), ind[i])-val);
				if (localError>globalMaxError) globalMaxError=localError;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Maximum error on Elements
////////////////////////////////////////////////////////////////////////////////

/**
 * This function computes maximum error of a grid function on an element by element loop. 
 *
 * @param[in] globalMaxError		reference to global max error
 * @param[in] spInterpolFunction	data providing exact solution
 * @param[out] spGridFct			interpolated grid function
 * @param[in] fct					symbolic name of function component
 * @param[in] si					subset
 * @param[in] time					time point
 */
template <typename TElem, typename TGridFunction>
void MaxErrorOnElements(
		number& globalMaxError,
		SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
		SmartPtr<TGridFunction> spGridFct, size_t fct, int si, number time)
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				ref_elem_type;
	const ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

//	dimension of reference element
	const int dim = ref_elem_type::dim;

//	domain type and position_type
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_type position_type;

//	get iterators
	typename TGridFunction::template traits<TElem>::const_iterator iterEnd, iter;
	iterEnd = spGridFct->template end<TElem>(si);
	iter = spGridFct->template begin<TElem>(si);

//	check if something to do:
	if(iter == iterEnd) return;

//	id of shape functions used
	LFEID id = spGridFct->local_finite_element_id(fct);

//	get trial space
	const LocalShapeFunctionSet<dim>& trialSpace =
			LocalFiniteElementProvider::get<dim>(roid, id);

//	number of dofs on element
	const size_t nsh = trialSpace.num_sh();

// 	load local positions of dofs for the trial space on element
	std::vector<MathVector<dim> > loc_pos(nsh);
	for(size_t i = 0; i < nsh; ++i)
		if(!trialSpace.position(i, loc_pos[i]))
			UG_THROW("MaxErrorOnElem: Cannot find meaningful"
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
		CollectCornerCoordinates(vCorner, *elem, *spGridFct->domain());

	//	update the reference mapping for the corners
		mapping.update(&vCorner[0]);

	//	get multiindices of element
		std::vector<DoFIndex> ind;
		spGridFct->dof_indices(elem, fct, ind);

	//	check multi indices
		if(ind.size() != nsh)
			UG_THROW("MaxErrorOnElem: On subset "<<si<<": Number of shapes is "
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
			(*spInterpolFunction)(val, glob_pos, time, si);
			
		//	compute error
			number localError = std::abs(DoFRef((*spGridFct), ind[i])-val);
			if (localError>globalMaxError) globalMaxError=localError;
		}
	}
}

/**
 * This function computes maximum error of grid function on an element by element loop. 
 *
 * @param[in] globalMaxError		reference to global max error
 * @param[in] spInterpolFunction	data providing exact solution
 * @param[out] spGridFct			grid function
 * @param[in] fct					symbolic name of function component
 * @param[in] ssGrp					subsets
 * @param[in] time					time point
 */
template <typename TGridFunction>
void MaxErrorOnElements(
			number& globalMaxError,
			SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                           SmartPtr<TGridFunction> spGridFct,
                           size_t fct,
                           number time,
                           const SubsetGroup& ssGrp)
{
//	loop subsets
	for(size_t i = 0; i < ssGrp.size(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	skip if function is not defined in subset
		if(!spGridFct->is_def_in_subset(fct, si)) continue;

	//	switch dimensions
		try
		{
		const int dim = ssGrp.dim(i);

		if(dim > TGridFunction::dim)
			UG_THROW("MaxErrorOnElements: Dimension of subset is "<<dim<<", but "
			         " World Dimension is "<<TGridFunction::dim<<". Cannot interpolate this.");

		switch(dim)
		{
		case DIM_SUBSET_EMPTY_GRID: break;
		case 0: /* \TODO: do nothing may be wrong */	break;
		case 1:
			MaxErrorOnElements<RegularEdge, TGridFunction>(globalMaxError,spInterpolFunction, spGridFct, fct, si, time);
			break;
		case 2:
			MaxErrorOnElements<Triangle, TGridFunction>(globalMaxError,spInterpolFunction, spGridFct, fct, si, time);
			MaxErrorOnElements<Quadrilateral, TGridFunction>(globalMaxError,spInterpolFunction, spGridFct, fct, si, time);
			break;
		case 3:
			MaxErrorOnElements<Tetrahedron, TGridFunction>(globalMaxError,spInterpolFunction, spGridFct, fct, si, time);
			MaxErrorOnElements<Hexahedron, TGridFunction>(globalMaxError,spInterpolFunction, spGridFct, fct, si, time);
			MaxErrorOnElements<Prism, TGridFunction>(globalMaxError,spInterpolFunction, spGridFct, fct, si, time);
			MaxErrorOnElements<Pyramid, TGridFunction>(globalMaxError,spInterpolFunction, spGridFct, fct, si, time);
			MaxErrorOnElements<Octahedron, TGridFunction>(globalMaxError,spInterpolFunction, spGridFct, fct, si, time);
			break;
		default: UG_THROW("MaxErrorOnElements: Dimension " <<dim<<
		                " not possible for world dim "<<3<<".");
		}
		}
		UG_CATCH_THROW("MaxErrorOnElements: Cannot interpolate on elements.");
	}
}

////////////////////////////////////////////////////////////////////////////////
// MaxError routine
////////////////////////////////////////////////////////////////////////////////

/// computes maximum error of a grid function on a subset
/**
 * This function computes the maximum error of a GridFunction. To evaluate the exact solution on every
 * point a functor must be passed.
 *
 * @param[in] spInterpolFunction	data providing exact solution
 * @param[out] spGridFct			interpolated grid function
 * @param[in] cmp					symbolic name of function component
 * @param[in] subsets				subsets (NULL = everywhere)
 * @param[in] time					time point
 * @returns	  globalMaxError        maximum norm of difference
 */
template <typename TGridFunction>
number MaxError(
	SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets, number time)
{
	number globalMaxError=0;
//	check, that values do not depend on a solution
	if(spInterpolFunction->requires_grid_fct())
		UG_THROW("MaxError: The interpolation values depend on a grid function."
				" This is not allowed in the current implementation. Use constant,"
				" lua-callback or vrl-callback user data only (even within linkers).");

//	get function id of name
	const size_t fct = spGridFct->fct_id_by_name(cmp);

//	check that function found
	if(fct > spGridFct->num_fct())
		UG_THROW("MaxError: Name of component '"<<cmp<<"' not found.");

//	check if fast P1 interpolation can be used
	// \TODO: This should be improved. Manifold admissible if space continuous
	bool bUseP1Interpolation = false;
	if(spGridFct->local_finite_element_id(fct).type() == LFEID::LAGRANGE &&
			spGridFct->local_finite_element_id(fct).order() == 1)
		bUseP1Interpolation = true;
	const bool bAllowManyfoldInterpolation =
			(spGridFct->local_finite_element_id(fct).type() == LFEID::LAGRANGE);

//	create subset group
	SubsetGroup ssGrp(spGridFct->domain()->subset_handler());
	if(subsets != NULL)
	{
		ssGrp.add(TokenizeString(subsets));
		if(!bAllowManyfoldInterpolation)
			if(!SameDimensionsInAllSubsets(ssGrp))
				UG_THROW("MaxError: Subsets '"<<subsets<<"' do not have same dimension."
						 "Can not integrate on subsets of different dimensions.");
	}
	else
	{
	//	add all subsets and remove lower dim subsets afterwards
		ssGrp.add_all();
		if(!bAllowManyfoldInterpolation)
			RemoveLowerDimSubsets(ssGrp);
	}

//	forward
	if(bUseP1Interpolation)
		MaxErrorOnVertices<TGridFunction>(globalMaxError,spInterpolFunction, spGridFct, fct, time, ssGrp);
	else
		MaxErrorOnElements<TGridFunction>(globalMaxError,spInterpolFunction, spGridFct, fct, time, ssGrp);

	//	adjust parallel storage state
#ifdef UG_PARALLEL
	spGridFct->set_storage_type(PST_CONSISTENT);
#endif
	return globalMaxError;
}

template <typename TGridFunction>
number MaxError(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 number time)
{
	return MaxError(spInterpolFunction, spGridFct, cmp, NULL, time);
}
template <typename TGridFunction>
number MaxError(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets)
{
	return MaxError(spInterpolFunction, spGridFct, cmp, subsets, 0.0);
}
template <typename TGridFunction>
number MaxError(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp)
{
	return MaxError(spInterpolFunction, spGridFct, cmp, NULL, 0.0);
}

///////////////
// const data
///////////////

template <typename TGridFunction>
number MaxError(number val,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets, number time)
{
	static const int dim = TGridFunction::dim;
	SmartPtr<UserData<number, dim> > sp =
			make_sp(new ConstUserNumber<dim>(val));
	return MaxError(sp, spGridFct, cmp, subsets, time);
}
template <typename TGridFunction>
number MaxError(number val,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 number time)
{
	return MaxError(val, spGridFct, cmp, NULL, time);
}
template <typename TGridFunction>
number MaxError(number val,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets)
{
	return MaxError(val, spGridFct, cmp, subsets, 0.0);
}
template <typename TGridFunction>
number MaxError(number val,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp)
{
	return MaxError(val, spGridFct, cmp, NULL, 0.0);
}

///////////////
// lua data
///////////////

#ifdef UG_FOR_LUA
template <typename TGridFunction>
number MaxError(const char* LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets, number time)
{
	static const int dim = TGridFunction::dim;
	SmartPtr<UserData<number, dim> > sp =
			LuaUserDataFactory<number, dim>::create(LuaFunction);
	return MaxError(sp, spGridFct, cmp, subsets, time);
}
template <typename TGridFunction>
number MaxError(const char* LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 number time)
{	
	return MaxError(LuaFunction, spGridFct, cmp, NULL, time);
}

template <typename TGridFunction>
number MaxError(const char* LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets)
{	
	return MaxError(LuaFunction, spGridFct, cmp, subsets, 0.0);
}
template <typename TGridFunction>
number MaxError(const char* LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp)
{	
	return MaxError(LuaFunction, spGridFct, cmp, NULL, 0.0);
}
#endif


} // namespace ug

#endif /*__H__UG__LIB_DISC__MAX__ERROR__*/
