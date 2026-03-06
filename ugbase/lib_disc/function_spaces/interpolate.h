/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACES__INTERPOLATE__
#define __H__UG__LIB_DISC__FUNCTION_SPACES__INTERPOLATE__

#include "grid_function_global_user_data.h"
#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_grid/tools/subset_group.h"

#include "lib_disc/domain_util.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/reference_element/reference_mapping.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Interpolate on Vertices only
////////////////////////////////////////////////////////////////////////////////

/**
 * This function interpolates a grid function on a vertex loop. Thus, it can only
 * be used if all degrees of freedom are located in the vertices only (e.g. P1
 * finite elements). In those cases it is faster than the element by element
 * loop.
 *
 * @param[in] spInterpolFunction	data providing interpolation values
 * @param[out] spGridFct			interpolated grid function
 * @param[in] fct					symbolic name of function component
 * @param[in] ssGrp					subsets, where to interpolate
 * @param[in] time					time point
 * @param[in] diff_pos				different vector between Interpolatin values and spGridFct.
 */

template <typename TGridFunction>
void InterpolateOnDiffVertices(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
        SmartPtr<TGridFunction> spGridFct,
        size_t fct,
        number time,
        const SubsetGroup& ssGrp,
		const MathVector<TGridFunction::dim> diff_pos)
{

	//	domain type and position_type
		typedef typename TGridFunction::domain_type domain_type;
		typedef typename domain_type::position_type position_type;

		//std::cout<<"Interpolate diff_vector: "<<diff_pos;

	// get position accessor (of interpolated grid function)
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
				position_type glob_pos = aaPos[vrt]; // position (of interpolated grid function)
				position_type rel_pos=glob_pos;
				rel_pos -=diff_pos;

			//	value at position
				number val;
				(*spInterpolFunction)(val, rel_pos, time, si);

			//	get multiindices of element
				spGridFct->dof_indices(vrt, fct, ind);

			// 	loop all dofs
				for(size_t i = 0; i < ind.size(); ++i)
				{
				//	set value
					DoFRef(*spGridFct, ind[i]) = val;
				}
			}
		}
}
//getting value of spInterpolFunction at position


template <typename TGridFunction>
number get_number_on_coords(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
	typename TGridFunction::domain_type::position_type pos,
	number time,
    const int si
){
	number val;
	(*spInterpolFunction)(val, pos, time, si);

	return val;
}

/**
 * This function interpolates a grid function on a vertex loop. Thus, it can only
 * be used if all degrees of freedom are located in the vertices only (e.g. P1
 * finite elements). In those cases it is faster than the element by element
 * loop.
 *
 * @param[in] spInterpolFunction	data providing interpolation values
 * @param[out] spGridFct			interpolated grid function
 * @param[in] fct					symbolic name of function component
 * @param[in] ssGrp					subsets, where to interpolate
 * @param[in] time					time point
 */
template <typename TGridFunction>
void InterpolateOnVertices(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                           SmartPtr<TGridFunction> spGridFct,
                           size_t fct,
                           number time,
                           const SubsetGroup& ssGrp)
{
	//	dimension of reference element
	const int dim = TGridFunction::dim;

	MathVector<dim> diff_pos(0.0);
	InterpolateOnDiffVertices<TGridFunction>(spInterpolFunction, spGridFct, fct, time, ssGrp, diff_pos);
}



////////////////////////////////////////////////////////////////////////////////
// Interpolate on Elements
////////////////////////////////////////////////////////////////////////////////

/**
 * This function interpolates a grid function on an element by element loop. On
 * each element the all associated (up to the boundary of the element) are
 * interpolated and the values are stored in the grid function.
 *
 * @param[in] spInterpolFunction	data providing interpolation values
 * @param[out] spGridFct			interpolated grid function
 * @param[in] fct					symbolic name of function component
 * @param[in] si					subset, where to interpolate
 * @param[in] time					time point
 */
template <typename TElem, typename TGridFunction>
void InterpolateOnDiffElements(
		SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
		SmartPtr<TGridFunction> spGridFct, size_t fct, int si, number time,
		const MathVector<TGridFunction::dim> diff_pos)
{
//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
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
			UG_THROW("InterpolateOnElem: Cannot find meaningful"
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
			UG_THROW("InterpolateOnElem: On subset "<<si<<": Number of shapes is "
					<<nsh<<", but got "<<ind.size()<<" multi indices.");

	// 	loop all dofs
		for(size_t i = 0; i < nsh; ++i)
		{
		//	global position
			position_type glob_pos;


		//  map local dof position to global position
			mapping.local_to_global(glob_pos, loc_pos[i]);

			position_type rel_pos=glob_pos;
			rel_pos -=diff_pos;

		//	value at position
			number val;
			(*spInterpolFunction)(val, rel_pos, time, si);

		//	set value
			DoFRef(*spGridFct, ind[i]) = val;
		}
	}
}

/**
 * This function interpolates a grid function on an element by element loop. On
 * each element the all associated (up to the boundary of the element) are
 * interpolated and the values are stored in the grid function.
 *
 * @param[in] spInterpolFunction	data providing interpolation values
 * @param[out] spGridFct			interpolated grid function
 * @param[in] fct					symbolic name of function component
 * @param[in] ssGrp					subsets, where to interpolate
 * @param[in] time					time point
 */
template <typename TGridFunction>
void InterpolateOnDiffElements(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                           SmartPtr<TGridFunction> spGridFct,
                           size_t fct,
                           number time,
                           const SubsetGroup& ssGrp,const MathVector<TGridFunction::dim> diff_pos)
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
			UG_THROW("InterpolateOnElements: Dimension of subset is "<<dim<<", but "
			         " World Dimension is "<<TGridFunction::dim<<". Cannot interpolate this.");

		// FIXME (at least for Lagrange, order > 1, parallel)
		// In a parallel scenario, the distribution CAN cause elements of of lower
		// dimension than the rest of their subset to be located disconnected from
		// the rest of the subset on a processor. For example, in 2D, think of a
		// (1D) boundary subset and a distribution where the boundary of a proc's
		// domain only touches the boundary subset in a vertex, but intersects with
		// the boundary subset in another place.
		// This vertex will not be considered during interpolation even though it
		// should be!
		switch(dim)
		{
		case DIM_SUBSET_EMPTY_GRID: break;
		case 0: /* \TODO: do nothing may be wrong */	break;
		case 1:
			InterpolateOnDiffElements<RegularEdge, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time,diff_pos);
			break;
		case 2:
			InterpolateOnDiffElements<Triangle, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time, diff_pos);
			InterpolateOnDiffElements<Quadrilateral, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time, diff_pos);
			break;
		case 3:
			InterpolateOnDiffElements<Tetrahedron, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time, diff_pos);
			InterpolateOnDiffElements<Hexahedron, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time, diff_pos);
			InterpolateOnDiffElements<Prism, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time, diff_pos);
			InterpolateOnDiffElements<Pyramid, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time, diff_pos);
			InterpolateOnDiffElements<Octahedron, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time, diff_pos);
			break;
		default: UG_THROW("InterpolateOnElements: Dimension " <<dim<<
		                " not possible for world dim "<<3<<".");
		}
		}
		UG_CATCH_THROW("InterpolateOnElements: Cannot interpolate on elements.");
	}
}

/**
 * This function interpolates a grid function on an element by element loop. On
 * each element the all associated (up to the boundary of the element) are
 * interpolated and the values are stored in the grid function.
 *
 * @param[in] spInterpolFunction	data providing interpolation values
 * @param[out] spGridFct			interpolated grid function
 * @param[in] fct					symbolic name of function component
 * @param[in] si					subset, where to interpolate
 * @param[in] time					time point
 */
template <typename TElem, typename TGridFunction>
void InterpolateOnElements(
		SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
		SmartPtr<TGridFunction> spGridFct, size_t fct, int si, number time)
{
	//	dimension of reference element
	const int dim = TGridFunction::dim;

	MathVector<dim> diff_pos(0.0);
	InterpolateOnDiffElements<TElem,TGridFunction>(spInterpolFunction, spGridFct, fct, si,time, diff_pos);
}

/**
 * This function interpolates a grid function on an element by element loop. On
 * each element the all associated (up to the boundary of the element) are
 * interpolated and the values are stored in the grid function.
 *
 * @param[in] spInterpolFunction	data providing interpolation values
 * @param[out] spGridFct			interpolated grid function
 * @param[in] fct					symbolic name of function component
 * @param[in] ssGrp					subsets, where to interpolate
 * @param[in] time					time point
 */
template <typename TGridFunction>
void InterpolateOnElements(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
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
			UG_THROW("InterpolateOnElements: Dimension of subset is "<<dim<<", but "
			         " World Dimension is "<<TGridFunction::dim<<". Cannot interpolate this.");

		// FIXME (at least for Lagrange, order > 1, parallel)
		// In a parallel scenario, the distribution CAN cause elements of of lower
		// dimension than the rest of their subset to be located disconnected from
		// the rest of the subset on a processor. For example, in 2D, think of a
		// (1D) boundary subset and a distribution where the boundary of a proc's
		// domain only touches the boundary subset in a vertex, but intersects with
		// the boundary subset in another place.
		// This vertex will not be considered during interpolation even though it
		// should be!
		switch(dim)
		{
		case DIM_SUBSET_EMPTY_GRID: break;
		case 0: /* \TODO: do nothing may be wrong */	break;
		case 1:
			InterpolateOnElements<RegularEdge, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time);
			break;
		case 2:
			InterpolateOnElements<Triangle, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time);
			InterpolateOnElements<Quadrilateral, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time);
			break;
		case 3:
			InterpolateOnElements<Tetrahedron, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time);
			InterpolateOnElements<Hexahedron, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time);
			InterpolateOnElements<Prism, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time);
			InterpolateOnElements<Pyramid, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time);
			InterpolateOnElements<Octahedron, TGridFunction>(spInterpolFunction, spGridFct, fct, si, time);
			break;
		default: UG_THROW("InterpolateOnElements: Dimension " <<dim<<
		                " not possible for world dim "<<3<<".");
		}
		}
		UG_CATCH_THROW("InterpolateOnElements: Cannot interpolate on elements.");
	}
}

////////////////////////////////////////////////////////////////////////////////
// Interpolate routine
////////////////////////////////////////////////////////////////////////////////

/// interpolates a function on a subset
/**
 * This function interpolates a GridFunction. To evaluate the function on every
 * point a functor must be passed.
 *
 * @param[in] spInterpolFunction	data providing interpolation values
 * @param[out] spGridFct			interpolated grid function
 * @param[in] cmp					symbolic name of function component
 * @param[in] subsets				subsets, where to interpolate (NULL = everywhere)
 * @param[in] time					time point
 */
template <typename TGridFunction>
void Interpolate(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets, number time)
{
	//	dimension of reference element
	const int dim = TGridFunction::dim;
	MathVector<dim> diff_pos(0.0);
	Interpolate(spInterpolFunction, spGridFct, cmp, subsets, time, diff_pos);
}

/// interpolates a function on a subset
/**
 * This function interpolates a GridFunction. To evaluate the function on every
 * point a functor must be passed.
 *
 * @param[in] spInterpolFunction	data providing interpolation values
 * @param[out] spGridFct			interpolated grid function
 * @param[in] fct					id of function component
 * @param[in] subsets				subsets, where to interpolate (NULL = everywhere)
 * @param[in] time					time point
 */
template <typename TGridFunction>
void Interpolate(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const size_t fct,
				 const SubsetGroup& ssGrp, number time, const MathVector<TGridFunction::dim> diff_pos)
{

	//	check, that values do not depend on a solution
	if(spInterpolFunction->requires_grid_fct())
		UG_THROW("Interpolate: The interpolation values depend on a grid function."
				" This is not allowed in the current implementation. Use constant,"
				" lua-callback or vrl-callback user data only (even within linkers).");

//	check if fast P1 interpolation can be used
	// \TODO: This should be improved. Manifold admissible if space continuous
	bool bUseP1Interpolation = false;
	if(spGridFct->local_finite_element_id(fct).type() == LFEID::LAGRANGE &&
			spGridFct->local_finite_element_id(fct).order() == 1)
		bUseP1Interpolation = true;

	//forward
	if(bUseP1Interpolation){
		InterpolateOnDiffVertices<TGridFunction>(spInterpolFunction, spGridFct, fct, time, ssGrp, diff_pos);
	}else{
		InterpolateOnDiffElements<TGridFunction>(spInterpolFunction, spGridFct, fct, time, ssGrp, diff_pos);
	}
	//	adjust parallel storage state
#ifdef UG_PARALLEL
	spGridFct->set_storage_type(PST_CONSISTENT);
#endif
}

/// interpolates a function on a subset
/**
 * This function interpolates a GridFunction. To evaluate the function on every
 * point a functor must be passed.
 *
 * @param[in] spInterpolFunction	data providing interpolation values
 * @param[out] spGridFct			interpolated grid function
 * @param[in] cmp					symbolic name of function component
 * @param[in] subsets				subsets, where to interpolate (NULL = everywhere)
 * @param[in] time					time point
 */
template <typename TGridFunction>
void Interpolate(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets, number time, const MathVector<TGridFunction::dim> diff_pos)
{


//	get function id of name
	const size_t fct = spGridFct->fct_id_by_name(cmp);

//	check that function found
	if(fct > spGridFct->num_fct())
		UG_THROW("Interpolate: Name of component '"<<cmp<<"' not found.");

	Interpolate(spInterpolFunction, spGridFct,  fct, subsets, time, diff_pos);
}


/// interpolates a function on a subset
/**
 * This function interpolates a GridFunction. To evaluate the function on every
 * point a functor must be passed.
 *
 * @param[in] spInterpolFunction	data providing interpolation values
 * @param[out] spGridFct			interpolated grid function
 * @param[in] cmp					symbolic name of function component
 * @param[in] subsets				subsets, where to interpolate (NULL = everywhere)
 * @param[in] time					time point
 */
template <typename TGridFunction>
void Interpolate(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const size_t fct,
                 const char* subsets, number time, const MathVector<TGridFunction::dim> diff_pos)
{
	const bool bAllowManyfoldInterpolation =
				(spGridFct->local_finite_element_id(fct).type() == LFEID::LAGRANGE);

	//	create subset group
	SubsetGroup ssGrp(spGridFct->domain()->subset_handler());
	if(subsets != NULL)
	{
		ssGrp.add(TokenizeString(subsets));
		if(!bAllowManyfoldInterpolation)
			if(!SameDimensionsInAllSubsets(ssGrp))
				UG_THROW("Interpolate: Subsets '"<<subsets<<"' do not have same dimension."
							"Can not integrate on subsets of different dimensions.");
	}
	else
	{
	//	add all subsets and remove lower dim subsets afterwards
		ssGrp.add_all();
		if(!bAllowManyfoldInterpolation)
			RemoveLowerDimSubsets(ssGrp);
	}

	Interpolate(spInterpolFunction, spGridFct,  fct, ssGrp, time, diff_pos);
}


template <typename TGridFunction>
void Interpolate(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 number time)
{
	Interpolate(spInterpolFunction, spGridFct, cmp, NULL, time);
}
template <typename TGridFunction>
void Interpolate(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets)
{
	Interpolate(spInterpolFunction, spGridFct, cmp, subsets, 0.0);
}
template <typename TGridFunction>
void Interpolate(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp)
{
	Interpolate(spInterpolFunction, spGridFct, cmp, NULL, 0.0);
}
//interpolate with diff_vector
template <typename TGridFunction>
void Interpolate(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp, const MathVector<TGridFunction::dim>& m_diff_pos)
{
	Interpolate(spInterpolFunction, spGridFct, cmp, NULL, 0.0, m_diff_pos);
}

///////////////
// const data
///////////////

template <typename TGridFunction>
void Interpolate(number val,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets, number time)
{
	static const int dim = TGridFunction::dim;
	SmartPtr<UserData<number, dim> > sp =
			make_sp(new ConstUserNumber<dim>(val));
	Interpolate(sp, spGridFct, cmp, subsets, time);
}
template <typename TGridFunction>
void Interpolate(number val,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 number time)
{Interpolate(val, spGridFct, cmp, NULL, time);}
template <typename TGridFunction>
void Interpolate(number val,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets)
{Interpolate(val, spGridFct, cmp, subsets, 0.0);}
template <typename TGridFunction>
void Interpolate(number val,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp)
{Interpolate(val, spGridFct, cmp, NULL, 0.0);}

//interpolate with diff_vector
template <typename TGridFunction>
void Interpolate(number val,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets, number time,const SmartPtr<CplUserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > m_diff_pos)
{
	static const int dim = TGridFunction::dim;
	SmartPtr<UserData<number, dim> > sp =
			make_sp(new ConstUserNumber<dim>(val));

	InterpolateDiff(sp, spGridFct, cmp, subsets, time,m_diff_pos);
}



///////////////
// lua data
///////////////

#ifdef UG_FOR_LUA
// function-name as string
template <typename TGridFunction>
void Interpolate(const char* LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets, number time)
{
	static const int dim = TGridFunction::dim;
	SmartPtr<UserData<number, dim> > sp =
			LuaUserDataFactory<number, dim>::create(LuaFunction);
	Interpolate(sp, spGridFct, cmp, subsets, time);
}
template <typename TGridFunction>
void Interpolate(const char* LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 number time)
{Interpolate(LuaFunction, spGridFct, cmp, NULL, time);}
template <typename TGridFunction>
void Interpolate(const char* LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets)
{Interpolate(LuaFunction, spGridFct, cmp, subsets, 0.0);}
template <typename TGridFunction>
void Interpolate(const char* LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp)
{Interpolate(LuaFunction, spGridFct, cmp, NULL, 0.0);}

// function as function handle
template <typename TGridFunction>
void Interpolate(LuaFunctionHandle LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets, number time)
{
	static const int dim = TGridFunction::dim;
	SmartPtr<UserData<number, dim> > sp =
			make_sp(new LuaUserData<number, dim>(LuaFunction));
	Interpolate(sp, spGridFct, cmp, subsets, time);
}
template <typename TGridFunction>
void Interpolate(LuaFunctionHandle LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 number time)
{Interpolate(LuaFunction, spGridFct, cmp, NULL, time);}
template <typename TGridFunction>
void Interpolate(LuaFunctionHandle LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets)
{Interpolate(LuaFunction, spGridFct, cmp, subsets, 0.0);}
template <typename TGridFunction>
void Interpolate(LuaFunctionHandle LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp)
{Interpolate(LuaFunction, spGridFct, cmp, NULL, 0.0);}

//interpolate with Diff-vector:
template <typename TGridFunction>
void InterpolateDiff(const char* LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets, number time, SmartPtr<CplUserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > m_diff_pos )
{
	static const int dim = TGridFunction::dim;
	SmartPtr<UserData<number, dim> > sp =
			LuaUserDataFactory<number, dim>::create(LuaFunction);
	InterpolateDiff(sp, spGridFct, cmp,subsets, time, m_diff_pos);
}

template <typename TGridFunction>
void InterpolateDiff(LuaFunctionHandle LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets, number time,SmartPtr<CplUserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > m_diff_pos)
{
	static const int dim = TGridFunction::dim;
	SmartPtr<UserData<number, dim> > sp =
			make_sp(new LuaUserData<number, dim>(LuaFunction));
	InterpolateDiff(sp, spGridFct, cmp, subsets, time, m_diff_pos);
}
template <typename TGridFunction>
void Interpolate(LuaFunctionHandle LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,const SmartPtr<CplUserData<MathVector<TGridFunction::dim>, TGridFunction::dim> > m_diff_pos)
{InterpolateDiff(LuaFunction, spGridFct, cmp, NULL, 0.0, m_diff_pos);}



#endif
#ifdef UG_PARALLEL
	/// interpolates a gridFunction across similiar geometries with different meshing, where distribution across processes is not
	/// consistent. This is inefficient and should not be used in hot loops.
	/**
	 * Creates GlobalGridFunctionNumber and transfers unknown values across different vertex configurations on different processes
	 *
	 * @param[in] spGfSource			data providing interpolation values
	 * @param[out] spGfTarget			interpolated grid function
	 * @param[in] fct					symbolic name of function component to be interpolated
	 * @param[in] time					time point
	 */
	template <typename TGridFunction>
	void InterpolateGlobalGridFunctionAcrossProcesses(SmartPtr<TGridFunction> spGfSource,
			SmartPtr<TGridFunction> spGfTarget,
			const char* cmp,
			number time)
{

	static const int dim  = TGridFunction::dim;
	SmartPtr<GlobalGridFunctionNumberData<TGridFunction> > data = make_sp(new GlobalGridFunctionNumberData<TGridFunction>(spGfSource, cmp)) ;
	//	domain type and position_type

	// check if we have multiple processes - if no we can just forward globalGfUserData to standard interpolate
	if (pcl::NumProcs() == 1) {
		Interpolate(data, spGfTarget, cmp, time);
		return;
	}
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_type position_type;


	// get position accessor (of interpolated grid function)
	const typename domain_type::position_accessor_type& aaPos
										= spGfTarget->domain()->position_accessor();


	std::vector<DoFIndex> ind;
	typename TGridFunction::template dim_traits<0>::const_iterator iterEnd, iter;
	const size_t fct = spGfSource->fct_id_by_name(cmp);


	pcl::ProcessCommunicator com;
	// each process iterates over its own local elements
	for (int i = 0; i < pcl::NumProcs(); ++i){
		// 	iterate over all elements
		position_type glob_pos;
		bool finished = false;
		// check for active Proc
		if (pcl::ProcRank() == i) {
			iterEnd = spGfTarget->template end<Vertex>();
			iter = spGfTarget->template begin<Vertex>();
			for(; iter != iterEnd; )
			{
				//	get vertex
				Vertex* vrt = *iter;

				//	global position (in case this position is not contained on the same proc for both distributions, we have to
				//	search all procs)
				glob_pos = aaPos[vrt]; // position (of interpolated grid function)
				// current active proc sends current evaluation position
				com.broadcast(glob_pos, i);


				//	value at position
				number val;
				// evaluate source grid function value for cmp and given position
				// evaluate_global() is already mpi - parallelized and thus needs to be called with the same
				// poosition from all procs - see logic for non active proc below
				data->evaluate_global(val,  glob_pos);

				// set local function value
				//	get multiindices of element
				spGfTarget->dof_indices(vrt, fct, ind);

				// 	loop all dofs
				for(size_t i = 0; i < ind.size(); ++i)
				{
					//	set value
					DoFRef(*spGfTarget, ind[i]) = val;
				}
				++iter;
				// check if active proc is finished with element loop and broadcast result
				finished = (iter ==	iterEnd);
				com.broadcast(finished, i);
			}


		}
		// if we are not the active proc, we wait for the active proc to finish iterating over its local elements
		// and answer the evaluate_global function call with the corresponding broadcasted position
		else {
			while (!finished){
				// while active proc has elements, we recieve current position and evaluate
				com.broadcast(glob_pos, i);

				number val;
				// matching evaluate_global call for active Proc
				data->evaluate_global(val,  glob_pos);
				com.broadcast(finished, i);
				// finish loop when active proc has no elements left
			}
		}


	}
	// resulting vector is consistent
	spGfTarget->set_storage_type(PST_CONSISTENT);

}
// implement as overload for Interpolate
template <typename TGridFunction>
void Interpolate(SmartPtr<TGridFunction> spGfSource,
					SmartPtr<TGridFunction> spGfTarget,
					const char* cmp,
					number time)
{InterpolateGlobalGridFunctionAcrossProcesses(spGfSource, spGfTarget,cmp, time);}
#endif

} // namespace ug

#endif /*__H__UG__LIB_DISC__FUNCTION_SPACES__INTERPOLATE__*/
