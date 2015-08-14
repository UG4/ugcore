/*
 * interpolate.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACES__INTERPOLATE__
#define __H__UG__LIB_DISC__FUNCTION_SPACES__INTERPOLATE__

#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_disc/domain_util.h"
#include "lib_disc/common/subset_group.h"
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
 */
template <typename TGridFunction>
void InterpolateOnVertices(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
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
			//	set value
				DoFRef(*spGridFct, ind[i]) = val;
			}
		}
	}
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
void InterpolateOnElements(
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

		//	value at position
			number val;
			(*spInterpolFunction)(val, glob_pos, time, si);

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
//	check, that values do not depend on a solution
	if(spInterpolFunction->requires_grid_fct())
		UG_THROW("Interpolate: The interpolation values depend on a grid function."
				" This is not allowed in the current implementation. Use constant,"
				" lua-callback or vrl-callback user data only (even within linkers).");

//	get function id of name
	const size_t fct = spGridFct->fct_id_by_name(cmp);

//	check that function found
	if(fct > spGridFct->num_fct())
		UG_THROW("Interpolate: Name of component '"<<cmp<<"' not found.");

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

//	forward
	if(bUseP1Interpolation)
		InterpolateOnVertices<TGridFunction>(spInterpolFunction, spGridFct, fct, time, ssGrp);
	else
		InterpolateOnElements<TGridFunction>(spInterpolFunction, spGridFct, fct, time, ssGrp);

	//	adjust parallel storage state
#ifdef UG_PARALLEL
	spGridFct->set_storage_type(PST_CONSISTENT);
#endif
}

template <typename TGridFunction>
void Interpolate(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 number time)
{Interpolate(spInterpolFunction, spGridFct, cmp, NULL, time);}
template <typename TGridFunction>
void Interpolate(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets)
{Interpolate(spInterpolFunction, spGridFct, cmp, subsets, 0.0);}
template <typename TGridFunction>
void Interpolate(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp)
{Interpolate(spInterpolFunction, spGridFct, cmp, NULL, 0.0);}

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

///////////////
// lua data
///////////////

#ifdef UG_FOR_LUA
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
#endif


} // namespace ug

#endif /*__H__UG__LIB_DISC__FUNCTION_SPACES__INTERPOLATE__*/
