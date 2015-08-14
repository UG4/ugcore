/*
 * interpolate_inner.h
 *
 *  Created on: 14.08.2015
 *      Author: mbreit
 */

#ifndef UG__LIB_DISC__FUNCTION_SPACES__INTERPOLATE_INNER_H
#define UG__LIB_DISC__FUNCTION_SPACES__INTERPOLATE_INNER_H

#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_disc/function_spaces/interpolate.h"
#include "lib_disc/domain_util.h"
#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/reference_element/reference_mapping.h"
#include "lib_disc/function_spaces/dof_position_util.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{

/**
 * This function interpolates a grid function on a specific subset and for
 * a specific element type. All inner dofs of all grid elements of the given
 * type in the given subset will be assigned an interpolated value.
 *
 * @param[in] spInterpolFunction	data providing interpolation values
 * @param[out] spGridFct			interpolated grid function
 * @param[in] fct					symbolic name of function component
 * @param[in] si					subset, where to interpolate
 * @param[in] time					time point
 */
template <typename TElem, typename TGridFunction>
void InterpolateOnElementsInner(
		SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
		SmartPtr<TGridFunction> spGridFct, size_t fct, int si, number time)
{
	// do nothing if no DoFs here
	if (!spGridFct->max_fct_dofs(fct, ROID_TETRAHEDRON, si))
		return;

	//	domain type and position_type
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_type position_type;

	//	get iterators
	typename TGridFunction::template traits<TElem>::const_iterator iterEnd, iter;
	iterEnd = spGridFct->template end<TElem>(si);
	iter = spGridFct->template begin<TElem>(si);

	//	check if something to do:
	if (iter == iterEnd) return;

	//	id of shape functions used
	LFEID id = spGridFct->local_finite_element_id(fct);

	// 	iterate over all elements
	for (; iter != iterEnd; ++iter)
	{
	//	get element
		TElem* elem = *iter;

	//	get multiindices of element
		std::vector<DoFIndex> ind;
		spGridFct->inner_dof_indices(elem, fct, ind);

	//	global positions of DoFs
		std::vector<position_type> glob_pos;
		InnerDoFPosition<domain_type>(glob_pos, elem, *spGridFct->domain(), id);

	//	check global positions size
		size_t ind_sz = ind.size();
		size_t gp_sz = glob_pos.size();
		if (ind_sz != gp_sz)
			UG_THROW("InterpolateOnElem: On subset " << si << ": Number of DoFs is "
					<< ind_sz << ", but number of DoF positions is "
					<< gp_sz << "." << std::endl);

	// 	loop all dofs
		for (size_t i = 0; i < ind_sz && i < gp_sz; ++i)
		{
		//	value at position
			number val;
			(*spInterpolFunction)(val, glob_pos[i], time, si);

		//	set value
			DoFRef(*spGridFct, ind[i]) = val;
		}
	}
}


/**
 * This function interpolates a grid function on a subset group. All inner
 * dofs of all grid elements in the subset group will be assigned an inter-
 * polated value.
 *
 * @param[in] spInterpolFunction	data providing interpolation values
 * @param[out] spGridFct			interpolated grid function
 * @param[in] fct					symbolic name of function component
 * @param[in] ssGrp					subsets, where to interpolate
 * @param[in] time					time point
 */
template <typename TGridFunction>
void InterpolateOnElementsInner
(
	SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
	SmartPtr<TGridFunction> spGridFct,
	size_t fct,
	number time,
	const SubsetGroup& ssGrp
)
{
//	loop subsets
	for(size_t i = 0; i < ssGrp.size(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	skip if function is not defined in subset
		if (!spGridFct->is_def_in_subset(fct, si)) continue;

	//	switch dimensions
		try
		{
			const int dim = ssGrp.dim(i);

			if (dim > TGridFunction::dim)
				UG_THROW("InterpolateOnElements: Dimension of subset is " << dim << ", but "
						 " World Dimension is " << TGridFunction::dim << ". Cannot interpolate this.");

			// In a parallel scenario, the distribution CAN cause elements of of lower
			// dimension than the rest of their subset to be located disconnected from
			// the rest of the subset on a processor. For example, in 2D, think of a
			// (1D) boundary subset and a distribution where the boundary of a proc's
			// domain only touches the boundary subset in a vertex, but intersects with
			// the boundary subset in another place.
			// Therefore, we need to interpolate on all dimensions < dim and we will
			// only consider inner indices then, thus we simply switch-case backwards
			// and without breaks.
			switch (dim)
			{
				case DIM_SUBSET_EMPTY_GRID: break;
				case 3:
					if (spGridFct->max_fct_dofs(fct, VOLUME, si))
					{
						InterpolateOnElementsInner<Tetrahedron, TGridFunction>
							(spInterpolFunction, spGridFct, fct, si, time);
						InterpolateOnElementsInner<Hexahedron, TGridFunction>
							(spInterpolFunction, spGridFct, fct, si, time);
						InterpolateOnElementsInner<Prism, TGridFunction>
							(spInterpolFunction, spGridFct, fct, si, time);
						InterpolateOnElementsInner<Pyramid, TGridFunction>
							(spInterpolFunction, spGridFct, fct, si, time);
						InterpolateOnElementsInner<Octahedron, TGridFunction>
							(spInterpolFunction, spGridFct, fct, si, time);
					}
				case 2:
					if (spGridFct->max_fct_dofs(fct, FACE, si))
					{
						InterpolateOnElementsInner<Triangle, TGridFunction>
							(spInterpolFunction, spGridFct, fct, si, time);
						InterpolateOnElementsInner<Quadrilateral, TGridFunction>
							(spInterpolFunction, spGridFct, fct, si, time);
					}
				case 1:
					if (spGridFct->max_fct_dofs(fct, EDGE, si))
						InterpolateOnElementsInner<RegularEdge, TGridFunction>
							(spInterpolFunction, spGridFct, fct, si, time);
				case 0:
					if (spGridFct->max_fct_dofs(fct, VERTEX, si))
						InterpolateOnElementsInner<RegularVertex, TGridFunction>
							(spInterpolFunction, spGridFct, fct, si, time);
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
 * This function interpolates a Lagrange type GridFunction.
 * In contrast to its "big brother" Interpolate() in interpolate.h,
 * this function will only write to _inner_ DoFs, but to _all_ of them,
 * on any element located in any of the given subsets.
 * This is more intuitive and secure than writing to all (not necessarily
 * inner) DoFs of an element located in a subset, as this takes into
 * account the unique subsets of sub-elements.
 *
 * At the moment, this is only meant for Lagrange-type elements.
 * Please feel free to add functionality for other types if you like.
 *
 * @param[in] spInterpolFunction	data providing interpolation values
 * @param[out] spGridFct			interpolated grid function
 * @param[in] cmp					symbolic name of function component
 * @param[in] subsets				subsets, where to interpolate (NULL = everywhere)
 * @param[in] time					time point
 */
template <typename TGridFunction>
void InterpolateInner(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets, number time)
{
	// check, that values do not depend on a solution
	if (spInterpolFunction->requires_grid_fct())
		UG_THROW("Interpolate: The interpolation values depend on a grid function."
				" This is not allowed in the current implementation. Use constant,"
				" lua-callback or vrl-callback user data only (even within linkers).");

	// get function id of name
	const size_t fct = spGridFct->fct_id_by_name(cmp);

	// check that function found
	if (fct > spGridFct->num_fct())
		UG_THROW("Interpolate: Name of component '"<< cmp <<"' not found.");

	// check that type is Lagrange
	if (spGridFct->local_finite_element_id(fct).type() != LFEID::LAGRANGE)
		UG_THROW("This interpolation only allows Lagrange-type elements.\n"
				 "Feel free to add functionality for other types as needed.");

	// check if fast P1 interpolation can be used
	const bool bUseP1Interpolation
		= spGridFct->local_finite_element_id(fct).order() == 1;

	// create subset group
	SubsetGroup ssGrp(spGridFct->domain()->subset_handler());
	if (subsets != NULL)
		ssGrp.add(TokenizeString(subsets));
	else
		ssGrp.add_all();

	// forward
	if (bUseP1Interpolation)
		InterpolateOnVertices<TGridFunction>(spInterpolFunction, spGridFct, fct, time, ssGrp);
	else
		InterpolateOnElementsInner<TGridFunction>(spInterpolFunction, spGridFct, fct, time, ssGrp);

	//	adjust parallel storage state
#ifdef UG_PARALLEL
	spGridFct->set_storage_type(PST_CONSISTENT);
#endif
}

template <typename TGridFunction>
void InterpolateInner(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 number time)
{InterpolateInner(spInterpolFunction, spGridFct, cmp, NULL, time);}
template <typename TGridFunction>
void InterpolateInner(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets)
{InterpolateInner(spInterpolFunction, spGridFct, cmp, subsets, 0.0);}
template <typename TGridFunction>
void InterpolateInner(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp)
{InterpolateInner(spInterpolFunction, spGridFct, cmp, NULL, 0.0);}

///////////////
// const data
///////////////

template <typename TGridFunction>
void InterpolateInner(number val,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets, number time)
{
	static const int dim = TGridFunction::dim;
	SmartPtr<UserData<number, dim> > sp =
			make_sp(new ConstUserNumber<dim>(val));
	InterpolateInner(sp, spGridFct, cmp, subsets, time);
}
template <typename TGridFunction>
void InterpolateInner(number val,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 number time)
{InterpolateInner(val, spGridFct, cmp, NULL, time);}
template <typename TGridFunction>
void InterpolateInner(number val,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets)
{InterpolateInner(val, spGridFct, cmp, subsets, 0.0);}
template <typename TGridFunction>
void InterpolateInner(number val,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp)
{InterpolateInner(val, spGridFct, cmp, NULL, 0.0);}

///////////////
// lua data
///////////////

#ifdef UG_FOR_LUA
template <typename TGridFunction>
void InterpolateInner(const char* LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets, number time)
{
	static const int dim = TGridFunction::dim;
	SmartPtr<UserData<number, dim> > sp =
			LuaUserDataFactory<number, dim>::create(LuaFunction);
	InterpolateInner(sp, spGridFct, cmp, subsets, time);
}
template <typename TGridFunction>
void InterpolateInner(const char* LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 number time)
{InterpolateInner(LuaFunction, spGridFct, cmp, NULL, time);}
template <typename TGridFunction>
void InterpolateInner(const char* LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets)
{InterpolateInner(LuaFunction, spGridFct, cmp, subsets, 0.0);}
template <typename TGridFunction>
void InterpolateInner(const char* LuaFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp)
{InterpolateInner(LuaFunction, spGridFct, cmp, NULL, 0.0);}
#endif


} // namespace ug

#endif // UG__LIB_DISC__FUNCTION_SPACES__INTERPOLATE_INNER_H
