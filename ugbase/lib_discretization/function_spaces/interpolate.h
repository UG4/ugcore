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
#include "lib_discretization/local_finite_element/local_shape_function_set.h"
#include <boost/function.hpp>

namespace ug{

/// interpolates a function on an element
template <typename TElem, typename TGridFunction>
bool InterpolateFunctionOnElem(
	boost::function<void (number& res,const MathVector<TGridFunction::domain_type::dim>& x, number time)> InterpolFunction,
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
		{
			UG_LOG("ERROR in 'InterpolateFunctionOnElem': Cannot find meaningful"
					" local positions of dofs.\n");
			return false;
		}

//	create a reference mapping
	ReferenceMapping<ref_elem_type, domain_type::dim> mapping;

//	get iterators
	typename geometry_traits<TElem>::const_iterator iterEnd, iter;
	iterEnd = u.template end<TElem>(si);
	iter = u.template begin<TElem>(si);

// 	iterate over all elements
	for( ; iter != iterEnd; ++iter)
	{
	//	get element
		TElem* elem = *iter;

	//	get all corner coordinates
		std::vector<position_type> vCorner;
		CollectCornerCoordinates(vCorner, *elem, u.get_domain());

	//	update the reference mapping for the corners
		mapping.update(&vCorner[0]);

	//	get multiindices of element
		typename TGridFunction::multi_index_vector_type ind;
		u.multi_indices(elem, fct, ind);

	//	check multi indices
		if(ind.size() != nsh)
		{
			UG_LOG("ERROR in 'InterpolateFunctionOnElem': Number of shapes is "
					<<nsh<<", but got "<<ind.size()<<" multi indices.\n");
			return false;
		}

	// 	loop all dofs
		for(size_t i = 0; i < nsh; ++i)
		{
		//	global position
			position_type glob_pos;

		//  map local dof position to global position
			mapping.local_to_global(glob_pos, loc_pos[i]);

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

template <int dim, typename TGridFunction>
struct InterpolationLoopHelp
{};

template <typename TGridFunction>
struct InterpolationLoopHelp<3, TGridFunction>{
static number invoke(	boost::function<void (
									number& res,
									const MathVector<TGridFunction::domain_type::dim>& x,
									number time)>
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
#ifdef UG_PARALLEL
		const int dim = ssGrp.dim(i, &u.get_process_communicator());
#else
		const int dim = ssGrp.dim(i);
#endif
		switch(dim)
		{
		case DIM_SUBSET_EMPTY_GRID: break;
		case 0: /* \TODO: do nothing may be wrong */	break;
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
		default: UG_LOG("InterpolateFunction: Dimension " <<dim<<
		                " not possible for world dim "<<3<<".\n"); return false;
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
};

template <typename TGridFunction>
struct InterpolationLoopHelp<2, TGridFunction>{
static number invoke(	boost::function<void (
									number& res,
									const MathVector<TGridFunction::domain_type::dim>& x,
									number time)>
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
#ifdef UG_PARALLEL
		const int dim = ssGrp.dim(i, &u.get_process_communicator());
#else
		const int dim = ssGrp.dim(i);
#endif
		switch(dim)
		{
		case DIM_SUBSET_EMPTY_GRID: break;
		case 0: /* \TODO: do nothing may be wrong */	break;
		case 1:
			bRes &= InterpolateFunctionOnElem<Edge, TGridFunction>(InterpolFunction, u, fct, si, time);
			break;
		case 2:
			bRes &= InterpolateFunctionOnElem<Triangle, TGridFunction>(InterpolFunction, u, fct, si, time);
			bRes &= InterpolateFunctionOnElem<Quadrilateral, TGridFunction>(InterpolFunction, u, fct, si, time);
			break;
		default: UG_LOG("InterpolateFunction: Dimension " << dim <<
		                " not possible for world dim "<<2<<".\n"); return false;
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
};

template <typename TGridFunction>
struct InterpolationLoopHelp<1, TGridFunction>{
static number invoke(	boost::function<void (
									number& res,
									const MathVector<TGridFunction::domain_type::dim>& x,
									number time)>
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
#ifdef UG_PARALLEL
		const int dim = ssGrp.dim(i, &u.get_process_communicator());
#else
		const int dim = ssGrp.dim(i);
#endif
		switch(dim)
		{
		case DIM_SUBSET_EMPTY_GRID: break;
		case 0: /* \TODO: do nothing may be wrong */	break;
		case 1:
			bRes &= InterpolateFunctionOnElem<Edge, TGridFunction>(InterpolFunction, u, fct, si, time);
			break;
		default: UG_LOG("InterpolateFunction: Dimension " <<dim<<
		                " not possible for world dim "<<1<<".\n"); return false;
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
bool InterpolateFunction(
		const boost::function<void (number& res, const MathVector<TGridFunction::domain_type::dim>& x, number time)>& InterpolFunction,
		TGridFunction& u, const char* name, number time, const char* subsets)
{
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
	return InterpolationLoopHelp<TGridFunction::dim, TGridFunction>::
								invoke(InterpolFunction, u, fct, time, ssGrp);
}

/// interpolates a function on the whole domain
template <typename TGridFunction>
bool InterpolateFunction(
		const boost::function<void (number& res, const MathVector<TGridFunction::domain_type::dim>& x, number time)>& InterpolFunction,
		TGridFunction& u, const char* name, number time)
{
//	forward
	return InterpolateFunction(InterpolFunction, u, name, time, NULL);
}


/// copy a gridfunction only on a subset
template <typename TGridFunction>
bool AssignP1GridFunctionOnSubset(TGridFunction& uDest, const TGridFunction& uSrc,
                                const char* subsetNames)
{
//	get Approximation Space
	const typename TGridFunction::approximation_space_type& approxSpace
												= uDest.get_approximation_space();
//	create subset group
	SubsetGroup ssGrp;
	ssGrp.set_subset_handler(*approxSpace.get_subset_handler());

//	read subsets
	ConvertStringToSubsetGroup(ssGrp, *approxSpace.get_subset_handler(), subsetNames);


//	loop subsets
	for(size_t i = 0; i < ssGrp.num_subsets(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];


	//	loop vertices
		typename geometry_traits<VertexBase>::const_iterator iter
									= uDest.template begin<VertexBase>(si);
		typename geometry_traits<VertexBase>::const_iterator iterEnd
									= uDest.template end<VertexBase>(si);

		for(; iter != iterEnd; ++iter)
		{
		//	get vertex
			VertexBase* vrt = *iter;

		//	get algebra indices
			typename TGridFunction::algebra_index_vector_type vInd;
			uDest.inner_algebra_indices(vrt, vInd);

		//	loop indices
			for(size_t k = 0; k < vInd.size(); ++k)
			{
			//	get algebra index
				const size_t index = vInd[k];

			//	assign values
				uDest[index] = uSrc[index];

			}
		}
	}

//	we're done
	return true;
}

/////////////////////////////////////////////////////
// Vertices only section
/////////////////////////////////////////////////////

/// interpolates a function on vertices
template <typename TGridFunction>
bool InterpolateFunctionOnVertices(
		boost::function<void (number& res, const MathVector<TGridFunction::domain_type::dim>& x, number time)> InterpolFunction,
		TGridFunction& u, size_t fct, int si, number time)
{
//	domain type and position_type
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_type position_type;

//	id of shape functions used
	LFEID id = u.local_finite_element_id(fct);

// get position accessor
	const typename domain_type::position_accessor_type& aaPos
										= u.get_domain().get_position_accessor();

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
		InterpolFunction(val, glob_pos, time);

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

//	we're done
	return true;
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
bool InterpolateFunctionOnVertices(
		const boost::function<void (number& res, const MathVector<TGridFunction::domain_type::dim>& x, number time)>& InterpolFunction,
		TGridFunction& u, const char* name, number time, const char* subsets)
{
//	get Function Pattern
	const typename TGridFunction::approximation_space_type& approxSpace
				= u.get_approximation_space();

//	get function id of name
	const size_t fct = approxSpace.fct_id_by_name(name);

//	check that function found
	if(fct == (size_t)-1)
	{
		UG_LOG("ERROR in InterpolateFunctionOnVertices: Name of function not found.\n");
		return false;
	}

//	check that function exists
	if(fct >= u.num_fct())
	{
		UG_LOG("ERROR in InterpolateFunctionOnVertices: Function space does not contain"
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

//	loop subsets
	for(size_t i = 0; i < ssGrp.num_subsets(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	skip if function is not defined in subset
		if(!u.is_def_in_subset(fct, si)) continue;

		if(!InterpolateFunctionOnVertices<TGridFunction>
							(InterpolFunction, u, fct, si, time))
		{
			UG_LOG("ERROR in InterpolateFunctionOnVertices: "
					"Cannot interpolate on elements.\n");
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


} // namespace ug

#endif /*__H__LIBDISCRETIZATION__FUNCTION_SPACES__INTERPOLATE__*/
