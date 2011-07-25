/*
 * integrate.h
 *
 *  Created on: 04.04.2011
 *      Author: kxylouris, avogel
 */

#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACES__INTEGRATE__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACES__INTEGRATE__

#include <cmath>

#include "common/common.h"

#include "lib_discretization/common/subset_group.h"
#include "lib_discretization/domain_util.h"
#include "lib_discretization/quadrature/quadrature.h"
#include "lib_discretization/local_finite_element/local_shape_function_set.h"
#include <boost/function.hpp>

namespace ug{

/// interpolates a function on an element
/**
 * This function adds to diffValSquared the contributions of all elements on
 * the subset si.
 */
template <typename TElem, typename TGridFunction>
bool DiffSquaredOnElems( number& diffValSquared,
						boost::function<void (
									number& res,
									const MathVector<TGridFunction::domain_type::dim>& x,
									number time)> ExactSolution,
								TGridFunction& u, size_t fct, int si, number time)
{
//	order of quadrature rule
//	\todo: generalize
	const int order = 1;

//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type
				ref_elem_type;

//	dimension of reference element
	static const int dim = ref_elem_type::dim;

//	domain type and position_type
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_type position_type;

//	id of shape functions used
	LFEID id = u.local_finite_element_id(fct);

//	get trial space
	const LocalShapeFunctionSet<ref_elem_type>& trialSpace =
			LocalShapeFunctionSetProvider::get<ref_elem_type>(id);

//	number of dofs on element
	const size_t num_sh = trialSpace.num_sh();

//	get quadrature Rule
	const QuadratureRule<dim>& rQuadRule
			= QuadratureRuleProvider<dim>::template get_rule<ref_elem_type>(order);

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

	//	get multiindices of element
		typename TGridFunction::multi_index_vector_type ind;
		u.multi_indices(elem, fct, ind);

	//	check multi indices
		if(ind.size() != num_sh)
		{
			UG_LOG("ERROR in 'L2ErrorOnElem': Wrong number of"
					" multi indices.\n");
			return false;
		}

	//	contribution of this element
		number intValElem = 0;

	//	loop integration points
		for(size_t ip = 0; ip < rQuadRule.size(); ++ip)
		{
		//	get local integration point
			const MathVector<dim>& locIP = rQuadRule.point(ip);

		//	compute global integration point
			position_type globIP;

		//  map local dof position to global position
			mapping.local_to_global(globIP, locIP);

		//	compute exact solution at integration point
			number exactSolIP;
			ExactSolution(exactSolIP, globIP, time);

		//	reset approx solution
			number approxSolIP = 0;

		// 	sum up contributions of all shape functions at ip
			for(size_t sh = 0; sh < num_sh; ++sh)
			{
			//	get value at shape point (e.g. corner for P1 fct)
				const number valSH = BlockRef(u[ind[sh][0]], ind[sh][1]);

			//	add shape fct at ip * value at shape
				approxSolIP += valSH * trialSpace.shape(sh, locIP);
			}

		//	get squared of difference
			number diffVal = (exactSolIP - approxSolIP);
			diffVal *= diffVal;

		//	get quadrature weight
			const number weightIP = rQuadRule.weight(ip);

		//	get determinate of mapping
			const number det = mapping.jacobian_det(locIP);

		//	add contribution of integration point
			intValElem += diffVal * weightIP * det;
		}

	//	add to global sum
		diffValSquared += intValElem;
	}

//	we're done
	return true;
}

template <int dim, typename TGridFunction>
struct L2ErrorLoopHelp
{};

template <typename TGridFunction>
struct L2ErrorLoopHelp<3, TGridFunction>{
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
//	difference squared on all elements
	number diffSquared = 0;

//	loop subsets
	for(size_t i = 0; i < ssGrp.num_subsets(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	skip if function is not defined in subset
		if(!u.is_def_in_subset(fct, si)) continue;


	//	switch dimensions
		bool bRes = true;
		switch(ssGrp.dim(i))
		{
		case 1:
			bRes &= DiffSquaredOnElems<Edge, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			break;
		case 2:
			bRes &= DiffSquaredOnElems<Triangle, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			bRes &= DiffSquaredOnElems<Quadrilateral, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			break;
		case 3:
			bRes &= DiffSquaredOnElems<Tetrahedron, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			bRes &= DiffSquaredOnElems<Hexahedron, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			bRes &= DiffSquaredOnElems<Prism, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			bRes &= DiffSquaredOnElems<Pyramid, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			break;
		default: UG_LOG("ERROR in L2ErrorHelp: Dimension "<<ssGrp.dim(i) <<
		                " not supported in world dimension "<<3<<".");
			throw(UGFatalError("Dimension not supported."));
		}

	//	check success
		if(!bRes)
			throw(UGFatalError("Error when summing up l2norm"));
	}

//	compute norm by taking root
	const number l2norm = sqrt(diffSquared);

//	we're done
	return l2norm;
}
};

template <typename TGridFunction>
struct L2ErrorLoopHelp<1, TGridFunction>{
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
//	difference squared on all elements
	number diffSquared = 0;

//	loop subsets
	for(size_t i = 0; i < ssGrp.num_subsets(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	skip if function is not defined in subset
		if(!u.is_def_in_subset(fct, si)) continue;


	//	switch dimensions
		bool bRes = true;
		switch(ssGrp.dim(i))
		{
		case 1:
			bRes &= DiffSquaredOnElems<Edge, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			break;
		default: UG_LOG("ERROR in L2ErrorHelp: Dimension "<<ssGrp.dim(i) <<
		                " not supported in world dimension "<<1<<".");
			throw(UGFatalError("Dimension not supported."));
		}

	//	check success
		if(!bRes)
			throw(UGFatalError("Error when summing up l2norm"));
	}

//	compute norm by taking root
	const number l2norm = sqrt(diffSquared);

//	we're done
	return l2norm;
}
};

template <typename TGridFunction>
struct L2ErrorLoopHelp<2, TGridFunction>{
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
//	difference squared on all elements
	number diffSquared = 0;

//	loop subsets
	for(size_t i = 0; i < ssGrp.num_subsets(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	skip if function is not defined in subset
		if(!u.is_def_in_subset(fct, si)) continue;


	//	switch dimensions
		bool bRes = true;
		switch(ssGrp.dim(i))
		{
		case 1:
			bRes &= DiffSquaredOnElems<Edge, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			break;
		case 2:
			bRes &= DiffSquaredOnElems<Triangle, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			bRes &= DiffSquaredOnElems<Quadrilateral, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			break;
		default: UG_LOG("ERROR in L2ErrorHelp: Dimension "<<ssGrp.dim(i) <<
		                " not supported in world dimension "<<2<<".");
			throw(UGFatalError("Dimension not supported."));
		}

	//	check success
		if(!bRes)
			throw(UGFatalError("Error when summing up l2norm"));
	}

//	compute norm by taking root
	const number l2norm = sqrt(diffSquared);

//	we're done
	return l2norm;
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
number L2Error(
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
		UG_LOG("ERROR in L2Error: Name of function not found.\n");
		return false;
	}

//	check that function exists
	if(fct >= u.num_fct())
	{
		UG_LOG("ERROR in L2Error: Function space does not contain"
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
	return L2ErrorLoopHelp<TGridFunction::dim, TGridFunction>::invoke(InterpolFunction, u, fct, time, ssGrp);
}

/// interpolates a function on the whole domain
template <typename TGridFunction>
number L2Error(
	const boost::function<void (number& res, const MathVector<TGridFunction::domain_type::dim>& x, number time)>& InterpolFunction,
	TGridFunction& u, const char* name, number time)
{
//	forward
	return L2Error(InterpolFunction, u, name, time, NULL);
}

} // namespace ug

#endif /*__H__LIBDISCRETIZATION__FUNCTION_SPACES__INTEGRATE__*/
