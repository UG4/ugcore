/*
 * integrate.h
 *
 *  Created on: 04.04.2011
 *      Author: kxylouris, avogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACES__INTEGRATE__
#define __H__UG__LIB_DISC__FUNCTION_SPACES__INTEGRATE__

#include <cmath>

#include "common/common.h"

#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "lib_disc/spatial_disc/ip_data/ip_data.h"
#include <boost/function.hpp>

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{


/// Abstract integrand interface (using Barton & Nackman trick)
template <class Derived, typename TElem, typename TGridFunction>
class Integrand {
private:
	Derived &asLeaf() {return static_cast<Derived&> (*this);}


public:
	/// initialize for element-wise evaluation
	void initElem(TElem *elem)
	{asLeaf().initElem(elem);}

	/// provide values at IP
	void getValue(number &diffVal,
				const MathVector<reference_element_traits<TElem>::reference_element_type::dim>& locIP,
				typename TGridFunction::domain_type::position_type& globIP)
	{ asLeaf().getValue(diffVal, locIP, globIP);}
};


/// This class provides the integrand for the L2 error of a grid function.
/** For each element, it returns the L2 error in the integration points.*/
template <typename TElem, typename TGridFunction>
class L2ErrorIntegrandTmp : Integrand< L2ErrorIntegrandTmp<TElem, TGridFunction>, TElem, TGridFunction>{

private:
	// grid function
	TGridFunction &u;
	size_t fct;
	LFEID id;
	const LocalShapeFunctionSet<typename reference_element_traits<TElem>::reference_element_type>  &trialSpace;

	// exact solution
	IPData<number, TGridFunction::domain_type::dim>& ExactSolution;
	number time; //<time
	int si; //< subset index

	// aux. index array
	typename std::vector<MultiIndex<2> > ind;

public:

	//typedef TElem element_type;
	/// constructor
	L2ErrorIntegrandTmp(TGridFunction &f, size_t fctid,
	            		IPData<number, TGridFunction::domain_type::dim>& e,
	            		number t, int _si) :
				//	grid function and id of shape functions used
				u(f), fct(fctid), id(u.local_finite_element_id(fct)),
				//	trial space
				trialSpace(LocalShapeFunctionSetProvider::get<typename reference_element_traits<TElem>::reference_element_type>(id)),
				// exact solution
				ExactSolution(e), time(t), si(_si)
	{};


	/// initialize for element-wise evaluation
	void initElem(TElem *elem) {

		//	get multiindices of element
		u.multi_indices(elem, fct, ind);

		//	check multi indices
		if(ind.size() != trialSpace.num_sh())
			UG_THROW("L2ErrorOnElem: Wrong number of multi indices.");
	};

	/// provide values at IP
	void getValue(number &diffVal,
			const MathVector<reference_element_traits<TElem>::reference_element_type::dim>& locIP,
			typename TGridFunction::domain_type::position_type& globIP)
	{


		//	number of dofs on element
		const size_t num_sh = trialSpace.num_sh();

		//	compute exact solution at integration point
		number exactSolIP;
		ExactSolution(exactSolIP, globIP, time, si);

		//	reset approx solution
		number approxSolIP = 0.0;

		// 	sum up contributions of all shape functions at ip
		for(size_t sh = 0; sh < num_sh; ++sh)
		{
			//	get value at shape point (e.g. corner for P1 fct)
			const number valSH = BlockRef(u[ind[sh][0]], ind[sh][1]);

			//	add shape fct at ip * value at shape
			approxSolIP += valSH * trialSpace.shape(sh, locIP);
		}

		//	get squared of difference
		diffVal = (exactSolIP - approxSolIP);
		diffVal *= diffVal;
	};
};


/// Integrate some function over a subset group si
/**
 * @param
 * @param si: index for subset group*/
template <typename TElem, typename TGridFunction>
void SumValuesForSubsetGroup(number& addValue,
                             IPData<number, TGridFunction::domain_type::dim>& ExactSolution,
                             TGridFunction& u, size_t fct, int si, number time)
{
//	order of quadrature rule
//	\todo: generalize
	const int order = 1;

//	get reference element type
	typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

//	dimension of reference element
	static const int dim = ref_elem_type::dim;

//	domain type and position_type
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_type position_type;

//	check if something to do
	typename TGridFunction::template traits<TElem>::const_iterator iterEnd, iter;
	iterEnd = u.template end<TElem>(si);
	iter = u.template begin<TElem>(si);
	if(iter==iterEnd) return;

	// initialize local class
	// NOTE: This is independent of the rest and may be provided as its own class!!!
	L2ErrorIntegrandTmp<TElem,TGridFunction> integrand =
			L2ErrorIntegrandTmp<TElem,TGridFunction>(u, fct, ExactSolution, time, si);

//	get quadrature Rule
	const QuadratureRule<dim>& rQuadRule
	= QuadratureRuleProvider<dim>::template get_rule<ref_elem_type>(order);

//	create a reference mapping
	ReferenceMapping<ref_elem_type, domain_type::dim> mapping;

// 	iterate over all elements
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		TElem* elem = *iter;

	//	get all corner coordinates
		std::vector<position_type> vCorner;
		CollectCornerCoordinates(vCorner, *elem, *u.domain());

	//	update the reference mapping for the corners
		mapping.update(&vCorner[0]);

	// initialize integrand for a given element
		integrand.initElem(elem);

	//	contribution of this element
		number intValElem = 0;

	//	loop integration points
		for(size_t ip = 0; ip < rQuadRule.size(); ++ip)
		{

		//	get local integration point and compute global coordinates
			const MathVector<dim>& locIP = rQuadRule.point(ip);
			position_type globIP;
			mapping.local_to_global(globIP, locIP);

		//	compute value in integration point
			number valIP;
			integrand.getValue(valIP, locIP, globIP);

		//	get quadrature weight
			const number weightIP = rQuadRule.weight(ip);

		//	get determinate of mapping
			const number det = mapping.jacobian_det(locIP);

		//	add contribution of integration point
			intValElem += valIP * weightIP * fabs(det);
		}

	//	add to global sum
		addValue += intValElem;
	}
}

template <int dim, typename TGridFunction>
struct IntegralElementLoopHelp
{};

template <typename TGridFunction>
struct IntegralElementLoopHelp<3, TGridFunction>{
static number invoke(IPData<number, TGridFunction::domain_type::dim>& InterpolFunction,
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
		try
		{
		switch(ssGrp.dim(i))
		{
		case 1:
			SumValuesForSubsetGroup<Edge, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			break;
		case 2:
			SumValuesForSubsetGroup<Triangle, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			SumValuesForSubsetGroup<Quadrilateral, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			break;
		case 3:
			SumValuesForSubsetGroup<Tetrahedron, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			SumValuesForSubsetGroup<Hexahedron, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			SumValuesForSubsetGroup<Prism, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			SumValuesForSubsetGroup<Pyramid, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			break;
		default: UG_THROW("L2ErrorHelp: Dimension "<<ssGrp.dim(i) <<
		                        " not supported in world dimension "<<3<<".");
		}
		}
		UG_CATCH_THROW("Error when summing up l2norm");
	}

//	compute norm by taking root
	const number l2norm = sqrt(diffSquared);

//	we're done
	return l2norm;
}
};

template <typename TGridFunction>
struct IntegralElementLoopHelp<1, TGridFunction>{
static number invoke(IPData<number, TGridFunction::domain_type::dim>& InterpolFunction,
                     TGridFunction& u, size_t fct, number time, const SubsetGroup& ssGrp)
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
		try
		{
		switch(ssGrp.dim(i))
		{
		case 1:
			SumValuesForSubsetGroup<Edge, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			break;
		default: UG_THROW("L2ErrorHelp: Dimension "<<ssGrp.dim(i) <<
		                        " not supported in world dimension "<<1<<".");
		}
		}
		UG_CATCH_THROW("Error when summing up l2norm");
	}

//	compute norm by taking root
	const number l2norm = sqrt(diffSquared);

//	we're done
	return l2norm;
}
};

template <typename TGridFunction>
struct IntegralElementLoopHelp<2, TGridFunction>{
static number invoke(IPData<number, TGridFunction::domain_type::dim>& InterpolFunction,
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
		try
		{
		switch(ssGrp.dim(i))
		{
		case 1:
			SumValuesForSubsetGroup<Edge, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			break;
		case 2:
			SumValuesForSubsetGroup<Triangle, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			SumValuesForSubsetGroup<Quadrilateral, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			break;
		default: UG_THROW("L2ErrorHelp: Dimension "<<ssGrp.dim(i) <<
		                " not supported in world dimension "<<2<<".");
		}
		}
		UG_CATCH_THROW("Error when summing up l2norm");
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
		IPData<number, TGridFunction::domain_type::dim>& InterpolFunction,
		TGridFunction& u, const char* name, number time, const char* subsets)
{
//	get function id of name
	const size_t fct = u.fct_id_by_name(name);

//	check that function found
	if(fct == (size_t)-1) UG_THROW("ERROR in L2Error: Name of function not found.");

//	check that function exists
	if(fct >= u.num_fct())
		UG_THROW("ERROR in L2Error: Function space does not contain"
				" a function with index " << fct);

//	create subset group
	SubsetGroup ssGrp; ssGrp.set_subset_handler(u.domain()->subset_handler());

//	read subsets
	if(subsets != NULL)
		ConvertStringToSubsetGroup(ssGrp, u.domain()->subset_handler(), subsets);
	else // add all if no subset specified
		ssGrp.add_all();


//	forward
	return IntegralElementLoopHelp<TGridFunction::dim, TGridFunction>::invoke(InterpolFunction, u, fct, time, ssGrp);
}

/// interpolates a function on the whole domain
template <typename TGridFunction>
number L2Error(
	IPData<number, TGridFunction::domain_type::dim>& InterpolFunction,
	TGridFunction& u, const char* name, number time)
{
//	forward
	return L2Error(InterpolFunction, u, name, time, NULL);
}

#ifdef UG_FOR_LUA
template <typename TGridFunction>
number L2Error(const char* InterpolFunction,
             TGridFunction& u, const char* name, number time)
{
	LuaUserData<number, TGridFunction::domain_type::dim> p(InterpolFunction);
	return L2Error(p, u, name, time);
}

template <typename TGridFunction>
number L2Error(const char* InterpolFunction,
             TGridFunction& u, const char* name, number time, const char* subsets)
{
	LuaUserData<number, TGridFunction::domain_type::dim> p(InterpolFunction);
	return L2Error(p, u, name, time, subsets);
}
#endif


} // namespace ug

#endif /*__H__UG__LIB_DISC__FUNCTION_SPACES__INTEGRATE__*/
