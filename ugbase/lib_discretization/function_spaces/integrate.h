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


/// Abstract integrand interface (using Barton & Nackman trick)
template <class Derived, typename TElem, typename TGridFunction>
class Integrand {
private:
	Derived &asLeaf() {return static_cast<Derived&> (*this);}


public:
	/// initialize for element-wise evaluation
	bool initElem(TElem *elem)
	{return asLeaf().initElem(elem);}

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
	boost::function<void (number&, const MathVector<TGridFunction::domain_type::dim>&, number)> ExactSolution;
	number time;

	// aux. index array
	typename TGridFunction::multi_index_vector_type ind;

public:

	//typedef TElem element_type;
	/// constructor
	L2ErrorIntegrandTmp(TGridFunction &f, size_t fctid,
			boost::function<void (number&, const MathVector<TGridFunction::domain_type::dim>&, number)> e,
			number t) :
				//	grid function and id of shape functions used
				u(f), fct(fctid), id(u.local_finite_element_id(fct)),
				//	trial space
				trialSpace(LocalShapeFunctionSetProvider::get<typename reference_element_traits<TElem>::reference_element_type>(id)),
				// exact solution
				ExactSolution(e), time(t)
	{};


	/// initialize for element-wise evaluation
	bool initElem(TElem *elem) {

		//	get multiindices of element
		u.multi_indices(elem, fct, ind);

		//	check multi indices
		if(ind.size() != trialSpace.num_sh())
		{
			UG_LOG("ERROR in 'L2ErrorOnElem': Wrong number of"
					" multi indices.\n");
			return false;
		}

		return true;
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
		ExactSolution(exactSolIP, globIP, time);

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
bool SumValuesForSubsetGroup( number& addValue,
					boost::function<void (number& res,
									const MathVector<TGridFunction::domain_type::dim>& x,
									number time)> ExactSolution,
								TGridFunction& u, size_t fct, int si, number time)

/*
 template <typename Derived, typename TElem, typename TGridFunction>
 bool SumValuesForSubsetGroup( number& addValue,
					Integrand<Derived, TElem, TGridFunction> integrand, int si)



 */
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
	typename geometry_traits<TElem>::const_iterator iterEnd, iter;
	iterEnd = u.template end<TElem>(si);
	iter = u.template begin<TElem>(si);
	if(iter==iterEnd) return true;

	// initialize local class
	// NOTE: This is independent of the rest and may be provided as its own class!!!
	L2ErrorIntegrandTmp<TElem,TGridFunction> integrand= L2ErrorIntegrandTmp<TElem,TGridFunction>(u, fct, ExactSolution, time);

//	get quadrature Rule
	const QuadratureRule<dim>& rQuadRule
	= QuadratureRuleProvider<dim>::template get_rule<ref_elem_type>(order);

//	create a reference mapping
	ReferenceMapping<ref_elem_type, domain_type::dim> mapping;

//	id of shape functions used
//	LFEID id = u.local_finite_element_id(fct);

//	get trial space

//	typedef LocalShapeFunctionSet<ref_elem_type> trialspace_type;
	//const LocalShapeFunctionSet<ref_elem_type>& trialSpace =
		//	LocalShapeFunctionSetProvider::get<ref_elem_type>(id);
//	const trialspace_type& trialSpace =
//			LocalShapeFunctionSetProvider::get<ref_elem_type>(id);


//	number of dofs on element
//	const size_t num_sh = trialSpace.num_sh();

	// initialize local class
	// NOTE: This is independent of the rest and may be provided as its own class!!!
	//L2ErrorIntegrandTmp<TElem,TGridFunction> local = L2ErrorIntegrandTmp<TElem,TGridFunction>(u, fct, ExactSolution);


// 	iterate over all elements
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		TElem* elem = *iter;

	//	get all corner coordinates
		std::vector<position_type> vCorner;
		CollectCornerCoordinates(vCorner, *elem, u.get_domain());

	//	update the reference mapping for the corners
		mapping.update(&vCorner[0]);
	//	get multiindices of element
	//	typename TGridFunction::multi_index_vector_type ind;
	//	u.multi_indices(elem, fct, ind);

		//	check multi indices
		//		if(ind.size() != num_sh)
		//		{
		//			UG_LOG("ERROR in 'L2ErrorOnElem': Wrong number of"
		//					" multi indices.\n");
		//			return false;
		//		}

	// initialize integrand for a given element
		if (integrand.initElem(elem) ==false)
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
			intValElem += valIP * weightIP * det;
		}

	//	add to global sum
		addValue += intValElem;
	}

//	we're done
	return true;
}

template <int dim, typename TGridFunction>
struct IntegralElementLoopHelp
{};

template <typename TGridFunction>
struct IntegralElementLoopHelp<3, TGridFunction>{
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
			// L2ErrorIntegrandTmp<Edge, TGridFunction> integrand(InterpolFunction, u, time);
			// bRes &= SumValuesForSubsetGroup<Edge, TGridFunction>(diffSquared, integrand, si);
			bRes &= SumValuesForSubsetGroup<Edge, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			break;
		case 2:
			bRes &= SumValuesForSubsetGroup<Triangle, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			bRes &= SumValuesForSubsetGroup<Quadrilateral, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			break;
		case 3:
			bRes &= SumValuesForSubsetGroup<Tetrahedron, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			bRes &= SumValuesForSubsetGroup<Hexahedron, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			bRes &= SumValuesForSubsetGroup<Prism, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			bRes &= SumValuesForSubsetGroup<Pyramid, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
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
struct IntegralElementLoopHelp<1, TGridFunction>{
static number invoke(boost::function<void (number& res, const MathVector<TGridFunction::domain_type::dim>& x,number time)> InterpolFunction,
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
		bool bRes = true;
		switch(ssGrp.dim(i))
		{
		case 1:
			bRes &= SumValuesForSubsetGroup<Edge, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
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
struct IntegralElementLoopHelp<2, TGridFunction>{
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
			bRes &= SumValuesForSubsetGroup<Edge, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			break;
		case 2:
			bRes &= SumValuesForSubsetGroup<Triangle, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
			bRes &= SumValuesForSubsetGroup<Quadrilateral, TGridFunction>(diffSquared, InterpolFunction, u, fct, si, time);
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
	return IntegralElementLoopHelp<TGridFunction::dim, TGridFunction>::invoke(InterpolFunction, u, fct, time, ssGrp);
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
