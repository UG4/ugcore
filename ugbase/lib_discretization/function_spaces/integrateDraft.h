/*
 * integrate.h
 *
 *  Created on: 04.04.2011
 *      Author: kxylouris, avogel
 */

#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACES__INTEGRATEDRAFT__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACES__INTEGRATEDRAFT__

#include <cmath>

#include "common/common.h"

#include "lib_discretization/common/subset_group.h"
#include "lib_discretization/domain_util.h"
#include "lib_discretization/quadrature/quadrature.h"
#include "lib_discretization/local_finite_element/local_shape_function_set.h"
#include <boost/function.hpp>

namespace ug{


/// Abstract integrand interface
template <int dim>
class IIntegrand
{
	public:
	///	returns the values of the integrand for a bunch of ips
	/**
	 * \param[in]	pElem	the element to integrate
	 * \param[out]	value	the value of the integrand at the ips
	 * \param[in]	vLocIP	local integration point positions
	 * \param[in]	vGlobIP	global integration point positions
	 * \param[in]	numIP	number of integration points
	 */
		virtual void getValues(GeometricObject* pElem,
		                       number value[],
		                       const MathVector<dim> vLocIP[],
		                       const MathVector<dim> vGlobIP[][],
		                       size_t numIP) = 0;
};


template <typename TGridFunction>
class L2ErrorIntegrand : IIntegrand<TGridFunction::dim>
{
	public:
	//	world dimension of grid function
		static const int dim = TGridFunction::dim;

	private:
	// grid function
		TGridFunction* pGridFct;

	//	component of function
		size_t fct;

	//	local finite element id
		LFEID id;

	//	local shape function set
		const DimLocalShapeFunctionSet<dim>* pTrialSpace;

	// exact solution
		typedef boost::function<void (number&, const MathVector<dim>&, number)> ExactSolFunctor;
		ExactSolFunctor ExactSolution;

	//	time
		number time;

	// 	aux. index array
		typename TGridFunction::multi_index_vector_type ind;

	public:

	/// constructor
		L2ErrorIntegrand() : pGridFct(NULL) {};

	///	sets grid function and component
		void set_grid_function(TGridFunction& gridFct, size_t cmp)
		{
			pGridFct = &gridFct;
			fct = cmp;
			id = pGridFct->local_finite_element_id(fct);
		}

	///	set exact solution
		void set_exact_solution(ExactSolFunctor& exactSol)
		{
			ExactSolution = exactSol;
		}

	///	set time
		void set_time(number time_) {time = time_;}

	/// \copydoc IIntegrand::getValues
		virtual void getValues(GeometricObject* pElem,
		                       number value[],
		                       const MathVector<dim> vLocIP[],
		                       const MathVector<dim> vGlobIP[][],
		                       size_t numIP)
		{
		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			ReferenceObjectID roid = (ReferenceObjectID) pElem->reference_object_id();

			try{
		//	get trial space
			const DimLocalShapeFunctionSet<dim>& TrialSpace =
							LocalShapeFunctionSetProvider::get<dim>(roid, id);

		//	number of dofs on element
			const size_t num_sh = trialSpace.num_sh();

		//	get multiindices of element
			u.multi_indices(pElem, fct, ind);

		//	check multi indices
			if(ind.size() != num_sh)
			{
				UG_LOG("ERROR in 'L2ErrorIntegrand::getValues': Wrong number of"
						" multi indices.\n");
				return false;
			}

		//	loop all integration points
			for(size_t ip = 0; ip < numIP; ++ip)
			{
			//	compute exact solution at integration point
				number exactSolIP;
				ExactSolution(exactSolIP, vGlobIP[ip], time);

			// 	compute approximated solution at integration point
				number approxSolIP = 0.0;
				for(size_t sh = 0; sh < num_sh; ++sh)
				{
				//	get value at shape point (e.g. corner for P1 fct)
					const number valSH = BlockRef(u[ind[sh][0]], ind[sh][1]);

				//	add shape fct at ip * value at shape
					approxSolIP += valSH * trialSpace.shape(sh, vLocIP[ip]);
				}

			//	get squared of difference
				value[ip] = (exactSolIP - approxSolIP);
				value[ip] *= value[ip];
			}

			}catch(UG_ERROR_LocalShapeFunctionSetNotRegistered& ex)
			{
				UG_LOG("ERROR in 'L2ErrorIntegrand::getValues': "<<ex.get_msg()<<"\n");
				return false;
			}
		};
};


/// integrates on the whole domain
/**
 * This function integrates an arbitrary integrand over the whole domain.
 * Note:
 *  - only grid elements of the same dimension as the world dimension of the
 *    domain are integrated. Thus, no manifolds.
 *  - The implementation is using virtual functions. Thus, there is a small
 *    performance drawback compared to hard coding everything, but we gain
 *    flexibility. In addition all virtual calls compute for the whole set of
 *    integration points to avoid many virtual calls, i.e. only one virtual
 *    call for all integration points is needed.
 *
 *    \todo: the selection of the elements should not be done using a grid
 *    		 function u, but instead a domain should be passed.
 *
 * \param[in]		u			grid function
 * \param[in]		integrand	Integrand
 * \param[in]		quadOrder	order of quadrature rule
 */
template <typename TGridFunction>
number Integrate(TGridFunction& u, IIntegrand<TGridFunction::dim>& integrand, int quadOrder)
{
//	reset the result
	number integral = 0;

//	get world (and reference) dimension
	static const int dim = TGridFunction::dim;

//	domain type and position_type
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_type position_type;

//	note: these iterators are for the base elements, e.g. Face and not
//			for the special type, e.g. Triangle, Quadrilateral
	typename domain_traits<dim>::const_iterator iterEnd, iter;

//	this is the base element type (e.g. Face). This is the type when the
//	iterators above are dereferenciated.
	typedef typename domain_traits<dim>::geometric_base_object geometric_base_object;

//	loop all subsets
	for(int si = 0; si < u.num_subsets(); ++si)
	{
	//	get iterators
		iterEnd = u.template end<TElem>(si);
		iter = u.template begin<TElem>(si);

	// 	iterate over all elements
		for(; iter != iterEnd; ++iter)
		{
		//	get element
			geometric_base_object* pElem = *iter;

		//	get reference object id (i.e. Triangle, Quadrilateral, Tetrahedron, ...)
			ReferenceObjectID roid = (ReferenceObjectID) pElem->reference_object_id();

		//	get quadrature Rule for reference object id and order
			try{
			const QuadratureRule<dim>& quadRule
					= QuadratureRuleProvider<dim>::get_rule(roid, orderQuad);

		//	number of integration points
			const size_t numIP = rQuadRule.size();

		//	get reference element mapping by reference object id
			try{
			DimReferenceMapping<dim, dim>& mapping
								= ReferenceMappingProvider::get<dim, dim>(roid);

		//	get all corner coordinates
			std::vector<position_type> vCorner;
			CollectCornerCoordinates(vCorner, *pElem, u.get_domain());

		//	update the reference mapping for the corners
			mapping.update(&vCorner[0]);

		//	compute global integration points
			std::vector<position_type> vGlobIP(numIP);
			mapping.local_to_global(&(vGlobIP[0]), rQuadRule.points(), numIP);

		//	compute integrand values at integration points
			std::vector<number> vValue(numIP);
			integrand.getValue(pElem, &(vValue[0]),
			                   rQuadRule.points(), &(vGlobIP[0]), numIP);

		//	reset contribution of this element
			number intValElem = 0;

		//	loop integration points
			for(size_t ip = 0; ip < numIP; ++ip)
			{
			//	get quadrature weight
				const number weightIP = rQuadRule.weight(ip);

			//	get determinate of mapping
				const number det = mapping.jacobian_det(locIP);

			//	add contribution of integration point
				intValElem += vValue[ip] * weightIP * det;
			}

		//	add to global sum
			integral += intValElem;

			}catch(UG_ERROR_QuadratureRuleNotRegistered& ex){
				UG_LOG("ERROR in 'SumValuesOnElems': " << ex.get_msg() << ".\n");
				return false;
			}
			}catch(UG_ERROR_ReferenceMappingMissing& ex){
				UG_LOG("ERROR in 'SumValuesOnElems': " << ex.get_msg() << ".\n");
				return false;
			}
		} // end elem
	} // end loop subsets

//	return the summed integral contributions of all elements
	return integral;
}

} // namespace ug

#endif /*__H__LIBDISCRETIZATION__FUNCTION_SPACES__INTEGRATEDRAFT__*/
